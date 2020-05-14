"""
@author Alberto Perez de Alba Ortiz
"""
import openpathsampling as paths

try:
    import mdtraj as md
except ImportError:
    HAS_MDTRAJ = False
else:
    HAS_MDTRAJ = True

try:
    import openpathsampling.engines.openmm as peng_omm
except ImportError:
    HAS_OPENMM = False
else:
    HAS_OPENMM = True

try:
    import plumed
except ImportError:
    HAS_PLUMED = False
else:
    HAS_PLUMED = True
    from openpathsampling.collectivevariables.plumed_wrapper\
         import PLUMEDCV, PLUMEDInterface

import pytest
from .test_helpers import data_filename

import numpy as np
import os
import copy


class TestPLUMED(object):
    def setup(self):
        if not HAS_PLUMED:
            pytest.skip("PLUMED module not installed.")
        if not HAS_OPENMM:
            pytest.skip("OpenMM module not installed.")
        if not HAS_MDTRAJ:
            pytest.skip("MDTraj module not installed.")
        self.topology = peng_omm.tools.topology_from_pdb(data_filename(
                                                         "plumed_wrapper/" +
                                                         "AD_initial_frame" +
                                                         ".pdb"))

        self.trajectory = peng_omm.trajectory_from_mdtraj(md.load(
                                                          data_filename(
                                                          "ala_small_traj" +
                                                          ".pdb")))

        if os.path.isfile("test.nc"):
            os.remove("test.nc")

    def teardown(self):
        if os.path.isfile("test.nc"):
            os.remove("test.nc")
        root = os.listdir(".")
        for item in root:
            if item.endswith("plumed.log"):
                os.remove(item)


class TestPLUMEDCV(TestPLUMED):

    def test_storage(self):
        plmd = PLUMEDInterface(self.topology)
        dist = PLUMEDCV("dist", plmd, "DISTANCE ATOMS=7,9")
        store = paths.Storage("test.nc", "w")
        store.save(self.trajectory[:2])
        store.save(dist)
        store.sync()
        store.close()

        store2 = paths.Storage('test.nc', 'r')
        dist2 = store2.cvs[dist.name]
        store2.close()

        np.testing.assert_almost_equal(dist(self.trajectory),
                                       dist2(self.trajectory), decimal=6)

    def test_storage_components(self):
        plmd = PLUMEDInterface(self.topology)
        dist = PLUMEDCV("dist", plmd, "DISTANCE ATOMS=7,9 COMPONENTS",
                        components=["x", "y", "z"])
        store = paths.Storage("test.nc", "w")
        store.save(self.trajectory[:2])
        store.save(dist)
        store.sync()
        store.close()

        store2 = paths.Storage('test.nc', 'r')
        dist2 = store2.cvs[dist.name]
        store2.close()

        np.testing.assert_almost_equal(dist(self.trajectory),
                                       dist2(self.trajectory), decimal=6)

    def test_distance_vs_mdtraj(self):
        plmd = PLUMEDInterface(self.topology)
        dist_md = paths.MDTrajFunctionCV("dist_md", md.compute_distances,
                                         self.topology, atom_pairs=[[6, 8]])
        dist_pl = PLUMEDCV("dist_pl", plmd, "DISTANCE ATOMS=7,9")

        np.testing.assert_almost_equal(dist_md(self.trajectory),
                                       dist_pl(self.trajectory), decimal=3)

    def test_angle_vs_mdtraj(self):
        plmd = PLUMEDInterface(self.topology)
        ang_md = paths.MDTrajFunctionCV("ang_md", md.compute_angles,
                                        self.topology,
                                        angle_indices=[[4, 6, 8]])
        ang_pl = PLUMEDCV("ang_pl", plmd, "ANGLE ATOMS=5,7,9")

        np.testing.assert_almost_equal(ang_md(self.trajectory),
                                       ang_pl(self.trajectory), decimal=3)

    def test_torsion_vs_mdtraj(self):
        plmd = PLUMEDInterface(self.topology)
        tor_md = paths.MDTrajFunctionCV("tor_md", md.compute_dihedrals,
                                        self.topology,
                                        indices=[[6, 8, 14, 16]])
        tor_pl = PLUMEDCV("tor_pl", plmd, "TORSION ATOMS=7,9,15,17")

        np.testing.assert_almost_equal(tor_md(self.trajectory),
                                       tor_pl(self.trajectory), decimal=3)

    def test_rmsd_vs_mdtraj(self):
        plmd = PLUMEDInterface(self.topology)
        rmsd_ref_file = data_filename(os.path.join("plumed_wrapper",
                                                   "AD_plumed_rmsd.pdb"))
        md_ref = md.load(rmsd_ref_file)
        rmsd_md = paths.MDTrajFunctionCV("rmsd_md", md.rmsd, self.topology,
                                         reference=md_ref, frame=0,
                                         atom_indices=range(22))
        rmsd_pl = PLUMEDCV("rmsd_pl", plmd, "RMSD REFERENCE="
                           + rmsd_ref_file + " TYPE=OPTIMAL")

        np.testing.assert_almost_equal(rmsd_md(self.trajectory),
                                       rmsd_pl(self.trajectory), decimal=3)

    def test_components(self):
        plmd = PLUMEDInterface(self.topology)
        dist_pl = PLUMEDCV("dist", plmd, "DISTANCE ATOMS=7,9")
        comp_pl = PLUMEDCV("comp", plmd, "DISTANCE ATOMS=7,9 COMPONENTS",
                           components=["x", "y", "z"])
        comp_dist_pl = (comp_pl(self.trajectory)[:, 0]**2 +
                        comp_pl(self.trajectory)[:, 1]**2 +
                        comp_pl(self.trajectory)[:, 2]**2)**0.5
        np.testing.assert_almost_equal(dist_pl(self.trajectory),
                                       comp_dist_pl, decimal=6)

    def test_combine(self):
        plmd = PLUMEDInterface(self.topology)
        phi_pl = PLUMEDCV("phi", plmd, "TORSION ATOMS=5,7,9,15")
        psi_pl = PLUMEDCV("psi", plmd, "TORSION ATOMS=7,9,15,17")
        sum_pl = phi_pl(self.trajectory) + psi_pl(self.trajectory)
        comb_pl = PLUMEDCV("sum", plmd, "COMBINE ARG=phi,psi PERIODIC=NO")
        np.testing.assert_almost_equal(sum_pl,
                                       comb_pl(self.trajectory), decimal=6)

    def test_group(self):
        plmd = PLUMEDInterface(self.topology)
        tor_pl = PLUMEDCV("tor", plmd, "TORSION ATOMS=7,9,15,17")
        plmd.set("group", "GROUP ATOMS=7,9,15,17")
        tor_group_pl = PLUMEDCV("tor_group", plmd, "TORSION ATOMS=group")

        np.testing.assert_almost_equal(tor_pl(self.trajectory),
                                       tor_group_pl(self.trajectory),
                                       decimal=6)

    def test_molinfo(self):
        plmd = PLUMEDInterface(self.topology,
                               molinfo=data_filename("plumed_wrapper/" +
                                                     "AD_initial_frame.pdb"))
        tor_md = paths.MDTrajFunctionCV("tor_md", md.compute_dihedrals,
                                        self.topology,
                                        indices=[[6, 8, 14, 16]])
        tor_pl = PLUMEDCV("tor_pl", plmd, "TORSION ATOMS=@psi-2")

        np.testing.assert_almost_equal(tor_md(self.trajectory),
                                       tor_pl(self.trajectory), decimal=3)

    '''
    def test_no_masses(self):
        plmd = PLUMEDInterface(self.topology)
        plmd.set("center","CENTER ATOMS=7,9,15,17")
        plmd.set("com","CENTER ATOMS=7,9,15,17 MASS")
        dist_pl = PLUMEDCV("dist",plmd,"DISTANCE ATOMS=center,com")

        np.testing.assert_array_equal(dist_pl(self.trajectory),
                                      np.zeros(self.trajectory))
    '''


class TestPLUMEDInterface(TestPLUMED):

    def test_storage(self):
        plmd = PLUMEDInterface(self.topology)
        plmd.set("", "UNITS LENGTH=nm")
        plmd.set("", "MOLINFO STRUCTURE=" +
                 data_filename("plumed_wrapper/AD_initial_frame.pdb"))
        store = paths.Storage("test.nc", "w")
        store.save(self.trajectory[:2])
        store.tags['a'] = plmd
        store.sync()
        store.close()

        store2 = paths.Storage('test.nc', 'r')
        plmd2 = store2.tags['a']
        store2.close()

        np.testing.assert_array_equal(plmd.get(), plmd2.get())
