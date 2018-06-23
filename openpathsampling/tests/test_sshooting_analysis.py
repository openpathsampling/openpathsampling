from .test_helpers import (data_filename)
from nose.tools import (assert_equal, assert_not_equal,
                        assert_true, assert_greater, assert_raises,
                        assert_almost_equal)
import openpathsampling as paths
import openpathsampling.engines.toy as toys
from openpathsampling.pathsimulators.sshooting_simulator import *
from openpathsampling.analysis.sshooting_analysis import *
import numpy as np
import os

class testSShootingAnalysis(object):
    def setup(self):
        # PES is one-dimensional negative slope (y(x) = -x)
        pes = toys.LinearSlope(m=[-1.0], c=[0.0])
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.02)
        options = {
            'integ' : integrator,
            'n_frames_max' : 1000,
            'n_steps_per_frame' : 1
        }
        self.engine = toys.Engine(options=options, topology=topology)
        # test uses snapshots with different velocities
        # (1) does not ever reach A
        # (2) goes form A to S to B
        # (3) goes from B to S to A
        self.initial_snapshots = [toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.1]]),
                                      velocities=np.array([[2.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[-0.2]]),
                                      velocities=np.array([[-2.0]]),
                                      engine=self.engine)]
        # trajectory length is set to 100 steps
        self.l = 100
        # reaction coordinate is just x coordinate
        rc = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        # state A: [-inf, -1]
        self.state_A = paths.CVDefinedVolume(rc, float("-inf"), -1.0)
        # state B: [1, inf]
        self.state_B = paths.CVDefinedVolume(rc, 1.0, float("inf"))
        # state S: [-0.5, 0.5]
        self.state_S = paths.CVDefinedVolume(rc, -0.5, 0.5)
        # define state labels
        self.state_labels = {
            "A" : self.state_A,
            "B" : self.state_B,
            "S" : self.state_S,
            "None" :~(self.state_A | self.state_B | self.state_S)}
        # velocities are not randomized
        randomizer = paths.NoModification()

        self.filename = data_filename("sshooting_test.nc")
        self.storage = paths.Storage(self.filename, mode="w")
        self.storage.save(self.initial_snapshots)

        self.simulation = SShootingSimulation(
                              storage=self.storage,
                              engine=self.engine,
                              state_S=self.state_S,
                              randomizer=randomizer,
                              initial_snapshots=self.initial_snapshots,
                              trajectory_length=self.l)
        self.simulation.output_stream = open(os.devnull, 'w')
        self.simulation.run(n_per_snapshot=1)
        self.analysis = None

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_analysis(self):
        self.storage = paths.Storage(self.filename, mode="r")
        self.analysis = SShootingAnalysis(steps=self.storage.steps,
                                          states=[self.state_A,
                                                  self.state_B,
                                                  self.state_S])

        # check wrong analyze_single_step() argument
        empty_dict = self.analysis.analyze_single_step(3.1415)
        assert_equal(empty_dict, {})

        # dictionary with three entries returned
        assert_equal(len(self.analysis), 3)

        # run C_AB() in advance, this triggers calculate_averages()
        cab = self.analysis.C_AB()
        assert_equal(len(cab), self.l+1)
        assert_greater(cab.sum(), 0)

        # get total and per-snapshot results
        M, Ns, I_Bt, Ns_Bt, hAhB_Bt, sres = self.analysis.calculate_averages()

        # check total number of generated subtrajectories
        assert_equal(M, 3 * (self.l + 1))
        # since we did not use a bias Ns_Bt should be 1.0
        assert_almost_equal(Ns_Bt, 1.0)
        # now do the same tests on a per-snapshot basis
        snaps = sres.keys()
        for snap in snaps:
            assert_equal(sres[snap]["M"], self.l+1)
            assert_almost_equal(sres[snap]["Ns_Bt"], 1.0)
            # check whether time correlation is only filled for case (2), for
            # case (1) and (3) the array should be empty.
            if snap.xyz[0][0] == 0.0:
                assert_almost_equal(np.array(sres[snap]["hAhB_Bt"]).sum(), 0)
            elif snap.xyz[0][0] == 0.1:
                assert_greater(np.array(sres[snap]["hAhB_Bt"]).sum(), 0)
            elif snap.xyz[0][0] == -0.2:
                assert_almost_equal(np.array(sres[snap]["hAhB_Bt"]).sum(), 0)

        # run C_AB() afterwards, this triggers another if clause
        cab = self.analysis.C_AB()
        assert_equal(len(cab), self.l+1)
        assert_greater(cab.sum(), 0)
        
        # try another analysis with different arguments for complete coverage
        def dummy_bias(x): 
            return 1.0 
        cv_b = paths.CoordinateFunctionCV(name="cv_b", f=dummy_bias) 
        rc = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        self.state_S2 = paths.CVDefinedVolume(rc, -0.55, 0.05)
        self.analysis2 = SShootingAnalysis(steps=self.storage.steps,
                                           states=[self.state_A,
                                                   self.state_B,
                                                   self.state_S2],
                                           bias=cv_b)

        assert_raises(NotImplementedError, self.analysis.committor)
        assert_raises(NotImplementedError,
                      self.analysis.committor_histogram)
