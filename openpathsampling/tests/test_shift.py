from openpathsampling.shifting import *
import openpathsampling as paths
import openpathsampling.tests.test_helpers as ops_test
from openpathsampling.tests.test_helpers import make_1d_traj
from openpathsampling.tests.test_helpers import CalvinistDynamics

from nose.plugins.skip import Skip, SkipTest

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_true, assert_in, assert_not_in)

from numpy.testing import assert_allclose

from openpathsampling.collectivevariable import FunctionCV
from openpathsampling.engines.trajectory import Trajectory
#from openpathsampling.ensemble import EnsembleFactory as ef
from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.pathmover import *
#from openpathsampling.pathmover import IdentityPathMover
from openpathsampling.sample import Sample, SampleSet
from openpathsampling.shooting import UniformSelector
from openpathsampling.volume import CVDefinedVolume
import openpathsampling.engines.toy as toys

from openpathsampling.high_level.move_scheme import MoveScheme, DefaultScheme
from openpathsampling.high_level.move_strategy import *
#from openpathsampling import VolumeFactory as vf

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


import collections


from openpathsampling.tests.test_helpers import (
    true_func, assert_equal_array_array, MoverWithSignature,
    setify_ensemble_signature, reorder_ensemble_signature
)


class TestShiftingMover(object):
	def setup(self):

		
		ensemble = paths.LengthEnsemble(5)
		traj = ops_test.make_1d_traj([0, 1, 2, 3, 4])
		self.traj=traj

		input_sample = paths.Sample(trajectory=traj,
			ensemble=ensemble, replica=0)
	
		engine = ops_test.CalvinistDynamics(
			predestination=[ -2, -1, 0, 4, 5, 6]
			)
		self.engine = engine

		self.pathmover = ShiftingMover(ensemble, 2, self.engine)

	
	def test_run_set(self):
	#test the length of the set of unique trajectory values is exactly 2, corresponding to the forward
	#and backward shifted trajectory

		traj_list=[]
		detail_list=[]		
		for i in range(100):
			x, d = self.pathmover._run(self.traj, 0)
			traj_list.append(x)
			detail_list.append(d)


		traj_set=set(traj_list)
		
		#y = collections.Counter(traj_list)
		assert_equal(len(traj_set), 2) #test the number of unique trajectory values is 2
		#print traj_list 
		assert_not_equal(len(traj_set), 1)
		

	
	def test_run_traj(self):
	#test if exact trajectory coordinates for a backward and forward shift are present in the generated
	#trajectories (by the predetermined Calvinist Dynamics Engine)

		traj_fwd  = ops_test.make_1d_traj([2,3,4,5,6])
		traj_bkwd = ops_test.make_1d_traj([-2,-1,0,1,2])

		traj_list = []
	
		for i in range(100):
			x, d = self.pathmover._run(self.traj, 0)
			traj_list.append(x.xyz.tostring())
		traj_set = set(traj_list)

		assert_in(traj_fwd.xyz.tostring(), traj_set)
		assert_in(traj_bkwd.xyz.tostring(), traj_set)
		

		

	def test_run_len_traj(self):
	#check the length of every trajectory in the generated tractory list is equal to 5
	#also added in test to check it matches the length of initial traj and the anticipated
	#fwd and bkwd traj all of which are length 5, this may be redundent but means test will work
	#if length of ensemble is changed.
		
		traj_fwd  = ops_test.make_1d_traj([2,3,4,5,6])
		traj_bkwd = ops_test.make_1d_traj([-2,-1,0,1,2])

		traj_list = []
	
		for i in range(50):
			x, d = self.pathmover._run(self.traj, 0)
			traj_list.append(x)
		#traj_set=set(traj_list)

		for j in traj_list:

			assert_equal(len(traj_fwd.xyz), len(j.xyz))
			assert_equal(len(traj_bkwd.xyz),len(j.xyz))
			assert_equal(len(self.traj.xyz),len(j.xyz))
			assert_equal(len(j.xyz),5)
		

		#print(len(traj_list[1].xyz))	

	def test_run_details(self):
	#check keys and values of the dictonary output from the _run function of the Mover
	#keys are the name of the move which has 1 choice - "Shift" 
	#values are the shift which has 2 choices, either +2 or -2 depending if fwd or bkwd shift.  
	#also test the value of the shift length is +/-2 in all cases
	
		
		value_list = []
		key_list = []		
		for i in range(100):
			x, d = self.pathmover._run(self.traj, 0)
			for keys,values in d.iteritems():
				key_list.append(keys)
				value_list.append(values)
					
		
		key_set = set(key_list)
		value_set = set(value_list)
		
		assert_equal(len(value_set), 2) #test the number of unique values is 2, either+/-2
		assert_equal(len(key_set), 1)  #test the number of unique keys is 1
		
		#assert_equal(str(key_list), str('Shift'))
		
		
		#for j in value_list:
		#	assert_in(j, 2)
		#	assert_in(j, -2)
		#	assert_not_in(j, 3)

		        #print(type(j))

		shift_list = []
		for k in key_list:
			if k in "Shift":
				shift_list.append(k)
				assert_equal(str(k),"Shift")

		assert_equal(len(shift_list), 100)

		
		fwd_list = []
		bwd_list = []

		for l in value_list:
			if l > 0:
				fwd_list.append(l)
				assert_equal(l, 2)
				assert_not_equal(l,1)
			elif l <0:
				bwd_list.append(l)
				assert_equal(l, -2)
				assert_not_equal(l, -1)

		assert_equal(len(fwd_list) + len(bwd_list), 100)

		#for x in fwd_list:
		#	print(x)

		#for y in bwd_list:
		#	print(y)

		#print(len(fwd_list)
		#print(len(bwd_list)

		
#Test the Shifting Strategy by creating Test Setup to include state definitions, interfaces and
# a network containing the ensembles of paths.
class ShiftStrategyTestSetup(object):
    def setup(self):
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
	self.stateA = paths.CVDefinedVolume(cvA, float("-inf"), -5)
        self.stateB = paths.CVDefinedVolume(cvB, float("-inf"), -5)
        interfacesA = paths.VolumeInterfaceSet(cvA, float("-inf"), 
                                               [ -3, -1])
        interfacesB = paths.VolumeInterfaceSet(cvB, float("-inf"), 
                                               [-3, -1])
        self.network = paths.MSTISNetwork(
            [(self.stateA, interfacesA),
             (self.stateB, interfacesB)],
              ms_outers=paths.MSOuterTISInterface.from_lambdas(
                {interfacesA: 0.0, interfacesB: 0.0}
		)
             )


class TestShiftingStrategy(ShiftStrategyTestSetup):
    def test_make_movers(self):
        strategy = ShiftingStrategy()
	self.strategy=strategy
        scheme = MoveScheme(self.network) 
        movers = strategy.make_movers(scheme)
        
	assert_equal(len(movers), 4) #setup contains 4 ensembles, test to check strategy applying to all of them
	#print(movers)

        for mover in movers:
           assert_equal(str(mover), str('Shifting')) #check Shifting Mover applied using Shifting Strategy
	
	
	

