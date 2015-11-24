"""
@author Jan-Hendrik Prinz
"""

from test_helpers import raises_with_message_like

import openpathsampling.pathsimulator as simulator

class testAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_volume(self):
        mover = simulator.PathSimulator()
