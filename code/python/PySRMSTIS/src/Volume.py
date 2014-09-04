'''
Created on 03.09.2014

@author: jan-hendrikprinz
'''

class Volume(object):
    def __init__(self):
        pass
    
    def inside(self, snapshot):
        '''
        Returns 1 if the given snapshot is part of the defined Region
        '''
        
        return False
    
    def outside(self, snapshot):
        '''
        Returns 1 if the given snapshot is not part of the defined Region
        '''
        return not self.inside(snapshot)
        
    def surface(self, snapshot):
        '''
        Returns a surface function of the given
        '''
        pass
    
    def __str__(self):
        return 'Volume'

    def __or__(self, other):
        if self is other:
            return self
        else:
            return VolumeCombination(self, other, fnc = lambda a,b : a or b, str_fnc = '{0} or {1}')

    def __xor__(self, other):
        if self is other:
            # A xor A is always Empty
            return EmptyVolume()
        else:
            return VolumeCombination(self, other, fnc = lambda a,b : a ^ b, str_fnc = '{0} xor {1}')

    def __and__(self, other):
        if self is other:
            return self
        else:
            return VolumeCombination(self, other, fnc = lambda a,b : a and b, str_fnc = '{0} and {1}')

    def __sub__(self, other):
        if self is other:
            # A \ A is always Empty
            return EmptyVolume()
        else:
            return VolumeCombination(self, other, fnc = lambda a,b : a and not b, str_fnc = '{0} and not {1}')

    
class VolumeCombination(Volume):
    def __init__(self, volume1, volume2, fnc, str_fnc):
        super(VolumeCombination, self).__init__()
        self.volume1 = volume1
        self.volume2 = volume2
        self.fnc = fnc
        self.sfnc = str_fnc

    def inside(self, snapshot):
        return self.fnc(self.volume1.inside(snapshot), self.volume2.inside(snapshot))
    
    def __str__(self):
        return '(' + self.sfnc.format(str(self.volume1), str(self.volume2)) + ')'
    
class EmptyVolume(Volume):
    def __init__(self):
        super(EmptyVolume, self).__init__()

    def inside(self, snapshot):
        return False
    
    def __str__(self):
        return 'empty'
    
class FullVolume(Volume):
    def __init__(self):
        super(EmptyVolume, self).__init__()

    def inside(self, snapshot):
        return True
    
    def __str__(self):
        return 'all'

class LambdaVolume(Volume):
    def __init__(self, orderparameter, lambda_min = 0.0, lambda_max = 1.0):
        super(LambdaVolume, self).__init__()
        self.orderparameter = orderparameter
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        
    def inside(self, snapshot):
        l = self.orderparameter(snapshot)
        return l >= self.lambda_min and l <= self.lambda_max

    def __str__(self):
        return '{{x|or(x) in [{0}, {1}]}}'.format( self.lambda_min, self.lambda_max)
    
class Surface(object):
    def __init__(self, volume1, volume2, fnc):
        super(VolumeCombination, self).__init__()
        pass
    
    def cross(self, left, right):
        '''
        Returns 1 if the given snapshot is part of the defined Region
        '''
        return False