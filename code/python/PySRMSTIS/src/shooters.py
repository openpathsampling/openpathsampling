

class ShootingPointSelector(object):
    def f(self, frame):
        return 1.0

    def select(self, path):
        totbias = 0.0
        for snapshot in path:
            totbias += self.f(snapshot)
        pass


class GaussianBias(ShootingPointSelector):
    def f(self, frame):
        return exp(-self.alpha*(frame.lambda-self.l0)**2)
            


class UniformSelector(ShootingPointSelector):
    def select(self, path):
        return (chosen_snapshot, 1.0)