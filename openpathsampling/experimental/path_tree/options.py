def canonicalize_mover(mover):
    if mover is not None:
        return mover.__class__.__name__


class OptionDict:
    def __init__(self, dictionary, default):
        self.dictionary = dictionary
        self.default = default

    def __getitem__(self, item):
        try:
            return self.dictionary[item]
        except KeyError:
            return self.default

    def __getattr__(self, attr):
        if attr in ['dictionary', 'default']:
            return object.__getattribute__(self, attr)
        return self[attr]

    def __setattr__(self, key, value):
        if key not in ['dictionary', 'default']:
            self.dictionary[key] = value
        else:
            object.__setattr__(self, key, value)

    def __repr__(self):
        return (f"OptionDict({repr(self.dictionary)}, "
                f"default={repr(self.default)}")


class PathTreeMoverOptions(object):
    def __init__(self, show='full', color='#0f0'):
        self.show = show
        self.color = color

        self.get_left_right_methods = {
            'full': self.get_left_right_full,
            'new': self.get_left_right_new,
        }

    @staticmethod
    def get_left_right_full(step):
        offset = step.offset if step.offset is not None else 0
        return [(offset, offset + len(step.trajectory))]

    @staticmethod
    def get_left_right_new(step):
        offset = step.offset if step.offset is not None else 0
        return [(offset + start, offset + stop)
                for start, stop in step.segments.new]

    def get_left_right(self, step):
        get_left_right = self.get_left_right_methods[self.show]
        return get_left_right(step)

    def __repr__(self):
        return (f"{self.__class__.__name__}(show='{self.show}', "
                f"color='{self.color})")


class PathTreeOptions(object):
    def __init__(self, movers, rejected, allowed_movers, forbidden_movers):
        self.movers = movers
        self.rejected = rejected
        self.allowed_movers = allowed_movers
        self.forbidden_movers = forbidden_movers

    def _skip_step(self, step):
        mover = canonicalize_mover(step.mover)
        include_rejected = self.rejected != 'hidden'
        skip_rejected = not include_rejected and not step.accepted
        skip_allowed = self.allowed_movers and mover not in allowed
        skip_forbidden = self.forbidden_movers and mover in forbidden
        return skip_rejected or skip_allowed or skip_forbidden

    def filter_tree_steps(self, tree_steps):
        for step in tree_steps:
            if not self._skip_step(step):
                yield step


def create_default_options():
    two_way_options = PathTreeMoverOptions(color='#9A63EE')
    movers = OptionDict(
        {
            'ForwardShootMover': PathTreeMoverOptions(show='new',
                                                      color='#00f'),
            'BackwardShootMover': PathTreeMoverOptions(show='new',
                                                       color='#f00'),
            'BackwardFirstTwoWayShootingMover': two_way_options,
            'ForwardFirstTwoWayShootingMover': two_way_options,
            None: PathTreeMoverOptions(color='#808080'),

        },
        default=PathTreeMoverOptions(color="#377E22")
    )
    return PathTreeOptions(
        movers=movers,
        rejected='hidden',
        allowed_movers=None,
        forbidden_movers=None,
    )
