import numpy as np


class NN(object):
    def __and__(self, other):
        return self

    def __or__(self, other):
        return self

    def __gt__(self, other):
        return all(l in self for l in other.lengths)

    def __lt__(self, other):
        return all(l in other for l in self.lengths)

    def __eq__(self, other):
        return self > other > self

    def __contains__(self, item):
        return False

    def __nonzero__(self):
        return self > NoNN()

    def __neq__(self):
        return self

    @property
    def lengths(self):
        return []

    @property
    def parts(self):
        return [self]

    def matrix_mult(self, matrix):
        return matrix


class SetNN(NN):
    def __init__(self, length_set):
        super(SetNN, self).__init__()
        self.length_set = length_set

    def __str__(self):
        return 'L{' + ','.join(map(str, self.length_set)) + '}'

    def matrix_mult(self, matrix):
        return np.sum([
            np.linalg.matrix_power(matrix, pw) for pw in self.length_set
        ])

    def __nonzero__(self):
        return bool(self.length_set)

    def __neg__(self):
        mi = min(self.length_set)
        ma = max(self.length_set)

        if mi > 0:
            return MultiNN([
                SliceNN(slice(0, mi)),
                SetNN(set(range(mi, ma + 1)) - self.length_set),
                SliceNN(slice(ma + 1, None))
            ])
        else:
            return MultiNN([
                SetNN(set(range(mi, ma + 1)) - self.length_set),
                SliceNN(slice(ma + 1, None))
            ])

    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        if isinstance(other, SingleNN):
            return NoNN()

        elif isinstance(other, SliceNN):
            if other.length_slice.stop is not None:
                return SetNN(
                    self.length_set &
                    set(range(
                        *other.length_slice.indices(other.length_slice.stop))))
            else:
                return SetNN(
                    self.length_set &
                    set(range(
                        *other.length_slice.indices(max(self.length_set) + 1))))
        elif isinstance(other, SetNN):
            return SetNN(
                self.length_set & other.length_set
            )
        elif isinstance(other, MultiNN):
            return other & self

    def __or__(self, other):
        if other < self:
            return self
        elif other > self:
            return other

        if isinstance(other, SingleNN):
            return SetNN({other.length} | self.length_set)

        elif isinstance(other, SliceNN):
            if other.length_slice.stop is not None:
                return MultiNN([
                    SetNN(
                        self.length_set -
                        set(range(
                            *other.length_slice.indices(
                                other.length_slice.stop)))),
                    other
                ])
            else:
                if min(self.length_set) >= other.length_slice.start:
                    return other
                else:
                    return MultiNN([
                        SetNN(
                            self.length_set -
                            set(range(
                                *other.length_slice.indices(
                                    1 + max(self.length_set))))),
                        other
                    ])
        elif isinstance(other, SetNN):
            return SetNN(
                self.length_set | other.length_set
            )
        elif isinstance(other, MultiNN):
            return other | self

    def __contains__(self, item):
        if isinstance(item, int):
            return item in self.length_set

        elif isinstance(item, slice):
            if item.stop is not None:
                return set(range(*item.indices(item.stop))) <= self.length_set
            else:
                return False
        elif isinstance(item, set):
            return item <= self.length_set

    @property
    def lengths(self):
        return [self.length_set]


class MultiNN(NN):
    def __init__(self, length_list):
        super(MultiNN, self).__init__()
        self._parts = [l for l in length_list if l]

    def __contains__(self, item):
        for l in self.parts:
            if item in l:
                return True

        return False

    def __neg__(self):
        return reduce(lambda x, y: x & y, self._parts)

    @property
    def parts(self):
        return self._parts

    @property
    def lengths(self):
        return sum(map(lambda p: p.lengths, self._parts), [])

    def __nonzero__(self):
        if len(self.parts) > 0:
            return any(self._parts)
        else:
            return False

    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        return MultiNN([
            other & l for l in self.parts
        ])

    def __or__(self, other):
        if other < self:
            return self
        elif other > self:
            return other

        return MultiNN([
            l | self for l in other.length_list
        ])

    def matrix_mult(self, matrix):
        return np.sum([
            l.matrix_mult(matrix) for l in self._parts
        ])

    def __str__(self):
        return 'L[' + ','.join(
            map(lambda x: x.__str__()[1:], self._parts)) + ']'


class SingleNN(NN):
    def __init__(self, length):
        super(SingleNN, self).__init__()
        self.length = length

    def __str__(self):
        return 'L(%d)' % self.length

    def __nonzero__(self):
        return True

    def __neg__(self):
        if self.length > 0:
            return MultiNN([
                SliceNN(slice(0, self.length)),
                SliceNN(slice(self.length + 1, None))
            ])
        else:
            return MultiNN([
                SliceNN(slice(self.length + 1, None))
            ])

    def __contains__(self, item):
        if isinstance(item, int):
            return item == self.length
        elif isinstance(item, slice):
            if self.length == 0:
                if item.start is None or item.start == 0 and item.stop == 1:
                    return True
            else:
                if item.start == self.length and item.stop == self.length + 1:
                    return True

            return False
        elif isinstance(item, set):
            return {self.length} == item

    def __and__(self, other):
        if other > self:
            return self

        # if the single state is not smaller than other it must be empty
        return NoNN()

    def __or__(self, other):
        if other < self:
            return self

        elif other > self:
            return other

        if isinstance(other, SingleNN):
            return SetNN({other.length, self.length})

        elif isinstance(other, SliceNN):
            if self.length < other.length_slice.start or (
                    other.length_slice.stop is not None and
                    self.length > other.length_slice.stop - 1
            ):
                return MultiNN([
                    self,
                    other
                ])
            else:
                return other
        elif isinstance(other, SetNN):
            return SetNN(
                {self.length} | other.length_set
            )
        elif isinstance(other, MultiNN):
            return other | self

    @property
    def lengths(self):
        return [self.length]

    def matrix_mult(self, matrix):
        return np.linalg.matrix_power(matrix, self.length)


class SliceNN(NN):
    def __init__(self, length_slice):
        super(SliceNN, self).__init__()
        self.length_slice = length_slice

    def __str__(self):
        return 'L[%s, ..., %s)' % \
            (str(self.length_slice.start), str(self.length_slice.stop))

    def matrix_mult(self, matrix):
        start = self.length_slice.start
        if start is None:
            start = 0

        stop = self.length_slice.stop
        if stop is None:
            return np.dot(
                np.linalg.matrix_power(matrix, start),
                np.linalg.inv(np.identity(len(matrix)) - matrix)
            )
        else:
            return np.dot(
                np.linalg.matrix_power(matrix, start) -
                np.linalg.matrix_power(matrix, stop),
                np.linalg.inv(np.identity(len(matrix)) - matrix)
            )

    def __neg__(self):
        mi = self.length_slice.start
        ma = self.length_slice.stop

        if mi is not None and mi > 0 and ma is not None:
            return MultiNN([
                SliceNN(slice(0, mi)),
                SliceNN(slice(ma + 1, None))
            ])
        elif mi is not None and mi > 0:
            return MultiNN([
                SliceNN(slice(0, mi))
            ])

        elif ma is not None:
            return MultiNN([
                SliceNN(slice(ma + 1, None))
            ])
        else:
            return NoNN()

    def __nonzero__(self):
        if self.length_slice.stop is None:
            return True
        elif self.length_slice.start is None:
            return self.length_slice.stop is not 0
        else:
            return self.length_slice.start < self.length_slice.stop

    def __contains__(self, item):
        if isinstance(item, int):
            if (
                    self.length_slice.start is None or
                    self.length_slice.start <= item) and \
                (
                    self.length_slice.stop is None or
                    self.length_slice.stop > item):
                return True

            return False
        elif isinstance(item, slice):
            if item.start is None and item.stop is None:
                return self.length_slice.start is None and \
                    self.length_slice.stop is None
            elif item.start is None:
                return self.length_slice.start is None and \
                    self.length_slice.stop >= item.stop
            elif item.stop is None:
                return self.length_slice.stop is None and \
                    self.length_slice.start <= item.start
            else:
                return self.length_slice.start <= item.start and \
                    self.length_slice.stop >= item.stop

        elif isinstance(item, set):
            return self.length_slice.start <= min(item) and \
                self.length_slice.stop > max(item)

    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        if isinstance(other, SingleNN):
            return other & self

        elif isinstance(other, SliceNN):
            if other.length_slice.stop is None:
                stop = self.length_slice.stop
            elif self.length_slice.stop is None:
                stop = other.length_slice.stop
            else:
                stop = min(self.length_slice.stop, other.length_slice.stop)

            if other.length_slice.start is None:
                start = self.length_slice.start
            elif self.length_slice.start is None:
                start = other.length_slice.start
            else:
                start = max(self.length_slice.start, other.length_slice.start)

            if start is not None and stop is not None:
                stop = max(start, stop)

            return SliceNN(slice(start, stop))
        elif isinstance(other, SetNN):
            return other & self
        elif isinstance(other, MultiNN):
            return other & self

    def __or__(self, other):

        if other < self:
            return self
        elif other > self:
            return other

        if isinstance(other, SingleNN):
            return other | self

        elif isinstance(other, SliceNN):
            if other.length_slice.stop is None:
                stop = None
                left = self.length_slice.stop
            elif self.length_slice.stop is None:
                stop = None
                left = other.length_slice.stop
            else:
                stop = max(self.length_slice.stop, other.length_slice.stop)
                left = min(self.length_slice.stop, other.length_slice.stop)

            if other.length_slice.start is None:
                start = None
                right = self.length_slice.start
            elif self.length_slice.start is None:
                start = None
                right = other.length_slice.start
            else:
                start = max(self.length_slice.start, other.length_slice.start)
                right = min(self.length_slice.start, other.length_slice.start)

            if left is not None and right is not None and left < right:
                return MultiNN([
                    self, other
                ])
            else:
                return SliceNN(slice(start, stop))
        elif isinstance(other, SetNN):
            return other | self
        elif isinstance(other, MultiNN):
            return other | self

    @property
    def lengths(self):
        return [self.length_slice]


class NoNN(NN):
    def __init__(self):
        super(NoNN, self).__init__()

    def __nonzero__(self):
        return False

    def __and__(self, other):
        return self

    def __or__(self, other):
        return other

    def __contains__(self, item):
        return False

    @property
    def lengths(self):
        return []

    def __str__(self):
        return 'L{}'

    def __neg__(self):
        return AllNN()


class AllNN(SliceNN):
    def __init__(self):
        super(AllNN, self).__init__(slice(None))

    def __nonzero__(self):
        return True

    def __neg__(self):
        return NoNN()

    def __and__(self, other):
        return other

    def __or__(self, other):
        return self

    def __contains__(self, item):
        return True

    def __str__(self):
        return 'L{...}'


class OneNN(SingleNN):
    def __init__(self):
        super(OneNN, self).__init__(1)


class RangeLength(SliceNN):
    def __init__(self, min_length, max_length):
        super(RangeLength, self).__init__(slice(min_length, max_length))


class UnionLength(MultiNN):
    def __init__(self, lengths):
        super(UnionLength, self).__init__(sum((l.lengths for l in lengths), []))
