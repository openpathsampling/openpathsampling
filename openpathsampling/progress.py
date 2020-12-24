"""
Wrappers for progress bars. Mainly intended to use tqdm if it is available
on the user's system, but in principle, other progress bars could be used as
well.
"""

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover
    HAS_TQDM = False
else:  # pragma: no cover
    HAS_TQDM = True


def silent_progress(iterable, *args, **kwargs):
    """No progress output: ignores all options, and just gives the iterator

    Parameters
    ----------
    iterable : iterable
        thing to iterate over
    """
    return iter(iterable)


def DefaultProgress(**kwargs):  # pragma: no cover
    """Factory for getting the default progress behavior.

    Currently, uses tqdm if it is available; silent progress otherwise. In
    the future, we may add global configuration (e.g., rcparams) to allow
    users to customize this behavior.
    """
    # TODO: eventually, add some sort of rcparams support for this
    if HAS_TQDM:
        return TqdmPartial(**kwargs)
    else:
        return silent_progress


class TqdmPartial(object):
    """Progress bar wrapper for tqdm-based progress meters.

    Note that additional tqdm parameters can be set *after* initialization.

    Parameters
    ----------
    desc : str
        label for the progress bar
    leave : bool
        whether to leave to progress bar visible after completion (default
        True)
    """
    def __init__(self, desc=None, leave=True, file=None):
        self.kwargs = {'desc': desc,
                       'leave': leave,
                       'file': file}

    def __getattr__(self, attr):
        if attr != 'kwargs':
            return self.kwargs[attr]
        else:
            return super(TqdmPartial, self).__getattr__(attr)

    def __setattr__(self, attr, value):
        if attr != 'kwargs':
            self.kwargs[attr] = value
        else:  # pragma: no cover
            super(TqdmPartial, self).__setattr__(attr, value)

    def __call__(self, iterable, **kwargs):
        tqdm_kwargs = self.kwargs
        tqdm_kwargs = {k: v for k, v in self.kwargs.items()}
        tqdm_kwargs.update(kwargs)
        return tqdm(iterable, **tqdm_kwargs)


class SimpleProgress(object):
    """Mix-in for classes that need a progress meter.

    Note that this will use the names ``progress`` and ``_progress``.
    """
    @property
    def progress(self):
        if not hasattr(self, '_progress'):
            self._progress = DefaultProgress()
        return self._progress

    @progress.setter
    def progress(self, value):
        if value == 'tqdm':
            value = TqdmPartial() if HAS_TQDM else silent_progress
        elif value in ['silent', None]:
            value = silent_progress
        # else we assume it's already an AnalysisProgress
        self._progress = value
