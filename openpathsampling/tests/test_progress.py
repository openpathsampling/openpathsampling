import pytest

from openpathsampling.progress import *


def test_silent_progress(capsys):
    out, err = capsys.readouterr()
    x = [1, 2, 3]
    silent_progress(x, desc='foo', leave=True)
    _ = [xx**2 for xx in x]
    out, err = capsys.readouterr()
    assert out == ""
    assert err == ""


class TestTqdmPartial(object):
    def setup(self):
        _ = pytest.importorskip('tqdm')
        self.tqdm_partial = TqdmPartial(desc="foo")

    def test_get_attr(self):
        assert self.tqdm_partial.desc == "foo"
        assert self.tqdm_partial.leave is True
        assert self.tqdm_partial.file is None

    def test_get_kwargs(self):
        assert self.tqdm_partial.kwargs == {'desc': 'foo', 'leave': True,
                                            'file': None}

    def test_set_attr(self):
        self.tqdm_partial.file = "bar"
        assert self.tqdm_partial.kwargs['file'] == "bar"

    def test_call(self, capsys):
        out, err = capsys.readouterr()
        _ = [x for x in self.tqdm_partial([1, 2, 3])]
        out, err = capsys.readouterr()
        assert out == ""
        assert err != ""


class TestSimpleProgress(object):
    def setup(self):
        self.progress = SimpleProgress()

    def test_progress_getter(self):
        assert not hasattr(self.progress, '_progress')
        prog = self.progress.progress
        assert prog is not None
        assert hasattr(self.progress, '_progress')

    def test_progress_setter(self):
        prog = self.progress.progress
        assert prog is not None
        self.progress.progress = 'silent'
        assert self.progress.progress is silent_progress
        if HAS_TQDM:  # starimport from paths.progress
            self.progress.progress = 'tqdm'
            assert isinstance(self.progress.progress, TqdmPartial)
            my_tqdm = TqdmPartial()
            self.progress.progress = my_tqdm
            assert self.progress.progress is my_tqdm

        self.progress.progress = None
        assert self.progress.progress is silent_progress
