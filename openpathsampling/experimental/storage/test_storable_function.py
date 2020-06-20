import pytest
try:
    from unittest import mock
except ImportError:
    import mock

import numpy as np

from openpathsampling.tests.test_helpers import make_1d_traj

from .serialization_helpers import get_uuid
from .storable_functions import *
_MODULE = "openpathsampling.experimental.storage.storable_functions"

def test_requires_lists_pre():
    assert requires_lists_pre([1]) == [[1]]
    assert requires_lists_pre([1,2]) == [[1,2]]

@pytest.mark.parametrize('array_input,expected', [
    ([[1], [2], [3]], [1, 2, 3]),
    ([1, 2, 3], [1, 2, 3]),
    ([[1, 2], [3, 4]], [[1, 2], [3, 4]]),
    ([[[1, 2]], [[3, 4]]], [[1, 2], [3, 4]]),
])
def test_scalarize_singletons(array_input, expected):
    np.testing.assert_array_equal(
        scalarize_singletons(np.array(array_input)),
        np.array(expected)
    )

def test_wrap_numpy():
    for inp in [1, [1, 2]]:
        assert isinstance(wrap_numpy(inp), np.ndarray)

class TestStorableFunctionConfig(object):
    def setup(self):
        self.config = StorableFunctionConfig(processors=[
            scalarize_singletons,
            wrap_numpy,
            requires_lists_pre,
            requires_lists_post
        ])

    @staticmethod
    def func(values):
        return np.array([s.xyz[:,0] for s in values])

    def test_register(self):
        assert len(self.config.processors) == 4
        names = ['scalarize_singletons', 'requires_lists_pre',
                 'requires_lists_post', 'wrap_numpy']
        for key in names:
            assert key in self.config.processor_dict
        assert self.config.item_preprocessors == []
        assert self.config.list_preprocessors == [requires_lists_pre]
        assert self.config.item_postprocessors == [scalarize_singletons]
        assert self.config.list_postprocessors == [wrap_numpy,
                                                   requires_lists_post]

        mock_wrap_numpy = Processor(name='wrap_numpy',
                                    stage='item-pre',
                                    func=lambda x: x)
        self.config.register(mock_wrap_numpy)
        proc_dict_wrap_numpy = self.config.processor_dict['wrap_numpy']
        assert len(self.config.processors) == 4
        assert proc_dict_wrap_numpy is mock_wrap_numpy
        assert proc_dict_wrap_numpy is not wrap_numpy
        assert mock_wrap_numpy in self.config.processors
        assert wrap_numpy not in self.config.processors
        assert self.config.item_preprocessors == [mock_wrap_numpy]
        assert self.config.list_preprocessors == [requires_lists_pre]
        assert self.config.item_postprocessors == [scalarize_singletons]
        assert self.config.list_postprocessors == [requires_lists_post]

    @pytest.mark.parametrize('style', ['obj', 'name'])
    def test_deregister(self, style):
        dereg = {'obj': wrap_numpy, 'name': 'wrap_numpy'}[style]
        assert len(self.config.processors) == 4
        self.config.deregister(dereg)
        assert len(self.config.processors) == 3
        assert 'wrap_numpy' not in self.config.processor_dict
        assert wrap_numpy not in self.config.processors

    def test_deregister_error(self):
        with pytest.raises(KeyError):
            self.config.deregister('foo')

    def test_deregister_no_error(self):
        # just run it to ensure it doesn't error out
        self.config.deregister('foo', error_if_missing=False)

    def test_func(self):
        # test of the internally used test func
        snap = make_1d_traj([5.0])[0]
        assert self.func([snap]) == [[5]]

    def test_list_preprocess(self):
        snap = make_1d_traj([5.0])[0]
        assert self.config.list_preprocess([snap]) == [[snap]]

    def test_item_preprocess(self):
        snap = make_1d_traj([5.0])[0]
        assert self.config.item_preprocess(snap) == snap

    def test_item_postprocess(self):
        np.testing.assert_array_equal(
            self.config.item_postprocess(np.array([[5.0]])),
            np.array([5.0])
        )

    def test_list_postprocess(self):
        snap = make_1d_traj([5.0])[0]
        values = self.func([snap])
        np.testing.assert_array_equal(self.config.list_postprocess(values),
                                      np.array([5.0]))

    def test_storable_function_integration(self):
        snap = make_1d_traj([5.0])[0]
        sf = StorableFunction(self.func, result_type='float',
                              func_config=self.config)
        assert sf(snap) == 5.0
        np.testing.assert_array_equal(sf([snap]), np.array([5.0]))


class TestStorableFunctionResults(object):
    def setup(self):
        self.cv = StorableFunction(lambda x: x)
        self.cv.__uuid__ = "funcUUID"
        self.mapping = {'UUID1': "foo",
                        'UUID2': "bar"}
        self.sfr = StorableFunctionResults(self.cv, "funcUUID")
        self.sfr.result_dict = self.mapping
        self.sfr.local_uuids = set(self.mapping.keys())

    def test_get_results_as_dict_cached(self):
        result, missing = self.sfr.get_results_as_dict({'UUID1': "object"})
        assert result == {'UUID1': "foo"}
        assert missing == {}

    def test_get_results_as_dict_missing(self):
        result, missing = self.sfr.get_results_as_dict({"UUID3": "object"})
        assert result == {}
        assert missing == {"UUID3": "object"}

    def test_get_results_as_dict_storage(self):
        pytest.skip()
        pass

    def test_update(self):
        new_sfr = StorableFunctionResults(self.cv, "funcUUID")
        new_sfr.result_dict = {'UUID3': "baz"}
        new_sfr.local_uuids = set(['UUID3'])
        self.sfr.update(new_sfr)
        assert len(self.sfr) == 3
        assert "UUID3" in self.sfr.local_uuids
        assert self.sfr.result_dict["UUID3"] == "baz"

    # TODO: test_cache_results_nonpure_function
    # if you try to cache results that don't match the original, you get an
    # error

    def test_cache_results(self):
        self.sfr.cache_results({"UUID3": "baz"})
        assert len(self.sfr) == 3
        assert "UUID3" in self.sfr.local_uuids
        assert self.sfr.result_dict["UUID3"] == "baz"

    def test_clear(self):
        assert len(self.sfr) != 0
        self.sfr.clear()
        assert len(self.sfr) == 0
        assert self.sfr.result_dict == {}
        assert self.sfr.local_uuids == set([])

    def test_len(self):
        assert len(self.sfr) == 2

    def test_to_dict_from_dict_cycle(self):
        pass


@mock.patch(_MODULE + '.get_uuid', lambda x: x)
@mock.patch(_MODULE + '.has_uuid', lambda x: isinstance(x, str))
class TestStorableFunction(object):
    def setup(self):
        def get_expected(uuid):
            expected = {'uuid': 'eval', 'uuid1': 'other'}
            return expected[uuid]

        self.func = StorableFunction(get_expected)

    def test_gets_source(self):
        pytest.skip()
        pass

    def test_no_source_warning(self):
        pytest.skip()
        pass

    def test_disk_cache_property(self):
        pytest.skip()
        pass

    @pytest.mark.parametrize('mode', ['no-caching', 'analysis',
                                      'production'])
    def test_mode(self, mode):
        self.func.mode = mode
        assert self.func.mode == mode
        if mode == 'no-caching':
            assert self.func.local_cache is None
        else:
            assert self.func.local_cache is not None

    def test_bad_mode(self):
        with pytest.raises(ValueError):
            self.func.mode = 'foo'

    @staticmethod
    def _set_cache(func, mode, found_in, expected):
        if found_in == 'cache':
            func.local_cache.cache_results(expected)
        elif mode == 'no-caching':
            pass
        else:
            func.local_cache.clear()
        pass

    @staticmethod
    def _set_storage(func, mode, found_in, expected):
        if found_in == 'storage':
            def get_storage(cv_uuid, uuids):
                missing = {uuid: uuids[uuid] for uuid in uuids
                           if uuid not in expected.keys()}
                found = {uuid: uuids[uuid] for uuid in uuids
                         if uuid in expected.keys()}
                return {uuid: expected[uuid] for uuid in found}, missing
        else:
            def get_storage(cv_uuid, uuids):
                return {}, dict(uuids)

        storage = mock.MagicMock(get_function_results=get_storage)
        func._handler = storage

    @pytest.mark.parametrize('mode, found_in', [
        ('analysis', 'storage'), ('analysis', 'cache'),
        ('analysis', 'eval'), ('production', 'cache'),
        ('production', 'eval'), ('no-caching', 'eval')
    ])
    def test_call(self, mode, found_in):
        # mode = 'analysis'
        # found_in = 'cache'
        # setup, depending on the parametrized parameters
        expected = {'uuid': 'eval'}
        get_expected = lambda x: expected[x]
        func = StorableFunction(get_expected)

        func.mode = mode
        self._set_cache(func, mode, found_in, expected={'uuid': 'cache'})
        self._set_storage(func, mode, found_in, expected={'uuid': 'storage'})

        # validation of correct behavior
        # NOTE: some of this testing is based on internal behavior, which
        # perhaps shouldn't be in the public-facing API
        if found_in != 'cache' and mode != 'no-caching':
            assert 'uuid' not in func.local_cache.result_dict

        assert func('uuid') == found_in
        if mode != 'no-caching':
            assert func.local_cache.result_dict['uuid'] == found_in

    @pytest.mark.parametrize("found_in_1, found_in_2", [
        ('storage', 'storage'), ('cache', 'cache'), ('eval', 'eval'),
        ('cache', 'eval')
    ])
    def test_call_multiple(self, found_in_1, found_in_2):
        # only test this in analysis
        expected_dict = {'uuid': found_in_1, 'other': found_in_2}
        expected = {
            level: {uuid: expected
                    for uuid, expected in expected_dict.items()
                    if expected == level}
            for level in ['eval', 'cache', 'storage']
        }
        get_expected = lambda x: expected['eval'][x]
        func = StorableFunction(get_expected)
        self._set_cache(func, 'analysis', 'cache',
                        expected=expected['cache'])
        self._set_storage(func, 'analysis', 'storage',
                          expected=expected['storage'])

        assert func(['uuid', 'other']) == [found_in_1, found_in_2]

    def test_to_dict_from_dict_cycle(self):
        pytest.skip()
        pass

    def test_full_serialization_cycle(self):
        pass

    @pytest.mark.parametrize('found_in', ['cache', 'storage', 'eval'])
    def test_analysis_mode_integration(self, found_in):
        pytest.skip()
        pass


class TestStorageFunctionHandler(object):
    def setup(self):
        pass

    def test_register(self):
        pytest.skip()
        pass

    def test_update(self):
        pytest.skip()
        pass

    def test_update_mock_parallel(self):
        pytest.skip()
        # a test to show how this should work if multiple results with same
        # CV come in
        pass

