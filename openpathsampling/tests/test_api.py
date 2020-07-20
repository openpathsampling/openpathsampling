import pytest
import importlib
import os.path

from .test_helpers import data_filename

# select the file you want to use, or change to the file for your own API
# (add more files to ensure that we support 1.0, 1.1, etc.)
API_FILES = ["ops1x0_api.txt"]
COMBO_ERROR = True

@pytest.mark.parametrize('api_file', API_FILES)
def test_api(api_file):
    def import_object(import_path):
        path = import_path.split('.')
        current = importlib.import_module(path[0])
        for mod in path[1:]:
            current = getattr(current, mod)
        return current

    api_file = os.path.join("apis", api_file)
    with open(data_filename(api_file), mode='r') as api_file:
        tests = api_file.read().splitlines()

    errors = []
    for test in tests:
        try:
            obj = import_object(test)
        except AttributeError as e:
            if COMBO_ERROR:
                errors.append(test)
            else:
                raise e

    if errors:
        raise AssertionError("Some API entries failed: " + str(errors))

