import pytest
from nose.tools import assert_equal
from openpathsampling.tools import *

import tempfile
import hashlib
import os
from .test_helpers import data_filename

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

def test_pretty_print_seconds():
    # test the basics with full reporting
    assert_equal(pretty_print_seconds(93784),
                 "1 day 2 hours 3 minutes 4 seconds")
    assert_equal(pretty_print_seconds(179990),
                 "2 days 1 hour 59 minutes 50 seconds")
    assert_equal(pretty_print_seconds(31),
                 "31 seconds")

    # test positive n_labels (including testing rounding)
    assert_equal(pretty_print_seconds(93784, n_labels=2),
                 "1 day 2 hours")
    assert_equal(pretty_print_seconds(179990, n_labels=2),
                 "2 days 2 hours")
    assert_equal(pretty_print_seconds(179990, n_labels=3),
                 "2 days 2 hours 0 minutes")
    assert_equal(pretty_print_seconds(31, n_labels=2),
                 "31 seconds")

    # test negative n_labels
    assert_equal(pretty_print_seconds(93784, n_labels=-1),
                 "1.09 days")
    assert_equal(pretty_print_seconds(93784, n_labels=-2),
                 "1 day 2.05 hours")
    assert_equal(pretty_print_seconds(31, n_labels=-2),
                 "31 seconds")

    # test separators
    assert_equal(pretty_print_seconds(93784, separator=", "),
                 "1 day, 2 hours, 3 minutes, 4 seconds")

def test_progress_string():
    assert_equal(progress_string(0, 100, 10),
                 "Starting simulation...\nWorking on first step\n")
    assert_equal(progress_string(1, 11, 9378.4),
                 "Running for 2 hours 36 minutes 18 seconds - "
                 + "9378.40 seconds per step\n"
                 + "Estimated time remaining: 1 day 2.05 hours\n")


def test_ensure_file_dne():
    # when the file doesn't exist and you provide contents, ensure_file
    # should create the missing file
    try:
        tmp_dir = tempfile.TemporaryDirectory()
    except AttributeError:
        # Py2: we'll just skip this test (and not worry when Py2 goes away)
        pytest.skip("Test approach only valid in Python 3")
    filename = os.path.join(tmp_dir.name, "foo.data")
    assert not os.path.exists(filename)
    old_contents = "foo bar baz qux"
    old_hash = hashlib.sha1(old_contents.encode('utf-8')).digest()
    contents, hashed = ensure_file(filename, old_contents, old_hash)
    assert os.path.exists(filename)
    assert contents == old_contents
    assert hashed == old_hash


def test_ensure_file_exists():
    # this is the case where everything works: file exists; hash is correct
    filename = data_filename("ala_small_traj.pdb")
    assert os.path.exists(filename)
    with open(filename, mode='r') as f:
        old_contents = f.read()
    old_hash = hashlib.sha1(old_contents.encode('utf-8')).digest()

    contents, hashed = ensure_file(filename, old_contents, old_hash)
    assert contents == old_contents
    assert hashed == old_hash


def test_ensure_file_hash_mismatch():
    # this is the case where the file exists, but the hash is wrong
    filename = data_filename("ala_small_traj.pdb")
    assert os.path.exists(filename)
    with open(filename, mode='r') as f:
        old_contents = f.read()

    old_hash = "foo"

    with pytest.raises(RuntimeError):
        ensure_file(filename, old_contents, old_hash)


def test_ensure_file_no_file_no_contents():
    # if the file doesn't exist and the content doesn't exist, raise error
    with pytest.raises(RuntimeError):
        ensure_file("foo_bad_file.badfile", None, None)
