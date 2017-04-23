from nose.tools import assert_equal
from openpathsampling.tools import *

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
