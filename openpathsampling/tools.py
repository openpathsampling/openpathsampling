import sys

__author__ = 'Jan-Hendrik Prinz'


def refresh_output(output_str, print_anyway=True, refresh=True):
    try:
        import IPython.display
    except ImportError:
        if print_anyway:
            print(output_str)
    else:
        if refresh:
            IPython.display.clear_output(wait=True)
        print(output_str)
    sys.stdout.flush()