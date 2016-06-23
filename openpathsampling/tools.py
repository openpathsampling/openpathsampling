import sys

__author__ = 'Jan-Hendrik Prinz'


def refresh_output(output_str, print_anyway=True, refresh=True,
                   output_stream=None):
    if output_stream is None:
        output_stream = sys.stdout
    try:
        import IPython.display
    except ImportError:
        if print_anyway:
            output_stream.write(output_str)
    else:
        if refresh:
            IPython.display.clear_output(wait=True)
        output_stream.write(output_str)
    sys.stdout.flush()
