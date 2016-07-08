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


# a little code snippet to wrap strings around for nicer output
# idea found @ http://www.saltycrane.com/blog/2007/09/python-word-wrap-function/

def word_wrap(string, width=80):
    lines = string.split('\n')

    lines = [x.rstrip() for x in lines]
    result = []
    for line in lines:
        while len(line) > width:
            marker = width - 1
            while not line[marker].isspace():
                marker -= 1

            result.append(line[0:marker])
            line = line[marker + 1:]

        result.append(line)
    return '\n'.join(result)