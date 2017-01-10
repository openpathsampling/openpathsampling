import sys

__author__ = 'Jan-Hendrik Prinz'


try:
    import IPython
    import IPython.display

    def in_ipynb():
        try:
            ipython = get_ipython()

            import IPython.terminal.interactiveshell
            import ipykernel.zmqshell

            if isinstance(ipython, IPython.terminal.interactiveshell.TerminalInteractiveShell):
                # we are running inside an IPYTHON console
                return False
            elif isinstance(ipython, ipykernel.zmqshell.ZMQInteractiveShell):
                # we run in an IPYTHON notebook
                return True
            else:
                return False
        except NameError:
            # No IPYTHON
            return False
        except:
            # No idea, but we should not fail because of that
            return False

    is_ipynb = in_ipynb()

except ImportError:
    is_ipynb = False

last_output = None


def refresh_output(output_str, print_anyway=True, refresh=True,
                   output_stream=None, ipynb_display_only=False):

    global last_output

    if output_stream is None:
        output_stream = sys.stdout

    if not output_str:
        return

    if refresh:
        if is_ipynb:
            IPython.display.clear_output(wait=True)
        else:
            if last_output is not None:
                lines = len(last_output.split('\n'))
                CURSOR_UP_ONE = '\x1b[1A'
                ERASE_LINE = '\x1b[2K'

                output_stream.write((CURSOR_UP_ONE + ERASE_LINE) * (lines - 1))

            last_output = output_str
    else:
        last_output = None

    output_stream.write(output_str)
    output_stream.flush()


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
