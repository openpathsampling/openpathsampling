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


def refresh_output(output_str, print_anyway=True, refresh=True,
                   output_stream=None, ipynb_display_only=False):
    if output_stream is None:
        output_stream = sys.stdout

    if is_ipynb or not ipynb_display_only or print_anyway:
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
