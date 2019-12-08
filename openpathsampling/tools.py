import sys
import os
import hashlib

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
        elif output_stream is sys.stdout:
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


def pretty_print_seconds(seconds, n_labels=0, separator=" "):
    """
    Converts a number of seconds to a readable time string.

    Parameters
    ----------
    seconds : float or int
        number of seconds to represent, gets rounded with int()
    n_labels : int
        number of levels of label to show. For example, if n_labels=1,
        result will show the first of days, hours, minutes, seconds with
        greater than one. This allows you to round, e.g, 1 day 2 hours
        3 minutes 4 seconds to 1 day 2 hours (with n_labels=2). Default
        value of 0 gives all levels. If n_labels is negative, then the last
        value is shown as a decimal, instead of an int, with 2 decimals of
        precision.
    separator : string
        separator between levels of the time decomposition
    """
    ordered_keys = ['day', 'hour', 'minute', 'second']
    divisors = {
        'day': 86400,
        'hour': 3600,
        'minute': 60,
        'second': 1
    }

    s = int(seconds)

    def decompose_seconds(s):
        parts = {}
        fractional_parts = {}
        for k in ordered_keys:
            fractional_parts[k] = float(s) / divisors[k]
            parts[k], s = divmod(s, divisors[k])
        return parts, fractional_parts

    def make_seconds(parts):
        return sum([parts[p] * divisors[p] for p in parts.keys()])

    parts, fractional_parts = decompose_seconds(s)

    decimalize_final = (n_labels < 0)

    first_key = "second"
    for key in ordered_keys:
        if parts[key] > 0:
            first_key = key
            break
    first_key_index = ordered_keys.index(first_key)

    n_labels_real = len(ordered_keys) - first_key_index
    if n_labels != 0 and abs(n_labels) < n_labels_real:
        n_labels_real = abs(n_labels)

    max_label_index = first_key_index + len(ordered_keys) - n_labels_real
    if max_label_index >= len(ordered_keys):
        max_label_index = len(ordered_keys) - 1
    max_label = ordered_keys[max_label_index]

    if first_key != "second" and n_labels > 0:
        # round it!
        if fractional_parts[max_label] - parts[max_label] >= 0.5:
            parts[max_label] += 1
        else:
            pass
        for key in ordered_keys[max_label_index + 1:]:
            parts[key] = 0

    new_s = make_seconds(parts)
    parts, frac_parts = decompose_seconds(new_s)

    part_labels = {k: k if parts[k] == 1 else k + "s"
                   for k in ordered_keys}

    label_count = 0
    output_str = ""
    for key in ordered_keys[first_key_index:]:
        part = parts[key]
        label_str = part_labels[key]
        frac = fractional_parts[key]
        if part > 0 and label_count < n_labels_real - 1:
            output_str += str(part) + " " + label_str + separator
            label_count += 1
        elif label_count == n_labels_real - 1:
            if decimalize_final and key != 'second':
                output_str += "%.2f %s" % (frac, key+'s')
            else:
                output_str += str(part) + " " + label_str
            label_count += 1

    return output_str


def progress_string(n_steps_completed, n_steps_total, time_elapsed):
    """
    String to report on simulation progress.

    Assumes that the average time per Monte Carlo/trajectory-level step is
    constant throughout the simulation. These are trajectory-level steps in
    to distinguish them from molecular dynamics (time) steps -- in path
    sampling these are Monte Carlo steps; in order simulations (committor)
    they might not be.

    Parameters
    ----------
    n_steps_completed : int
        number of (Monte Carlo/trajectory-level) steps already completed
    n_steps_total : int
        total number of (Monte Carlo/trajectory-level) step in simulation
    time_elapsed : float-like
        time elapsed in the simulation, in seconds

    Returns
    -------
    str :
        string to output describing simulation progress, including estimated
        time remaining
    """
    try:
        time_per_step = time_elapsed / n_steps_completed
    except ZeroDivisionError:
        return "Starting simulation...\nWorking on first step\n"
    time_to_finish = (n_steps_total - n_steps_completed) * time_per_step
    output_str = (
        "Running for " + pretty_print_seconds(time_elapsed) + " - "
        + "%5.2f seconds per step\n" % (time_per_step)
        + "Estimated time remaining: %s\n" % (
            pretty_print_seconds(time_to_finish, n_labels=-2)
        )
    )
    return output_str


def ensure_file(filename, old_contents=None, old_hash=None):
    """Ensure that the existing file matches the old contents.

    If the file exists and we don't know the old contents/hash, trust the
    file (probably first initialization).
    If the file does not exist, this write the file. If the file exists,
    check that its contents (based on hash digest). If these match the
    original contents, we're fine. If not, raise an error.

    Parameters
    ----------
    filename : str
        filename
    old_contents : Union[str, None]
        expected file contents; if not given, assume we trust whatever
        content is in the file
    old_hash : Union[str, None]
        expected hash; if not given, we generate a hash from the old
        contents

    Returns
    -------
    contents : str
        file contents
    hashed : str
        hash of the file contents
    """
    hash_function = lambda text: hashlib.sha1(text.encode('utf-8')).digest()

    if old_hash is None and old_contents is not None:
        old_hash = hash_function(old_contents)

    if not os.path.exists(filename):
        # write the file if it doesn't exist
        if old_contents is not None:
            with open(filename, 'w') as f:
                f.write(old_contents)
        else:
            raise RuntimeError("No contents to write missing file " +
                               str(filename))

    with open(filename, mode='r') as f:
        contents = f.read()

    hashed = hash_function(contents)

    if old_hash and hashed != old_hash:
        raise RuntimeError("Existing file " + str(filename) + " does not"
                           + " match stored file.")

    return contents, hashed
