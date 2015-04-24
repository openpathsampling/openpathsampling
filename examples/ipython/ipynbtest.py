#!/usr/bin/env python
"""
Simple example script for running and testing IPython notebooks.

usage: ipnbdoctest.py [-h] [--timeout TIMEOUT] [--strict] [--fail-if-timeout]
                      [--show-diff]
                      file.ipynb

Run all cells in an ipython notebook as a test and check whether these
successfully execute and compares their output to the one inside the notebook

positional arguments:
  file.ipynb         the notebook to be checked

optional arguments:
  -h, --help         show this help message and exit
  --timeout TIMEOUT  the default timeout time in seconds for a cell
                     evaluation. Default is 300s.
  --strict           if set to true then the default test is that cell have to
                     match otherwise a diff will not be considered a failed
                     test
  --fail-if-timeout  if set to true then a timeout is considered a failed test
  --show-diff        if set to true differences in the cell are shown in
                     `diff` style

Each cell is submitted to the kernel, and the outputs are compared with those
stored in the notebook.

This version need IPython 3.0.0 and makes use of some nice features. It can
handle notebooks of version 3 (IPython 2) and version 4 (IPython 3)

The original is found in a gist under https://gist.github.com/minrk/2620735
"""

import os,sys
import re
import argparse
 
from Queue import Empty
import difflib

# IPython 3.0.0
from IPython.kernel.manager import KernelManager
import IPython

# Allows to read from all notebook versions
from IPython.nbformat.reader import reads
from IPython.nbformat import NotebookNode

import uuid

class TravisConsole(object):
    """
    A wrapper class to allow easier output to the console especially for travis
    """
    def __init__(self):
        self.stream = sys.stdout
        self.linebreak = '\n'
        self.fold_count = dict()
        self.fold_stack = dict()

        # very complicated, maybe too? complicated
        self.fold_uuid = str(uuid.uuid4()).split('-')[1]

    def fold_open(self, name):
        """
        open the travis fold with given name

        Parameters
        ----------
        name : string
            the name of the fold
        """
        if name not in self.fold_count:
            self.fold_count[name] = 0
            self.fold_stack[name] = []

        self.fold_count[name] += 1
        fold_name = self.fold_uuid + '.' + name.lower() \
                    + '.' + str(self.fold_count[name])

        self.fold_stack[name].append(fold_name)

        self.writeln("travis_fold:start:" + fold_name)

    def fold_close(self, name):
        """
        close the travis fold with given name

        Parameters
        ----------
        name : string
            the name of the fold
        """
        fold_name = self.fold_uuid + '.' + name.lower() \
                    + '.' + str(self.fold_count[name])
        self.writeln("travis_fold:end:" + fold_name)

    def _indent(self, s, num = 4):
        lines = s.splitlines(True)
        lines = map(lambda s: ' ' * num + s, lines)
        return ''.join(lines)

    def writeln(self, s, indent = 0):
        """write a string with linebreak at the end

        Parameters
        ----------
        s : string
            the string to be written to travis console
        indent : int, default 0
            if non-zero add the number of space before each line
        """
        self.write(s, indent)
        if s[-1] != '\n':
            # make a line break if there is none present
            self.br()

    def br(self):
        """
        Write a linebreak
        """
        self.stream.write(self.linebreak)

    def write(self, s, indent = 0):
        """write a string to travis output with possible indention

        Parameters
        ----------
        s : string
            the string to be written to travis console
        indent : int, default 0
            if non-zero add the number of space before each line
        """
        if indent > 0:
            self.stream.write(self._indent(s, indent))
        else:
            self.stream.write(s)

        self.stream.flush()

    def red(self, s):
        """format a string to be red in travis output

        Parameters
        ----------
        s : string
            string to be colored

        Returns
        -------
        string
            the colored string
        """
        RED = '\033[31m'
        DEFAULT = '\033[39m'
        return RED + s + DEFAULT

    def green(self, s):
        """format a string to be green in travis output

        Parameters
        ----------
        s : string
            string to be colored

        Returns
        -------
        string
            the colored string
        """
        GREEN = '\033[32m'
        DEFAULT = '\033[39m'
        return GREEN + s + DEFAULT

    def blue(self, s):
        """format a string to be blue in travis output

        Parameters
        ----------
        s : string
            string to be colored

        Returns
        -------
        string
            the colored string
        """
        BLUE = '\033[36m'
        DEFAULT = '\033[39m'
        return BLUE + s + DEFAULT

    def format_diff(self, diff):
        """format a list of diff commands for travis output

        this will remove empty lines, lines starting with `?` and
        add coloring depending on whether a line starts with `+` or `-`

        Parameters
        ----------
        diff : list of diff (string)
            the list of diff commands to be formatted

        Returns
        -------
        string
            a string representation of all diffs
        """
        colored_diffs = []
        for line in diff:
            if line[0] == '-':
                colored_diffs.append(self.red(line))
            elif line[0] == '+':
                colored_diffs.append(self.green(line))
            else:
                colored_diffs.append(line)

        # remove unnecessary linebreaks
        colored_diffs = [ line.replace('\n', '') for line in colored_diffs]

        # remove line we do not want
        colored_diffs = [ line for line in colored_diffs
                            if len(line) > 0 and line[0] != '?']

        return '\n'.join(colored_diffs)


class IPyTestConsole(TravisConsole):
    """
    Add support for different output results
    """
    def __init__(self):
        super(IPyTestConsole, self).__init__()

        self.default_results = {
            'success' : True,       # passed without differences
            'kernel' : False,       # kernel (IPYTHON) error occurred
            'error' : False,        # errors during execution
            'timeout' : True,       # kernel run timed out
            'diff' : True,          # passed, but with differences in the output
            'skip' : True,          # cell has been skipped
            'ignore' : True         # cell has been executed, but not compared
        }

        self.pass_count = 0
        self.fail_count = 0

        self.reset()

        self.last_fail = False

    def reset(self):
        self.result_count = { key : 0 for key in self.default_results.keys() }
        self.pass_count = 0
        self.fail_count = 0


    def write_result(self, result, okay_list = None):
        """write final result of test

        this keeps track of the result types
        """
        my_list = self.default_results.copy()
        if okay_list is not None:
            my_list.update(okay_list)

        if my_list[result]:
            self.write(self.green('ok'))
            self.pass_count += 1
            self.last_fail = False
        else:
            self.write(self.red('fail'))
            self.fail_count += 1
            self.last_fail = True

        self.writeln(' [' + result + ']')
        self.result_count[result] += 1


class IPyKernel(object):
    """
    A simple wrapper class to run cells in an IPython Notebook.

    Notes
    -----
    - Use `with` construct to properly instantiate
    - IPython 3.0.0+ is assumed for this version

    """

    def __init__(self, console = None, nb_version=4):
        # default timeout time is 60 seconds
        self.default_timeout = 60
        self.extra_arguments = ['--pylab=inline']
        self.nb_version = nb_version

    def __enter__(self):
        self.km = KernelManager()
        self.km.start_kernel(
            extra_arguments=self.extra_arguments,
            stderr=open(os.devnull, 'w')
        )

        self.kc = self.km.client()
        self.kc.start_channels()

        self.iopub = self.kc.iopub_channel
        self.shell = self.kc.shell_channel

        # run %pylab inline, because some notebooks assume this
        # even though they shouldn't

        self.shell.send("pass")
        self.shell.get_msg()
        while True:
            try:
                self.iopub.get_msg(timeout=1)
            except Empty:
                break

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.kc.stop_channels()
        self.km.shutdown_kernel()
        del self.km

    def run(self, cell, timeout = None):
        """
        Run a notebook cell in the IPythonKernel

        Parameters
        ----------
        cell : IPython.notebook.Cell
            the cell to be run
        timeout : int or None (default)
            the time in seconds after which a cell is stopped and assumed to
            have timed out. If set to None the value in `default_timeout`
            is used

        Returns
        -------
        list of outs
            a list of NotebookNodes of the returned types. This is
            similar to the list of outputs generated when a cell is run
        """

        use_timeout = self.default_timeout

        if timeout is not None:
            use_timeout = timeout

        if hasattr(cell, 'input'):
            self.kc.execute(cell.input)
        elif hasattr(cell, 'source'):
            self.kc.execute(cell.source)
        else:
            raise AttributeError('No source/input key')

        self.shell.get_msg(timeout=use_timeout)
        outs = []

        while True:
            try:
                msg = self.iopub.get_msg(timeout=0.5)
            except Empty:
                break
            msg_type = msg['msg_type']
            if msg_type in ('status', 'pyin', 'execute_input'):
                continue
            elif msg_type == 'clear_output':
                outs = []
                continue

            content = msg['content']
            out = NotebookNode(output_type=msg_type)

            if msg_type == 'stream':
                out.name = content['name']
                out.text = content['text']
            elif msg_type in ('display_data', 'pyout', 'execute_result'):
                if hasattr(content, 'execution_count'):
                    out['execution_count'] = content['execution_count']
                else:
                    out['execution_count'] = None
                out['data'] = content['data']
                out['metadata'] = content['metadata']

            elif msg_type == 'error':
                out.ename = content['ename']
                out.evalue = content['evalue']
                out.traceback = content['traceback']
            else:
                print "unhandled iopub msg:", msg_type

            outs.append(out)

        return outs

    def sanitize(self, s):
        """sanitize a string for comparison.

        fix universal newlines, strip trailing newlines, and normalize likely
        random values (memory addresses and UUIDs)
        """
        if not isinstance(s, basestring):
            return s
        # normalize newline:
        s = s.replace('\r\n', '\n')

        # ignore trailing newlines (but not space)
        s = s.rstrip('\n')

        # normalize hex addresses:
        s = re.sub(r'0x[a-f0-9]+', '0xFFFFFFFF', s)

        # normalize UUIDs:
        s = re.sub(r'[a-f0-9]{8}(\-[a-f0-9]{4}){3}\-[a-f0-9]{12}', 'U-U-I-D', s)

        # fix problem with

        return s

    def compare_outputs(
            self,
            test,
            ref,
            skip_compare=('traceback', 'latex', 'execution_count')
    ):
        """
        Compare two lists of `NotebookNode`s

        Parameters
        ----------
        test : list of `NotebookNode`
            the list of be tested generated by the kernel
        ref : list of `NotebookNode`
            the reference list read from the notebook
        skip_compare : list of str
            a list of strings that name node types that are not to be tested

        Returns
        -------
        bool
            is True if both lists are different
        list of diff
            a list of diff (str) the represent the differences
        """
        diff = False
        diff_list = []

        if self.nb_version == 4:
            for key in ref:
                if key not in test:
                    return True, [ "missing key: %s != %s" %
                                        (test.keys(), ref.keys()) ]

                elif key not in skip_compare:
                    if key == 'data':
                        for data_key in test[key]:
                            my_diff = self.do_diff(
                                data_key,
                                test[key],
                                ref[key])

                            if my_diff is not None:
                                diff_list += my_diff
                                diff = True

                    else:
                        # can this happen?
                        my_diff = self.do_diff(key, test, ref)
                        if my_diff is not None:
                            diff_list += my_diff
                            diff = True

        return diff, diff_list

    def do_diff(self, key, test_cell, ref_cell):
        """
        Compare the key of two dicts

        Parameters
        ----------
        key : string
            the key to be compared
        test_cell : dict
            a dict with `key` as a key of string value
        ref_cell : dict
            a dict with `key` as a key of string value

        Returns
        -------
        list of diff (str)
            a list of diff representing the differences
        """

        if hasattr(ref_cell, key):
            s1 = self.sanitize(ref_cell[key])
        else:
            s1 = ''

        if hasattr(test_cell, key):
            s2 = self.sanitize(test_cell[key])
        else:
            s2 = ''

        if key in ['image/png', 'image/svg', 'image/svg+xml']:
            if s1 != s2:
                return ['>>> diff in %s (size new : %d vs size old : %d )' %
                            (key, len(s1), len(s2) )]
        else:
            if s1 != s2:
                expected=s1.splitlines(1)
                actual=s2.splitlines(1)
                diff=difflib.ndiff(expected, actual)

                return [ '>>> diff in ' + key ] + list(diff)

        return None

    def get_commands(self, cell):
        """
        Extract potential commands from the first line of a cell

        if a code cell starts with the hashbang `#!` it can be followed by
        a comma separated list of commands. Each command can be
        a single key `skip`
        or
        a key/value pair separated by a colon `timeout:[int]`

        Parameters
        ----------
        cell : a NotebookCell
            the cell to be examined

        Returns
        -------
        dict
            a dict of key/value pairs. For a single command the value is `True`
        """
        commands = {}
        source = self.get_source(cell)
        if source is not None:
            lines = source.splitlines()
            if len(lines) > 0:
                first_line = lines[0]
                if first_line.startswith('#!'):
                    txt = first_line[2:].strip()

                    parts = txt.split(',')
                    for part in parts:
                        subparts = part.split(':')
                        if len(subparts) == 1:
                            commands[subparts[0].strip().lower()] = True
                        elif len(subparts) == 2:
                            commands[subparts[0].strip().lower()] = subparts[1]

        return commands

    def get_source(self, cell):
        """
        get the source code of a cell

        Notes
        -----
        This is legacy of IPython 2/3 conversion.

        Parameters
        ----------
        cell : a NotebookCell
            the cell to be examined

        Returns
        -------
        string
            the source code

        """
        if cell.cell_type == 'code':
            if hasattr(cell, 'input'):
                return cell.input
            elif hasattr(cell, 'source'):
                return cell.source
            else:
                return None

    def is_empty_cell(self, cell):
        """
        Check if a cell has no code

        Parameters
        ----------
        cell : a NotebookCell
            the cell to be examined

        Returns
        -------
        bool
            True if the cell has no code, False otherwise
        """
        source = self.get_source(cell)
        if source is None or source == '':
            return True
        else:
            return False

###############################################################################
##  MAIN
###############################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run all cells in an ipython notebook as a test and ' +
                    'check whether these successfully execute and ' +
                    'compares their output to the one inside the notebook')

    parser.add_argument(
                    'file',
                    metavar='file.ipynb',
                    help='the notebook to be checked',
                    type=str)

    parser.add_argument('--timeout', dest='timeout',
                    type=int, default=300,
                    help='the default timeout time in seconds for a cell ' +
                        'evaluation. Default is 300s.')

    parser.add_argument('--rerun-if-timeout', dest='rerun',
                    type=int, default=2, nargs='?',
                    help='if set then a timeout in a cell will cause to run ' +
                         'the. Default is 2 (means make upto 3 attempts)')

    parser.add_argument('--restart-if-fail', dest='restart',
                    type=int, default=2, nargs='?',
                    help='if set then a fail in a cell will cause to restart ' +
                         'the full notebook!. Default is 0 (means NO rerun).' +
                         'Use this with care.')

    parser.add_argument('--strict', dest='strict',
                    action='store_true',
                    default=False,
                    help='if set to true then the default test is that cell ' +
                        'have to match otherwise a diff will not be ' +
                        'considered a failed test')

    parser.add_argument(
                    '--fail-if-timeout',
                    dest='no_timeout', action='store_true',
                    default=False,
                    help='if set to true then a timeout is considered a ' +
                         'failed test')

    parser.add_argument(
                    '--show-diff',
                    dest='show_diff',
                    action='store_true',
                    default=False,
                    help='if set to true differences in the cell are shown ' +
                         'in `diff` style')

    parser.add_argument(
                    '-v, --verbose',
                    dest='verbose', action='store_true',
                    default=False,
                    help='if set then text output is send to the ' +
                         'console.')


    args = parser.parse_args()
    ipynb = args.file
    verbose = args.verbose

    tv = IPyTestConsole()

    if args.strict:
        tv.default_results['diff'] = False

    if args.no_timeout:
        tv.default_results['timeout'] = False

    tv.writeln('testing ipython notebook : "%s"' % ipynb)
    tv.fold_open('ipynb')

    timeout_rerun = args.rerun
    fail_restart = args.restart

    with open(ipynb) as f:
        nb = reads(f.read())
        # Convert all notebooks to the format IPython 3.0.0 uses to
        # simplify comparison
        nb = IPython.nbformat.convert(nb, 4)

    notebook_restart = True
    notebook_run_count = 0

    while (notebook_restart):
        notebook_restart = False
        notebook_run_count += 1

        tv.reset()
        tv.write("starting kernel ... ")
        with IPyKernel() as ipy:
            ipy.default_timeout = args.timeout
            tv.writeln("ok")

            nbs = ipynb.split('.')
            nb_class_name = nbs[1] + '.' + nbs[0].replace(" ", "_")

            tv.br()

            if hasattr(nb, 'worksheets'):
                ws = nb.worksheets[0]
            else:
                ws = nb

            for cell in ws.cells:

                if notebook_restart:
                    # if we restart anyway skip all remaining cells
                    continue

                if cell.cell_type == 'markdown':
                    for line in cell.source.splitlines():
                        # only tv.writeln(headlines in markdown
                        # TODO: exclude # comments in code blocks
                        if line.startswith('#') and line[1] != '!':
                            tv.writeln(line)

                if cell.cell_type == 'heading':
                    tv.writeln('#' * cell.level + ' ' + cell.source)

                if cell.cell_type != 'code':
                    continue

                # if code cell then continue

                if ipy.is_empty_cell(cell):
                    # empty cell will not be tested
                    continue

                if hasattr(cell, 'prompt_number'):
                    tv.write(nb_class_name + '.' + 'In [%3i]' %
                             cell.prompt_number + ' ... ')
                elif hasattr(cell, 'execution_count') and \
                                cell.execution_count is not None:
                    tv.write(nb_class_name + '.' + 'In [%3i]' %
                             cell.execution_count + ' ... ')
                else:
                    tv.write(nb_class_name + '.' + 'In [???]' + ' ... ')

                commands = ipy.get_commands(cell)

                result = 'success'

                timeout = ipy.default_timeout

                if 'skip' in commands:
                    tv.write_result('skip')
                    continue

                cell_run_count = 0
                cell_run_again = True
                cell_passed = True

                while cell_run_again:
                    cell_run_count += 1
                    cell_run_again = False
                    cell_passed = True

                    try:
                        if 'timeout' in commands:
                            outs = ipy.run(cell, timeout=int(commands['timeout']))
                        else:
                            outs = ipy.run(cell)

                    except Exception as e:
                        # Internal IPython error occurred (might still be
                        # that the cell did not execute correctly)
                        if 'ignore' not in commands:
                            cell_passed = False
                            if repr(e) == 'Empty()':
                                # Assume it has been timed out!
                                if cell_run_count <= timeout_rerun:
                                    cell_run_again = True
                                    tv.write('timeout [retry #%d] ' % cell_run_count)
                                else:
                                    tv.write_result('timeout')
                                # tv.writeln('>>> TimeOut (%is)' % args.timeout)
                            else:
                                tv.write_result('kernel')
                                tv.fold_open('ipynb.kernel')
                                tv.writeln('>>> ' + out.ename + ' ("' + out.evalue + '")')
                                tv.writeln(repr(e), indent=4)
                                tv.fold_close('ipynb.kernel')
                        else:
                            if repr(e) == 'Empty()':
                                # Assume it has been timed out!
                                tv.write('timeout / ')
                                tv.write_result('ignore')
                            else:
                                tv.write('kernel / ')
                                tv.write_result('ignore')
                                tv.fold_open('ipynb.kernel')
                                tv.writeln('>>> ' + out.ename + ' ("' + out.evalue + '")')
                                tv.writeln(repr(e), indent=4)
                                tv.fold_close('ipynb.kernel')

                if not cell_passed:
                    if tv.last_fail and notebook_run_count <= fail_restart:
                        notebook_restart = True

                    continue

                failed = False
                diff = False
                diff_str = ''
                out_str = ''
                for out, ref in zip(outs, cell.outputs):
                    if out.output_type == 'error':
                        # An python error occurred. Cell is not completed correctly
                        if 'ignore' not in commands:
                            tv.write_result('error')
                        else:
                            tv.write('error / ')
                            tv.write_result('ignore')

                        tv.fold_open('ipynb.error')
                        tv.writeln('>>> ' + out.ename + ' ("' + out.evalue + '")')

                        for idx, trace in enumerate(out.traceback[1:]):
                            tv.writeln(trace, indent=4)

                        tv.fold_close('ipynb.error')
                        failed = True
                    else:
                        this_diff, this_str = ipy.compare_outputs(out, ref)
                        if 'verbose' in commands or verbose:
                            if 'data' in out:
                                for key,value in out.data.iteritems():
                                    if 'text' in key:
                                        out_str += value + '\n'

                        if this_diff:
                            # Output is different than the one in the notebook.
                            diff_str += tv.format_diff(this_str)
                            diff = True

                if diff and not failed:
                    if 'ignore' not in commands:
                        if 'strict' in commands:
                            # strict mode means a difference will fail the test
                            tv.write_result('diff', okay_list={ 'diff' : False })
                        elif 'lazy' in commands:
                            # lazy mode means a difference will pass the test
                            tv.write_result('diff', okay_list={ 'diff' : True })
                        else:
                            # use defaults
                            tv.write_result('diff')
                    else:
                        tv.write('diff / ')
                        tv.write_result('ignore')

                    if args.show_diff:
                        tv.fold_open('ipynb.diff')
                        tv.writeln(diff_str)
                        tv.fold_close('ipynb.diff')

                if not failed and not diff:
                    tv.write_result('success')

                if tv.last_fail and notebook_run_count <= fail_restart:
                    # we had a fail so restart the whole notebook
                    notebook_restart = True

                if out_str != '' and 'quiet' not in commands:
                    tv.fold_open('ipynb.out')
                    tv.writeln(out_str)
                    tv.fold_close('ipynb.out')

            tv.br()
            tv.writeln("  testing results")
            tv.writeln("  ===============")
            if tv.pass_count > 0:
                tv.writeln("    %3i cells passed [" %
                           tv.pass_count + tv.green('ok') + "]" )
            if tv.fail_count > 0:
                tv.writeln("    %3i cells failed [" %
                           tv.fail_count + tv.red('fail') + "]" )

            tv.br()
            tv.writeln("  %3i cells successfully replicated [success]" %
                       tv.result_count['success'])
            tv.writeln("  %3i cells had mismatched outputs [diff]" %
                       tv.result_count['diff'])
            tv.writeln("  %3i cells timed out during execution [time]" %
                       tv.result_count['timeout'])
            tv.writeln("  %3i cells ran with python errors [error]" %
                       tv.result_count['error'])
            tv.writeln("  %3i cells have been run without comparison [ignore]" %
                       tv.result_count['ignore'])
            tv.writeln("  %3i cells failed to even run (IPython error) [kernel]" %
                       tv.result_count['kernel'])
            tv.writeln("  %3i cells have been skipped [skip]" %
                       tv.result_count['skip'])

            if notebook_restart:
                tv.br()
                tv.writeln(
                    tv.red("  attempt #%d of max %d failed, restarting notebook!" % (notebook_run_count, fail_restart + 1))
                )

            tv.br()
            tv.write("shutting down kernel ... ")


        tv.writeln('ok')

    tv.fold_close('ipynb')

    if tv.fail_count != 0:
        tv.writeln(tv.red('some tests not passed.'))
        exit(1)
    else:
        tv.writeln(tv.green('all tests passed.'))
        exit(0)