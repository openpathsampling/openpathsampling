#!/usr/bin/env python
"""
simple example script for running and testing notebooks.

Usage: `ipnbdoctest.py foo.ipynb [bar.ipynb [...]]`

Each cell is submitted to the kernel, and the outputs are compared with those stored in the notebook.

This is found in a gist under https://gist.github.com/minrk/2620735
"""
 
import os,sys,time
import base64
import re
 
from collections import defaultdict
from Queue import Empty

import difflib

try:
    from IPython.kernel import KernelManager
except ImportError:
    from IPython.zmq.blockingkernelmanager import BlockingKernelManager as KernelManager
 
from IPython.nbformat.current import reads, NotebookNode

def fold(name, block):
  return "travis_fold:start:" + name + "\r" + block + "travis_fold:end:" + name + "\r"

def fold_open(name):
    print "travis_fold:start:" + name

def fold_close(name):
    print "travis_fold:end:" + name

def indent(s, num = 4):
    lines = s.splitlines(True)
    lines = map(lambda s: ' ' * num + s, lines)
    return ''.join(lines)

def compare_png(a64, b64):
    """compare two b64 PNGs (incomplete)"""
    try:
        import Image
    except ImportError:
        pass
    adata = base64.decodestring(a64)
    bdata = base64.decodestring(b64)
    return True
 
def sanitize(s):
    """sanitize a string for comparison.
    
    fix universal newlines, strip trailing newlines, and normalize likely random values (memory addresses and UUIDs)
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
    
    return s
 
 
def consolidate_outputs(outputs):
    """consolidate outputs into a summary dict (incomplete)"""
    data = defaultdict(list)
    data['stdout'] = ''
    data['stderr'] = ''
    
    for out in outputs:
        if out.type == 'stream':
            data[out.stream] += out.text
        elif out.type == 'pyerr':
            print 'ERROR'
            data['pyerr'] = dict(ename=out.ename, evalue=out.evalue)
        else:
            for key in ('png', 'svg', 'latex', 'html', 'javascript', 'text', 'jpeg',):
                if key in out:
                    data[key].append(out[key])
    return data
 
 
def compare_outputs(test, ref, skip_compare=('png', 'traceback', 'latex', 'prompt_number', 'svg', 'html')):
    for key in ref:
        if key not in test:
#            print "missing key: %s != %s" % (test.keys(), ref.keys())
            return False
        elif key not in skip_compare:
            s1 = sanitize(test[key])
            s2 = sanitize(ref[key])
            if s1 != s2:

#                print "mismatch %s:" % key
                expected=s1.splitlines(1)
                actual=s2.splitlines(1)
                diff=difflib.unified_diff(expected, actual)

#                print ''.join(diff)
                return False
    return True
 

def run_cell(shell, iopub, cell):
    # print cell.input
    shell.execute(cell.input)
    # wait for finish, maximum 20s
    shell.get_msg(timeout=100)
    outs = []
    
    while True:
        try:
            msg = iopub.get_msg(timeout=0.2)
        except Empty:
            break
        msg_type = msg['msg_type']
        if msg_type in ('status', 'pyin'):
            continue
        elif msg_type == 'clear_output':
            outs = []
            continue
        
        content = msg['content']
        # print msg_type, content
        out = NotebookNode(output_type=msg_type)
        
        if msg_type == 'stream':
            out.stream = content['name']
            out.text = content['data']
        elif msg_type in ('display_data', 'pyout'):
            out['metadata'] = content['metadata']
            for mime, data in content['data'].iteritems():
                attr = mime.split('/')[-1].lower()
                # this gets most right, but fix svg+html, plain
                attr = attr.replace('+xml', '').replace('plain', 'text')
                setattr(out, attr, data)
            if msg_type == 'pyout':
                out.prompt_number = content['execution_count']
        elif msg_type == 'pyerr':
            out.ename = content['ename']
            out.evalue = content['evalue']
            out.traceback = content['traceback']
        else:
            print "unhandled iopub msg:", msg_type
        
        outs.append(out)
    return outs
    
 
def test_notebook(nb, nbname = ''):

    print "starting kernel ...",
    km = KernelManager()
    km.start_kernel(extra_arguments=['--pylab=inline'], stderr=open(os.devnull, 'w'))
    print "ok"

    try:
        kc = km.client()
        kc.start_channels()
        iopub = kc.iopub_channel
    except AttributeError:
        # IPython 0.13
        kc = km
        kc.start_channels()
        iopub = kc.sub_channel
    shell = kc.shell_channel
    
    # run %pylab inline, because some notebooks assume this
    # even though they shouldn't
    shell.execute("pass")
    shell.get_msg()
    while True:
        try:
            iopub.get_msg(timeout=1)
        except Empty:
            break
    
    successes = 0 # passed without differences
    failures = 0 # not passed, but executed with Python error
    errors = 0 # IPython errors (not even executed)
    differences = 0 # passed but with differences

    nbname = nbname.replace(" ", "_").lower()

    show_md = True

    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == 'markdown':
                for line in cell.source.splitlines():
                    # only print headlines in markdown
                    if line.startswith('#'):
                        print line

            if cell.cell_type == 'heading':
                print '#' * cell.level, cell.source

            if cell.cell_type != 'code':
                continue

            # If code cell then continue with checking it

            if hasattr(cell, 'prompt_number'):
                print nbname + '.' + 'In [%3i]' % cell.prompt_number, '...',
            else:
                print nbname + '.' + 'In [???]', '...',
            command = None
            try:
                first_line = cell.input.splitlines()[0]
                if first_line.startswith('#!'):
                    command = first_line[2:].strip()
                outs = run_cell(shell, iopub, cell)
            except Exception as e:
                # Internal IPython error occurred (might still be that the cell did not execute correctly)
                print 'ERROR'
                print repr(e)
                errors += 1
                continue

#            print outs

            if command == 'skip':
                print 'pass'
                continue
            
            failed = False
            diff = False
            for out, ref in zip(outs, cell.outputs):
                if out.output_type == 'pyerr':
                    # An python error occurred. Cell is not completed correctly
                    print "FAIL"
                    fold_open('FAIL')
                    print '>>> ' + out.ename + ' ("' + out.evalue + '")'
                    for idx, trace in enumerate(out.traceback):
                        print indent(trace, 4)
                    fold_close('FAIL')
                    failed = True
                else:
                    if not compare_outputs(out, ref):
                        # Output is different than the one in the notebook. This might be okay.
                        diff = True
            if failed:
                failures += 1

            if diff:
                differences += 1
                print 'DIFF'

            if not failed and not diff:
                successes += 1
                print 'ok'

    print
    print "tested notebook %s" % nb.metadata.name
    print "    %3i cells successfully replicated [ok]" % successes
    print "    %3i cells mismatched output [DIFF]" % differences
    print "    %3i cells ran with python errors [FAIL]" % failures
    print "    %3i cells failed to even execute (IPython error) [ERROR]" % errors
    print
    print "shutting down kernel ...",
    kc.stop_channels()
    km.shutdown_kernel()
    print "ok"
    del km

    if errors != 0:
        # There is an error. Stop testing the next files and quit with an error
        print 'errors occurred - exiting.'
        exit(1)
 
if __name__ == '__main__':
    for ipynb in sys.argv[1:]:
        print 'testing ipython notebook : "%s"' % ipynb
        with open(ipynb) as f:
            nb = reads(f.read(), 'json')
        test_notebook(nb, nbname = ipynb)

    print 'done.'

    exit(0)