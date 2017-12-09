#!/usr/bin/env python

import os
import re

docs_prefix = "docs"

# We should have thing_to_wrap.md, wrap it in rst
things_to_wrap = [
    #"guides/absolute_beginners", "guides/power_users", "guides/user_levels"
]

def clean_filename(fname, strip_path=True):
    fname = re.sub('\_', '-', fname)
    if strip_path:
        fname = re.sub('.*\/', '', fname)
    return fname

def wrap_file(basename):
    md_file = basename + ".md"
    assert(os.path.isfile(md_file))
    clean_name = clean_filename(basename, strip_path=True)
    rst_file = basename + ".rst"
    lines = [".. _"+clean_name+":\n",
             ".. markdown:: " + docs_prefix + "/" + md_file]
    rst = open(rst_file, "w")
    for line in lines:
        rst.write(line)

if __name__ == "__main__":
    for thing in things_to_wrap:
        wrap_file(thing)
