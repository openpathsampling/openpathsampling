#!/usr/bin/env python

# FileReader: an object for reading general data files. Basically, the
# assumptions are that you have columns, and comments use the octothorpe
# convention. The point is for this to be used in many, many scripts. This
# is to replace my crappy old read_acf file with a more uniform script that
# I'm not embarrassed to show the public.

# Written by David W.H. Swenson. Current version: 2013-12-**
import sys
class FileReader(object):
    
    def __init__(self, fname="", form=None, form_default="f"):
        self.fname = fname
        self.f = sys.stdin if self.fname == "" else open(self.fname, 'r')
        self.form_default = form_default
        self.form = []
        if (form != None):
            self.form = list(form)
        # form is the format string: each character corresponds to the type
        # for the associated column:
        #   s/c : string
        #   d/i : integer
        #   f/g : float
        # If not specified, we assume that you want the format given by
        # form_default
        self.formkey = {    "s" : str   , "c" : str,
                            "d" : int   , "i" : int,
                            "f" : float , "g" : float  }
        return

    def parse_line(self, line, form=""):
        from re import sub, split
        line = sub('\#.*', '', line)
        line = sub('^\s*', '', line)
        line = sub('\s*$', '', line)

        splitter = split('\s+', line)
        reslist = []
        if (len(splitter) > 1):
            for index in range(len(splitter)):
                if (len(self.form) > index):
                    # convert to the specified format
                    val = self.formkey[self.form[index]](splitter[index])
                else:
                    # conver to the default format
                    val = self.formkey[self.form_default](splitter[index])
                reslist.append(val)
        return reslist

    def next_line_as_array(self):
        myline = self.f.readline()
        splitter = self.parse_line(myline)
        return splitter

    def file_as_array(self, start=0):
        mydata = []
        if (self.f != sys.stdin):
            self.f.seek(start) # rewind the file, just in case
        for line in self.f:
            splitter= self.parse_line(line)
            if (len(splitter) > 0):
                mydata.append(splitter)
        return mydata

    def file_as_dict(self, nkeys=1):
        pass

    def columns_as_array(self, start=0):
        if (self.f != sys.stdin):
            self.f.seek(start)
        mydata = []
        ncol=-1
        for line in self.f:
            splitter = self.parse_line(line)
            if (ncol < 0):
                ncol = len(splitter)
                for col in range(ncol):
                    mydata.append([])

            for col in range(ncol):
                mydata[col].append(splitter[col])
        return mydata


