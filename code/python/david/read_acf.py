# This is rapidly becoming the overall file-reading file. And eventually, we
# need to support input from stdin.

def parse_line(line, nonfloat=0):
    from re import sub, split
    line = sub('\#.*', '', line)
    line = sub('^\s*', '', line)
    line = sub('\s*$', '', line)

    splitter = split('\s+', line)
    split_float = []
    if (len(splitter) > 1):
        for index in range(nonfloat):
            split_float.append(splitter[index])
        for index in range(nonfloat,len(splitter)):
            split_float.append(float(splitter[index]))

    return split_float

def read_file_as_array(fname, nonfloat=0):
    mydata = []
    myf = sys.stdin if fname=="" else open(fname, 'r')
    for line in myf:
        splitter = parse_line(line,nonfloat)
        if (len(splitter) > 0):
            mydata.append(splitter)
    return mydata

def next_line_as_array(f,nonfloat=0):
    myline = f.readline()
    splitter = parse_line(myline,nonfloat)
    return splitter


def read_file_as_dict(fname,nonfloat=0):
    file = open(fname, 'r')
    dict = { }
    for line in file:
        splitter = parse_line(line,nonfloat)
        if (len(splitter) > 0):
            key = splitter[0]
            vals = []
            for i in range(1,len(splitter)):
                vals.append(splitter[i])

            dict[key] = vals

    return dict

# this is for the case that the first N items are a combined key
def read_file_as_Nkey_list(fname, n):
    myf = open(fname, 'r')
    keys = []
    vals = []
    for line in myf:
        splitter = parse_line(line)
        if (len(splitter) > 0):
            key_i = []
            val_i = []
            for j in range(n):
                key_i.append(splitter[j])
            for j in range(n,len(splitter)):
                val_i.append(splitter[j])
            keys.append(key_i)
            vals.append(val_i)
    return (keys, vals)

## reads in an autocorrelation file
#@param filename autocorrelation file name
def read_acf(filename):
    from re import sub, split, search

    # open file
    myfile = open(filename, "r")
    
    time = []
    acf = []
    for line in myfile:
        # clear comments; clear initial and final whitespace
        line = sub('\#.*', '', line)
        line = sub('^\s*', '', line)
        line = sub('\s*$', '', line)

        splitter = split('\s+', line)
        if (len(splitter) == 2):
            time.append(float(splitter[0]))
            acf.append(float(splitter[1]))
        if (len(splitter) == 3):
            time.append(float(splitter[0]))
            acf.append(complex(float(splitter[1]), float(splitter[2])))
    
    myfile.close() 

    return time, acf

import os,sys
if __name__ == "__main__":
    (time, acf) = read_acf(sys.argv[1])
    for i in range(len(time)):
        print time[i], acf[i]
