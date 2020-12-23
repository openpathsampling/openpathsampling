import pathlib
import yaml
import argparse
import shutil
import urllib.request
import tempfile

def parse_conf():
    parser = argparse.ArgumentParser()
    parser.add_argument('conf')
    parser.add_argument('--copy', action='store_true')
    parser.add_argument('--delete', action='store_true')
    opts = parser.parse_args()
    if opts.copy == opts.delete:
        raise RuntimeError("Must set exactly one of --copy or --delete")
    action = 'copy' if opts.copy else 'delete'
    with open(opts.conf, mode='r') as f:
        conf = yaml.load(f.read(), yaml.SafeLoader)
    return conf, action

def copy_local(source, dest, files):
    source = pathlib.Path(source)
    dest = pathlib.Path(dest)
    for f in files:
        print("cp " + str(source / f) + " " + str(dest / f))
        shutil.copy(source / f, dest / f)

def copy_gitlab(source, dest, files):
    middle_str = "/-/raw/master/"
    dest = pathlib.Path(dest)
    for f in files:
        url = source + middle_str + f
        print(url)
        urllib.request.urlretrieve(url, dest / f)

DISPATCH = {
    'local': copy_local,
    'gitlab': copy_gitlab
}

def make_copies(fileset):
    dispatch_keys = set(DISPATCH)
    set_keys = set(fileset)
    set_types = list(dispatch_keys & set_keys)
    if len(set_types) != 1:
        raise RuntimeError("Bad yaml")
    dirs = fileset.get('mkdirs', [])
    for d in dirs:
        pathlib.Path(d).mkdir(parents=True)
    source = set_types[0]
    DISPATCH[source](fileset[source], fileset['target'], fileset['files'])

def delete_files(fileset):
    dest = pathlib.Path(fileset['target'])
    for f in fileset['files']:
        print(dest / f)
        (dest /f).unlink()
    for d in reversed(fileset['mkdirs']):
        print(d)
        try:
            pathlib.Path(d).rmdir()
        except OSError:
            pass  # guess we can't delete this one yet


def main():
    conf, runtype = parse_conf()
    func = {'copy': make_copies, 'delete': delete_files}[runtype]
    for fileset in conf['copyfiles']:
        func(fileset)


if __name__ == "__main__":
    main()





