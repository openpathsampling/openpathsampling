# This file is vendored from Autorelease
import os
import ast
import sys

import fnmatch  # Py 2

from setuptools import setup

def _glob_glob_recursive(directory, pattern):
    # python 2 glob.glob doesn't have a recursive keyword
    # this implements for the specific case that we want an exact match
    # See also https://stackoverflow.com/a/2186565
    matches = []
    for root, dirname, filenames in os.walk(directory):
        matches.extend([os.path.join(root, filename)
                        for filename in fnmatch.filter(filenames, pattern)])
    return matches


class VersionPyFinder(object):
    _VERSION_PY_FUNCTIONS = ['get_git_version', 'get_setup_cfg']
    def __init__(self, filename='version.py', max_depth=2):
        self.filename_base = filename
        self.max_depth = max_depth
        self.depth = None
        self.filename = os.getenv("AUTORELEASE_VERSION_PY",
                                  self._first_eligible())
        self.functions = self._get_functions(self.filename)

    def _find_files(self):
        # all_files = glob.glob("**/" + self.filename_base, recursive=True)
        all_files = _glob_glob_recursive('.', self.filename_base)
        meets_depth = [fname for fname in all_files
                       if len(fname.split(os.sep)) <= self.max_depth + 1]
        return meets_depth

    def _is_eligible(self, filename):
        with open(filename, mode='r') as f:
            contents = f.read()

        tree = ast.parse(contents)
        # we requrie that our functions be defined at module level -- we
        # know that's how we wrote them, at least!
        all_functions = [node.name for node in tree.body
                         if isinstance(node, ast.FunctionDef)]
        return all(func in all_functions
                   for func in self._VERSION_PY_FUNCTIONS)

    def _first_eligible(self):
        all_files = self._find_files()
        for fname in all_files:
            if self._is_eligible(fname):
                return fname
        return None

    @property
    def version_setup_depth(self):
        def get_depth(fname):
            return len(os.path.abspath(fname).split(os.sep))

        # we assume thta setup.py is in the same dir as setup.cfg
        diff = get_depth(self.filename) - get_depth(__file__)
        return diff

    def _get_functions(self, filename):
        with open(self.filename, mode='r') as f:
            contents = f.read()

        tree = ast.parse(contents)

        class MakeImportError(ast.NodeTransformer):
            """converts a from x import y into an import error"""
            def __init__(self, import_name):
                self.import_name = import_name

            def visit_ImportFrom(self, node):
                if node.module == self.import_name:
                    replacement = ast.Raise(exc=ast.Call(
                        func=ast.Name(id='ImportError', ctx=ast.Load()),
                        args=[],
                        keywords=[],
                    ), cause=None)
                    return ast.copy_location(replacement, node)
                else:
                    return node

        import_remover = MakeImportError("_installed_version")
        tree = import_remover.visit(tree)
        ast.fix_missing_locations(tree)

        locs = dict(globals())
        exec(compile(tree, filename="version.py", mode='exec'), locs)
        return {f: locs[f] for f in self._VERSION_PY_FUNCTIONS}


def write_installed_version_py(filename="_installed_version.py",
                               src_dir=None):
    version_finder = VersionPyFinder()
    directory = os.path.dirname(version_finder.filename)
    depth = version_finder.version_setup_depth
    get_git_version = version_finder.functions['get_git_version']
    get_setup_cfg = version_finder.functions['get_setup_cfg']

    installed_version = os.path.join(directory, "_installed_version.py")
    content = "_installed_version = '{vers}'\n"
    content += "_installed_git_hash = '{git}'\n"
    content += "_version_setup_depth = {depth}\n"

    # question: if I use the __file__ attribute in something I compile from
    # here, what is the file?
    my_dir = os.path.abspath(os.path.dirname(__file__))
    conf = get_setup_cfg(directory=my_dir, filename='setup.cfg')
    # conf = get_setup_cfg(directory=my_dir, filename='new_setup.cfg')
    version = conf.get('metadata', 'version')
    git_rev = get_git_version()

    if src_dir is None:
        src_dir = conf.get('metadata', 'name')

    with open (os.path.join(src_dir, filename), 'w') as f:
        f.write(content.format(vers=version, git=git_rev, depth=depth))

if __name__ == "__main__":
    # TODO: only write version.py under special circumstances
    write_installed_version_py()
    # write_version_py(os.path.join('autorelease', 'version.py'))
    setup()

