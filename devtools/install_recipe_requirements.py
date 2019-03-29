#!/usr/bin/env python
from __future__ import print_function

import sys
import tempfile
import argparse

# This requires that you have already installed conda and pyyaml
import conda.cli
import yaml

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('recipe_file', default=None, nargs='?',
                        help="recipe file (use stdin if not given)")
    parser.add_argument('--dry', action='store_true', default=False,
                        help="don't install; print to stdout instead")
    parser.add_argument('--individual', action='store_true', default=False,
                        help="install required packages one at a time")
    parser.add_argument('--env', default=None, type=str,
                        help="if given, an environment file will be created")
    return parser.parse_args()

def read_input(filename):
    """Read data from filename if given, stdin otherwise."""
    if filename is None:
        file_data = sys.stdin.read()
    else:
        with open(filename, 'r') as f:
            file_data = f.read()
    return file_data

def recipe_to_requirements(recipe_yaml):
    reqs = yaml.load(recipe_yaml)['requirements']
    required_packages = set(reqs['build'] + reqs['run'])
    return required_packages

if __name__ == "__main__":
    # convert the conda recipe to file of dependencies
    args = parse_arguments()
    recipe_yaml = read_input(args.recipe_file)
    required_packages = recipe_to_requirements(recipe_yaml)
    req_file_str = "\n".join(required_packages)

    if args.dry:
        print(req_file_str)
    else:
        for install in required_packages:
            channels = ['conda-forge', 'omnia']
            channel_str = sum([['-c', channel] for channel in channels], [])
            sys.argv = ['conda', 'install'] + channel_str + [install]
            # print(sys.argv)
            conda.cli.main()

        # with tempfile.NamedTemporaryFile() as tmp:
            # tmp.write(req_file_str)
            # print(tmp.name)
            # sys.argv = ['conda', 'install', '--file', tmp.name]
            # conda.cli.main('conda', 'install', '--file', tmp.name)
