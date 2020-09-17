try:
    import configparser
except ImportError:
    import ConfigParser as configparser

import argparse
import os.path

def make_parser():
    default_setup = os.path.join(os.path.dirname(__file__), "..",
                                 "setup.cfg")
    parser = argparse.ArgumentParser()
    parser.add_argument('setup_cfg', nargs='?', default=default_setup)
    parser.add_argument('--extra', default=None)
    return parser

def clean(string):
    return " ".join(string.split("\n"))

def main(setup_cfg, extra):
    config = configparser.ConfigParser()
    config.read(setup_cfg)
    # replace get with the __getitem__ syntax when we drop py27
    if extra is None:
        # reqs = config['options']['install_requires']
        reqs = config.get('options', 'install_requires')
    else:
        # reqs = config['options.extras_require'][extra]
        reqs = config.get('options', 'extras_require', extra)
    return clean(reqs)

if __name__ == "__main__":
    argparser = make_parser()
    args = argparser.parse_args()
    print(main(args.setup_cfg, args.extra))

