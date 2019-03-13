from __future__ import print_function
import os
import pkg_resources
import tempfile
import subprocess
import openpathsampling.version

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clobber', action='store_true',
                    help='completely overwrite S3 bucket (default false)')
parser.add_argument('--dry', action='store_true')
opts = parser.parse_args()

CLOBBER = opts.clobber

BUCKET_NAME = 'openpathsampling.org'
if not openpathsampling.version.release:
    PREFIX = 'latest'
else:
    PREFIX = openpathsampling.version.short_version

def is_s3cmd_installed():
    dists = pkg_resources.working_set
    if not any(d.project_name == 's3cmd' for d in dists):
        raise ImportError('The s3cmd package is required. '
                          'try $ pip install s3cmd')

def run_cmd(template, config_filename, dry=False):
    cmd = template.format(
        config=config_filename,
        bucket=BUCKET_NAME,
        prefix=PREFIX)
    if dry:
        return_val = cmd
    else:
        is_s3cmd_installed()
        return_val = subprocess.call(cmd.split())
    return return_val


# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
with tempfile.NamedTemporaryFile('w') as f:
    f.write('''[default]
access_key = {AWS_ACCESS_KEY_ID}
secret_key = {AWS_SECRET_ACCESS_KEY}
'''.format(**os.environ))
    f.flush()

    s3cmd_core = "s3cmd --config {config} --no-mime-magic --guess-mime-type"

    if CLOBBER:
        s3cmd = s3cmd_core + " put --recursive "
    else:
        s3cmd = s3cmd_core + " sync "

    # get the relevant docs onto S3
    template = s3cmd + "docs/_build/html/ s3://{bucket}/{prefix}/"
    print(run_cmd(template, f.name, dry=opts.dry))
    #cmd = template.format(
            #config=f.name,
            #bucket=BUCKET_NAME,
            #prefix=PREFIX)
    #return_val = subprocess.call(cmd.split())

    # get the index file (main redirect)
    template = s3cmd + 'devtools/ci/index.html s3://{bucket}/'
    print(run_cmd(template, f.name, dry=opts.dry))
    #cmd = template.format(
            #config=f.name,
            #bucket=BUCKET_NAME)
    #return_val = subprocess.call(cmd.split())

