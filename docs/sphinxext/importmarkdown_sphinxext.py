# Copied from the yt_project, commit e8fb57e
# yt/doc/extensions/notebook_sphinxext.py
#  https://bitbucket.org/yt_analysis/yt/src/e8fb57e66ca42e26052dadf054a5c782740abec9/doc/extensions/notebook_sphinxext.py?at=yt
#
# I want to rewrite it do the same for Markdown files.
# create a directive ::markdown [file.md] that add

from __future__ import print_function

import os, shutil
from sphinx.util.compat import Directive
from docutils.parsers.rst import directives

from docutils.frontend import OptionParser
from docutils.utils import new_document
from docutils.parsers.rst import Parser
import subprocess

class NotebookDirective(Directive):
    """Insert a markdown (.md) file converted to .rst into the document

    This uses pandoc and does no testing if it is installed. So make sure of
    it or disable the extension in the conf.py file

    Use this using `.. nodebook:: [markdownfile.md]`
    """
    required_arguments = 1
    optional_arguments = 1
    option_spec = {'skip_exceptions' : directives.flag}
    final_argument_whitespace = True

    def run(self): # check if there are spaces in the notebook name
        md_path = self.arguments[0]
        if ' ' in md_path: raise ValueError(
            "Due to issues with docutils stripping spaces from links, white "
            "space is not allowed in notebook filenames '{0}'".format(md_path))

        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        # get path to markdown file
        md_filename = self.arguments[0]
        md_dir = os.path.join(setup.confdir, '..')
        md_abs_path = os.path.abspath(os.path.join(md_dir, md_filename))

        p = subprocess.Popen(['pandoc',
                                '--from=markdown', '--to=rst', md_abs_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        out, err = p.communicate()

        print(err)

        settings = OptionParser(components=(Parser,)).get_default_values()
        parser = Parser()
        document = new_document('DOC', settings)
        parser.parse(out, document)

        return [node for node in document]

def formatted_link(path):
    return "`%s <%s>`__" % (os.path.basename(path), path)

def visit_notebook_node(self, node):
    self.visit_raw(node)

def depart_notebook_node(self, node):
    self.depart_raw(node)

def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

#    app.add_stylesheet('css/md.css')

    app.add_directive('markdown', NotebookDirective)

