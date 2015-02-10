#### Building the documentation

run
	
	make html
	
to generate the new documentation using sphinx.
For this project we agreed on the numpy style docstrings. See e.g.

http://pyfact.readthedocs.org/en/latest/code-doc_np_HOWTO_DOCUMENT.html
	
for a guide.

*Requires*: 

* [Sphinx](http://www.sphinx-doc.org/)
* [sphinx_rtd_theme](https://github.com/snide/sphinx_rtd_theme) : `pip install sphinx_rtd_theme`
* [numpydoc](https://pypi.python.org/pypi/numpydoc) : `pip install numpydoc`