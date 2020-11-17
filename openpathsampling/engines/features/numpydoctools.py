class NumpyDocTools(object):
    """
    Helper class to combine several orthogonal numpydocs

    It can read numpy docstrings for multiple sources and combine these into one numpydoc
    string. This mainly used for th features mixin to make sure the docstring of the
    combined class is accurate.

    Attributes
    ----------
    sections : dict str to list of str
        a list of content per section name
    order : list of str
        the final use order of sections
    """

    known_sections = [
        'HEAD',
        'attributes',
        'parameters',
        'other parameters',
        'returns',
        'yields',
        'raises',
        'notes',
        'see also',
        'references',
    ]

    def __init__(self):
        self.sections = dict()
        self.order = None

    def clear(self):
        """
        Clear the parser and remove all previous loaded docstrings

        """
        self.sections = dict()

    @staticmethod
    def _remove_trailing_spaces(lines):
        # determine overall shift in spaces from the left

        if len(lines) == 0:
            return []

        min_space = 10000
        for line in lines:
            left = line.lstrip()
            if len(left) > 0:
                min_space = min(len(line) - len(left), min_space)

        lines = [line[min_space:] for line in lines]

        return lines

    def add_docs_from(self, obj, *args, **kwargs):
        """

        Parameters
        ----------
        obj : object
            the object from which the `__doctring__` property should be parsed
        section : str
            the optional name of the initial section
        keep_only : list of str
            a list of sections that should be kept. All other section are discarded
        translate : dict str: str
            a dictionary of from_setion: to_section type. The from_sections are interpreted
            as the to_section.
        """
        if obj.__doc__ is not None:
            self.add_docs(obj.__doc__, *args, **kwargs)

    def add_docs(self, docs, section=None, keep_only=None, translate=None):
        """

        Parameters
        ----------
        docs : str
            the docstring to be parsed
        section : str
            the optional name of the initial section
        keep_only : list of str
            a list of sections that should be kept. All other section are discarded
        translate : dict str: str
            a dictionary of from_setion: to_section type. The from_sections are interpreted
            as the to_section.
        """
        if section is None:
            current_section = 'HEAD'
        else:
            current_section = section

        new_section = None

        lines = docs.split('\n')

        sec = list()

        for raw_line in lines:
            line = raw_line.rstrip()
            lower = line.strip().lower()

            if lower in self.known_sections:
                new_section = lower
            elif new_section is not None and (line[-1] in ['-', '='] or len(line) == 0):
                self.add_block(current_section, sec, keep_only, translate)
                sec = list()

                current_section = new_section
                new_section = None
            else:
                new_section = None
                sec += [line]

        self.add_block(current_section, sec, keep_only, translate)

    @property
    def attributes(self):
        """
        List of all attributes found in the combined docstring

        Returns
        -------
        list of str
            the list of all attributes found

        """
        attr = []
        if 'attributes' in self.sections:
            for line in self.sections['attributes']:
                l = line.rstrip()
                if l and l[0:1] != ' ':
                    attr += [line.split()[0]]

        return attr

    def add_block(self, section, lines, keep_only=None, translate=None):
        """
        Add a section of a numpy docstring

        Parameters
        ----------
        section : str
            name of the section type, e.g. Attributes, Notes, ...
        lines : list of str
            the lines to be added
        keep_only : list of str
            a list of sections that should be kept. All other section are discarded
        translate : dict str: str
            a dictionary of from_setion: to_section type. The from_sections are interpreted
            as the to_section.

        """
        if not lines:
            return

        if keep_only is not None and type(keep_only) is str:
            keep_only = [keep_only]

        lines = self._remove_trailing_spaces(lines)

        while lines and lines[-1] == '':
            lines.pop()

        while lines and lines[0] == '':
            lines.pop(0)

        if translate is not None and section in translate:
            section = translate[section]

        if keep_only is None or section in keep_only:
            if section not in self.sections:
                self.sections[section] = list()

            self.sections[section] += lines + ['']

    @property
    def _docs(self):
        """
        Recreate the final list of docstrings in the set order

        Returns
        -------
        list of str
            the list of docstrings

        """
        if self.order is None:
            order = [section for section in self.known_sections if section in self.sections]
        else:
            order = self.order

        docs = list()

        for section in order:
            doc = self.sections[section]

            while doc and doc[-1] == '':
                doc.pop()

            while doc and doc[0] == '':
                doc.pop(0)

            if section != 'HEAD':
                cap = ' '.join([s[0].upper() + s[1:].lower() for s in section.split(' ')])
                doc = [
                    cap,
                    '-' * len(section)
                ] + doc

            doc += ['']

            docs += doc

        return docs

    def get_docstring(self):
        """
        Get the final docstring as a string

        Returns
        -------
        str
            the combined docstring
        """
        doc = '\n'.join(self._docs)

        doc = doc.replace('\n\n\n', '\n\n')

        return doc