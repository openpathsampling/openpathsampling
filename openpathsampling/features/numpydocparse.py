class NumpyDocParser(object):

    known_sections = [
        'Parameters',
        'Attributes',
        'See Also',
        'Raises',
    ]

    def __init__(self):
        self.sections = {}

    def add_docs(self, docs):
        lines = docs.split()
        if

    def