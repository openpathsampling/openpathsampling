from openpathsampling.netcdfplus import DelayedLoader
from numpydoctools import NumpyDocTools
# from openpathsampling.snapshot import FeatureSnapshot


def set_features(*features):
    """
    Select snapshot features
    """

    # create a parser that can combine numpy docstrings
    parser = NumpyDocTools()

    # if this is set to true than the reversed counterpart of a Snapshot
    # is saved a lazy pointer otherwise we create a full copy
    use_lazy_reversed = True

    def _decorator(cls):
        """
        Class decorator that will attach function for compiled features

        This function will use a list of features and create `__init__` and
        copy functions based on the structure of features for performance.

        It will also take care of creating a joined docstring and the corect
        signature of the `__init__` function

        A attribute `__features__` will be added that contains information
        about the used features their structure. It is a dictionary with the
        following keys

        classes : dict of used features
        lazy : dict of string
            names of features that are treated as lazy loaded object
        properties : dict of string
            names of features that are treated as properties
        lazy : dict of string
            names of features that are treated as lazy loaded object
        reversal : dict of string
            names of features that are treated as being reversible
        minus : dict of string
            names of features that are treated as being reversible and
            should be multiplied by -1.0
        flip : dict of string
            names of features that are treated as being reversible and
            should be negated `~`
        attributes : dict of string
            names of features that are treated as being class attributes
        parameters : dict of string
            names of features that attributes but not properties and hence
            possible parameters for creation
        numpy : dict of string
            names of features that can use numpy for faster copying, etc.

        Parameters
        ----------
        cls : the `class` to the modified

        Returns
        -------
        class
            the modified class

        """

        # important for compile to work properly
        import openpathsampling as paths

        parser.clear()

        # create and fill `__features__` with values from feature structures
        __features__ = dict()
        __features__['classes'] = features
        for name in ['attributes', 'minus', 'reversal', 'properties', 'flip', 'numpy', 'lazy']:
            __features__[name] = []

        if use_lazy_reversed:
            __features__['lazy'] = ['_reversed']

        # loop over all the features
        for feature in features:

            # add properties
            for prop in feature.attributes:
                if hasattr(feature, prop) and callable(getattr(feature, prop)):
                    __features__['properties'] += [prop]
                    setattr(cls, prop, property(getattr(feature, prop)))

            # copy specific attribtue types
            for name in ['attributes', 'minus', 'lazy', 'flip', 'numpy']:
                if hasattr(feature, name):
                    content = getattr(feature, name)
                    if type(content) is str:
                        content = [content]

                    __features__[name] += content
                    
        __features__['reversal'] = [
            attr for attr in __features__['attributes']
            if attr not in __features__['minus']
            and attr not in __features__['properties']
            and attr not in __features__['flip']
        ]

        __features__['parameters'] = [
            attr for attr in __features__['attributes']
            if attr not in __features__['properties']
        ]

        # add descriptors that can handle lazy loaded objects
        for attr in __features__['lazy']:
            setattr(cls, attr, DelayedLoader())

        cls.__features__ = __features__

        # update the docstring to be a union of docstrings from the class
        # and the features

        # get docstring from class
        parser.add_docs_from(cls)

        # from top of features
        for feat in __features__['classes']:
            parser.add_docs_from(feat)

            # from properties
            for prop in __features__['properties']:
                if hasattr(feat, prop):
                    if prop not in parser.attributes:
                        parser.add_docs_from(
                            getattr(feat, prop),
                            keep_only=['attributes'],
                            translate={'returns': 'attributes'}
                        )

        # set new docstring. This is only possible since our class is created
        # using a Metaclass for abstract classes `abc`. Normal classes cannot
        # have thier docstring changed.
        cls.__doc__ = parser.get_docstring()

        # compile the function for .copy()

        # def copy(self):
        #     this = cls.__new__(cls)
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1

        code = []
        code += [
            "def copy(self):",
            "    this = cls.__new__(cls)",
        ]

        if __features__['lazy']:
            code += [
                "    this._lazy = {",
            ]
            code += [
                "       cls.{0} : self._lazy[cls.{0}],".format(lazy)
                for lazy in __features__['lazy']
            ]
            code += [
                "    }"
            ]

        code += map(
            "    this.{0} = self.{0}".format,
            filter(
                lambda x : x not in __features__['lazy'],
                __features__['parameters']
            )
        )

        code += [
            "    return this"
        ]

        # compile the code and register the new function
        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.copy = copy

        except RuntimeError as e:
            print e
            pass

        # compile the function for .copyto(target)

        # def copyto(self, target):
        #     this = target
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1
        #     return this

        code = []
        code += [
            "def copyto(self, target):",
            "    this = target",
        ]

        if __features__['lazy']:
            code += [
                "    this._lazy = {",
            ]
            code += [
                "       cls.{0} : self._lazy[cls.{0}],".format(lazy)
                for lazy in __features__['lazy']
            ]
            code += [
                "    }"
            ]

        code += map(
            "    this.{0} = self.{0}".format,
            filter(
                lambda x : x not in __features__['lazy'],
                __features__['parameters']
            )
        )

        code += [
            "    return this"
        ]

        # compile the code and register the new function
        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.copyto = copyto

        except RuntimeError as e:
            print e
            pass

        # compile the function for .copyto(target)

        # def create_reversed(self):
        #     this = cls.__new__(cls)
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1
        #     this.feature2 = - self.feature2  # minus feature
        #     this.feature3 = ~ self.feature3  # flip features
        #     return this

        code = []
        code += [
            "def create_reversed(self):",
            "    this = cls.__new__(cls)",
        ]

        if __features__['lazy']:
            if use_lazy_reversed:
                code += [
                    "    this._lazy = {cls._reversed : self}"
                ]
            else:
                code += [
                    "    this._lazy = dict()"
                ]

        if not use_lazy_reversed:
            code += [
                "    this._reversed = self"
            ]

        code += map("    this.{0} = self.{0}".format, __features__['reversal'])
        code += map("    this.{0} = - self.{0}".format, __features__['minus'])
        code += map("    this.{0} = ~ self.{0}".format, __features__['flip'])

        code += [
            "    return this"
        ]

        # compile the code and register the new function
        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.create_reversed = create_reversed

        except RuntimeError as e:
            print e
            pass

        # compile the function for __init__

        # def __init__(self, attribute1=None, ... ):
        #     self._lazy = { ... }
        #     self.feature1 = self.feature1
        #     return this

        # we use as signature all feature names in parameters
        parameters = ['engine=None']
        for feat in __features__['parameters']:
            if feat in __features__['flip']:
                parameters += ['{0}=False'.format(feat)]
            else:
                parameters += ['{0}=None'.format(feat)]

        signature = ', '.join(parameters)
        code = []
        code += [
            "def __init__(self, %s):" % signature,
        ]

        if __features__['lazy']:
            if use_lazy_reversed:
                code += [
                    "    self._lazy = {cls._reversed : None}"
                ]
            else:
                code += [
                    "    self._lazy = dict()"
                ]

        if not use_lazy_reversed:
            code += [
                "    self._reversed = None"
            ]

        code += map("    self.{0} = {0}".format, __features__['parameters'])

        # compile the code and register the new function
        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.__init__ = __init__

        except RuntimeError as e:
            print e
            pass

        return cls

    return _decorator