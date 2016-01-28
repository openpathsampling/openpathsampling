from openpathsampling.netcdfplus import DelayedLoader
from numpydocparse import NumpyDocParser

def has(attr):
    def _has(func):
        def inner(self, *args, **kwargs):
            if hasattr(self, attr) and getattr(self, attr) is not None:
                return func(self, *args, **kwargs)
            else:
                return None

        return inner

    return _has


def set_features(*features):
    """
    Select snapshot features
    """

    parser = NumpyDocParser()
    use_lazy_reversed = True

    def _decorator(cls):

        # important for compile to work properly
        import openpathsampling as paths

        parser.clear()

        __features__ = dict()
        __features__['classes'] = features
        for name in ['attributes', 'minus', 'reversal', 'properties', 'flip']:
            __features__[name] = list()

        if use_lazy_reversed:
            __features__['lazy'] = ['_reversed']
        else:
            __features__['lazy'] = list()

        for feature in features:
            # loop over all the features

            # add properties to class
            for prop in feature.attributes:
                if hasattr(feature, prop) and callable(getattr(feature, prop)):
                    __features__['properties'] += [prop]
                    setattr(cls, prop, property(getattr(feature, prop)))

            for name in ['attributes', 'minus', 'lazy', 'flip']:
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

        # add lazy decorators

        for attr in __features__['lazy']:
            setattr(cls, attr, DelayedLoader())

        cls.__features__ = __features__

        # update docstring

        # from class
        parser.add_docs_from(cls)

        # from top of features
        for feat in __features__['classes']:
            parser.add_docs_from(feat)

            # from properties ???
            for prop in __features__['properties']:
                if hasattr(feat, prop):
                    if prop not in parser.attributes:
                        parser.add_docs_from(
                            getattr(feat, prop),
                            keep_only=['attributes'],
                            translate={'returns': 'attributes'}
                        )

        cls.__doc__ = parser.get_docstring()

        # print '+++++++++++++++++++++++++++'
        # print cls.__doc__
        # print '+++++++++++++++++++++++++++'
        # print

        # compile copy()
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

        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.copy = copy

        except RuntimeError as e:
            print e
            pass


        # compile create_reversed
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

        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.create_reversed = create_reversed

        except RuntimeError as e:
            print e
            pass

        # compile __init__

        # we use as signature all attributes
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

        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.__init__ = __init__

        except RuntimeError as e:
            print e
            pass

        return cls

    return _decorator