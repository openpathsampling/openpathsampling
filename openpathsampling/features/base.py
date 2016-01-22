from openpathsampling.netcdfplus import DelayedLoader

def set_features(*features):
    """
    Select snapshot features
    """

    use_lazy_reversed = True

    def _decorator(cls):

        # important for compile to work properly
        import openpathsampling as paths

        __features__ = dict()
        __features__['classes'] = features
        for name in ['attributes', 'minus', 'reversal', 'properties']:
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

            for name in ['attributes', 'minus', 'lazy']:
                if hasattr(feature, name):
                    content = getattr(feature, name)
                    if type(content) is str:
                        content = [content]

                    __features__[name] += content
                    
        __features__['reversal'] = [
            attr for attr in __features__['attributes']
            if attr not in __features__['minus']
            and attr not in __features__['properties']
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

        docs = cls.__doc__.split('\n')

        docs += [
            '',
            'Attributes',
            '----------'
        ]

        for feat in __features__:
            if feat.__doc__ is not None:
                docs += [feat.__doc__]

        cls.__doc__ = '\n'.join(map(lambda x : x.strip(), docs))
        cls.__doc__ = cls.__doc__.replace('\n\n\n', '\n\n')

        # compile create_reversed
        code = []
        code += [
            "def create_reversed(self):",
            "    this = cls.__new__(cls)",
            "    this._is_reversed = True",
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
        parameters += map('{0}=None'.format, __features__['parameters'])
        signature = ', '.join(parameters)
        code = []
        code += [
            "def __init__(self, %s):" % signature,
            "    self._is_reversed = False",
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