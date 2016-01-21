def set_features(*features):
    """
    Select snapshot features
    """
    def _decorator(cls):

        import openpathsampling as paths

        cls.__features__ = features
        cls._feature_attributes = list()
        cls._feature_attributes_minus = list()
        cls._feature_attributes_not = list()

        for feature in features:
            # loop over all the features

            # add properties to class
            for prop in feature.properties:
                if hasattr(feature, prop):
                    setattr(cls, prop, property(getattr(feature, prop)))

            # update list of all attributes to be stored
            if type(feature.attributes) is str:
                cls._feature_attributes += [feature.attributes]
            else:
                cls._feature_attributes += feature.attributes

            # update list of all attributes to multiplied with -1
            if type(feature.attributes_minus) is str:
                cls._feature_attributes_minus += [feature.attributes_minus]
            else:
                cls._feature_attributes_minus += feature.attributes_minus

            # update list of all attributes to be flipped
            if type(feature.attributes_not) is str:
                cls._feature_attributes_not += [feature.attributes_not]
            else:
                cls._feature_attributes_not += feature.attributes_not

        use_lazy_reversed = True

        # update docstring

        docs = cls.__doc__.split('\n')

        docs += [
            '',
            'Attributes',
            '----------'
        ]

        for feat in cls.__features__:
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
            "    this.engine = self.engine",
        ]

        if use_lazy_reversed:
            code += [
                "    this._lazy = {cls._reversed : self}"
            ]
        else:
            code += [
                "    this._reversed = self"
            ]

        code += map("    this.{0} = self.{0}".format,
                    [attr for attr in cls._feature_attributes
                     if attr not in cls._feature_attributes_minus and
                     attr not in cls._feature_attributes_not
                     ]
                    )
        code += map("    this.{0} = - self.{0}".format, cls._feature_attributes_minus)
        code += map("    this.{0} = ~ self.{0}".format, cls._feature_attributes_not)

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
        parameters += map('{0}=None'.format, cls._feature_attributes)
        signature = ', '.join(parameters)
        code = []
        code += [
            "def __init__(self, %s):" % signature,
            "    self._is_reversed = False",
            "    self.engine = engine",
        ]

        if use_lazy_reversed:
            code += [
                "    self._lazy = {cls._reversed : None}"
            ]
        else:
            code += [
                "    self._reversed = None"
            ]

        code += map("    self.{0} = {0}".format, cls._feature_attributes)

        try:
            cc = compile('\n'.join(code), '<string>', 'exec')
            exec cc in locals()

            cls.__init__ = __init__

        except RuntimeError as e:
            print e
            pass

        return cls

    return _decorator