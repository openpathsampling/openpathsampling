def set_features(*features):
    """
    Select snapshot features
    """
    def _decorator(cls):

        cls.__features__ = features
        cls._feature_attributes = list()

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

            # update __init__

        _super_init = cls.__init__

        def _init(self, **kwargs):

            _super_init(self, **kwargs)

            for attr in self._feature_attributes:
                if attr in kwargs:
                    setattr(self, attr, kwargs[attr])

        cls.__init__ = _init
        return cls

    return _decorator