def setstorage(cls):
    def decorate(func):
        def call(this, storage = None, *args, **kwargs):
            if hasattr(kwargs, 'storage') and getattr(kwargs, 'storage') is not None:
                assert(isinstance(kwargs['storage'], Storage))
                old_storage = this.storage
                this.storage = kwargs['storage']

            result = func(*args, **kwargs)

            if hasattr(kwargs, 'storage') and getattr(kwargs, 'storage') is not None:
                this.storage = old_storage
            return result
        return call
    return decorate

def defaultstorage(cls):
    def decorate(func):
        def call(*args, **kwargs):
            this = args[0]
            if not hasattr(kwargs, 'storage') or getattr(kwargs, 'storage') is None:
                kwargs['storage'] = this.default_storage

            result = func(*args, **kwargs)
            return result
        return call
    return decorate
