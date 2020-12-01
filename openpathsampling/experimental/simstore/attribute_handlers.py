from .my_types import parse_ndarray_type


class HandlerFactory(object):
    def is_my_type(self, type_str):
        pass

    def serializer(self, type_str):
        pass

    def deserializer(self, type_str):
        pass

class NDArrayHandlerFactory(HandlerFactory):
    def is_my_type(self, type_str):
        as_ndarray = parse_ndarray_type(type_str)

    def serializer(self, type_str):
        as_ndarray = self.is_my_type(type_str)
        if as_ndarray:
            dtype, shape = as_ndarray
            handler = lambda data, _: \
                    np.fromstring(data, dtype=dtype).reshape(shape)
            return handler

    def deserializer(self, type_str):
        as_ndarray = self.is_my_type(type_str)
        if as_ndarray:
            dtype, shape = as_ndarray
            return lambda arr: \
                    arr.astype(dtype=dtype, copy=False).tostring()



