def tuple_keys_to_dict(to_dict, attr_name):
    """to-dict for attributes that are dicts with tuples as keys"""
    def tuple_keys_dict_to_dict(name, tuple_keys_dict):
        keys = list(tuple_keys_dict.keys())
        values = list(tuple_keys_dict.values())
        dct = {name + "_tuple_keys": keys,
               name + "_values": values}
        return dct

    def inner(self):
        dct = to_dict(self)
        dct[attr_name] = tuple_keys_dict_to_dict(attr_name, dct[attr_name])
        return dct
    return inner

def tuple_keys_from_dict(from_dict, attr_name):
    """from-dict for attributes that are dicts with tuples as keys"""
    def tuple_keys_dict_from_dict(name, dct):
        keys = dct[name + "_tuple_keys"]
        values = dct[name + "_values"]
        return {tuple(k): v for k, v in zip(keys, values)}

    def inner(cls, dct):
        dct = dict(dct)  # copy
        dct[attr_name] = tuple_keys_dict_from_dict(attr_name, dct[attr_name])
        return from_dict(dct)
    return inner
