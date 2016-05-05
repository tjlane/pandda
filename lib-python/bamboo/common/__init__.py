
class Meta(object):
    """Object for storing random data - can be edited and added to as required"""
    _initialized = False
    def __init__(self, args=None):
        if isinstance(args, dict):
            for k in args:
                self.__dict__[k] = args[k]
        elif isinstance(args, list):
            for l in args:
                self.__dict__[l] = None
        elif args is not None:
            raise Exception('args must be dict or list')
        _initialized = True

    def summary(self):
        out = []
        for k in self.__dict__:
            out.append('{!s}: {!s}'.format(k, self.__dict__[k]))
        return '\n'.join(out)

class Info(Meta):
    """Same as Meta, but cannot change variables after setting _initialized=True"""
    def __setattr__(self, name, value):
        if self._initialized and (not hasattr(self, name)):
            raise AttributeError('Cannot set new attributes after initialisation: {!s}'.format(name))
        object.__setattr__(self, name, value)

def compare_dictionaries(dict1, dict2):
    """Check differences between dicts. Returns tuple(differences in the keys, and differences in the values)."""

    uniq_keys = list(set(dict1.keys()).symmetric_difference(set(dict2.keys())))
    uniq_summ = [(d_key, (d_key in dict1.keys()), (d_key in dict2.keys())) for d_key in uniq_keys]

    comm_keys = list(set(dict1.keys()).intersection(set(dict2.keys())))
    comm_summ = [(c_key, dict1[c_key], dict2[c_key]) for c_key in comm_keys if dict1[c_key] != dict2[c_key]]

    return uniq_summ, comm_summ

