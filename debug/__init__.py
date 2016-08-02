class DebugDict(dict):
    """Debug a dictionary to make sure only specified keys are being accessed"""

    def __init__(self, copy=None, safe_keys=None):
        super(DebugDict, self).__init__()
        if isinstance(copy, dict):
            for (k, v) in copy.iteritems():
                self[k] = v
        elif copy is not None:
            raise AttributeError('Cannot invoke copy constructor on argument')
        if isinstance(safe_keys, dict):
            self.safe_keys = safe_keys.copy()
        elif safe_keys is None:
            self.safe_keys = {}
        else:
            raise KeyError('Unknown type for safe keys')
        for (k, v) in safe_keys.iteritems():
            if isinstance(v, dict):
                if isinstance(self[k], dict):
                    self[k] = DebugDict(self[k], v)
                else:
                    raise KeyError("Cannot recursively make DebugDict for attribute {}".format(k))

    def get(self, item, **kwargs):
        if item not in self.safe_keys:
            raise KeyError("Unsafe key: {}".format(item))
        return super(DebugDict, self).get(item, **kwargs)

    def __getitem__(self, item):
        if item not in self.safe_keys:
            raise KeyError("Unsafe key: {}".format(item))
        return super(DebugDict, self).__getitem__(item)
