def get_val(nested_dict, key):
    """
    Get the value to the key from the (possibly) nested dict values
    NOTE: this only works for the context we currently have, as it only descends down
    the first subdict in each dictionary recursively, so that it is possible that the key
    is in one of the dictionaiers, that gets never visited
    """
    if not key in nested_dict:
        for subdict in nested_dict.values():
            if not hasattr(subdict, '__iter__'):
                continue
            return get_val(subdict, key)
    else:
        return nested_dict[key]

    # the recursion has fallen through here, the key was nowhere to be found in
    # the nested dict (resp. not in the ones we looked at)
    return None
