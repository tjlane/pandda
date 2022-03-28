import giant.logs as lg
logger = lg.getLogger(__name__)

import json, collections

def make_sorted_dict(dictionary, recursive=False):

    new_dict = collections.OrderedDict()

    for k, v in sorted(dictionary.items()):

        # If is dict, make sorted dict
        if (recursive is True) and hasattr(v, 'keys'):
            v = make_sorted_dict(v, recursive=recursive)
            
        new_dict[k] = v

    return new_dict

def show_dict(dictionary, logger=None):

    if logger is None: 
        logger = lg.getLogger(__name__)
        
    logger(
        json.dumps(
            dictionary,
            indent=2,
            )
        )

def merge_dicts(
    master_dict,
    merge_dict,
    dict_class = collections.OrderedDict,
    overwrite = False,
    ):

    for j_key, j_val in merge_dict.items():

        if hasattr(j_val, 'keys'):

            m_dict = master_dict.setdefault(
                j_key, 
                dict_class(),
                )

            merge_dicts(
                master_dict = m_dict,
                merge_dict = j_val,
                dict_class = dict_class,
                overwrite = overwrite,
                )
        else: 

            if (overwrite is False) and (j_key in master_dict):
                show_dict(master_dict)
                show_dict(merge_dict)
                raise ValueError('error joining "{}"'.format(j_key))

            master_dict[j_key] = j_val

    return master_dict

def pretty_format_list(values, key=None, line_width=80):

    s_ = ""

    if (key is not None): 
        key_ = '{k} : '.format(k=str(key))
        s_ += key_

    s_ += '[\n'

    ll = 0
    for v in values: 
        v_ = (str(v) + ', ')
        l = len(v_)
        if (ll>0) and (ll+l > line_width): # would go over line width
            s_ += '\n'
            ll = 0
        s_ += v_
        ll += l

    s_ += '\n]'

    return s_.strip().strip(',')
