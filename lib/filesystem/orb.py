""" new functions for filesystem stuff
"""


def orbital_restriction(spc_info, thy_level):
    """ orbital restriction logical
    """
    mul = spc_info[2]
    if thy_level[3] == 'RR':
        orb_restr = True
    elif thy_level[3] == 'UU':
        orb_restr = False
    elif thy_level[3] == 'RU':
        orb_restr = bool(mul == 1)
    return orb_restr


def orb_rest(spc_info, thy_info):
    """ append to the theory level the orb restricted stuff
    """
    orb_restr = orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)
    return thy_level
