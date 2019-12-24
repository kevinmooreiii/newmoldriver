""" new functions for filesystem stuff
"""


def orb_rest(spc_info, thy_info):
    """ append to the theory level the orb restricted stuff
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_level = thy_info[0:3]
    thy_level.append(orb_restr)
    return thy_level
