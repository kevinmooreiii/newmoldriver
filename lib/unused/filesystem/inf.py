"""
Build paths and file systesm given species and theory
'info' objects
"""

import autofile
import moldr


def get_spc_run_path(run_prefix, spc_info):
    """ create species run path
    """
    spc_run_fs = autofile.fs.species(run_prefix)
    spc_run_fs.leaf.create(spc_info)
    spc_run_path = spc_run_fs.leaf.path(spc_info)
    return spc_run_path


def get_spc_save_path(save_prefix, spc_info):
    """ create species save path
    """
    spc_save_fs = autofile.fs.species(save_prefix)
    spc_save_fs.leaf.create(spc_info)
    spc_save_path = spc_save_fs.leaf.path(spc_info)
    return spc_save_path


def get_thy_run_path(run_prefix, spc_info, thy_info):
    """ create theory run path
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_run_path = get_spc_run_path(run_prefix, spc_info)
    thy_run_fs = autofile.fs.theory(spc_run_path)
    thy_run_fs.leaf.create(thy_lvl)
    thy_run_path = thy_run_fs.leaf.path(thy_lvl)
    return thy_run_path


def get_thy_save_fs(save_prefix, spc_info, thy_info):
    """ create theory save filesystem
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    return thy_save_fs, thy_lvl


def get_thy_save_path(save_prefix, spc_info, thy_info):
    """ create theory save path
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_info)
    thy_lvl = thy_info[0:3]
    thy_lvl.append(orb_restr)
    spc_save_path = get_spc_save_path(save_prefix, spc_info)
    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs.leaf.create(thy_lvl)
    thy_save_path = thy_save_fs.leaf.path(thy_lvl)
    return thy_save_path


def get_thy_info(lvldic):
    """ convert theory level dictionary to theory information array
    """
    err_msg = ''
    info = ['program', 'method', 'basis', 'orb_res']
    for i, inf in enumerate(info):
        if inf in lvldic:
            info[i] = lvldic[inf]
        else:
            err_msg = inf
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return info


def get_spc_info(spc_dct_i):
    """ convert species dictionary to species_info array
    """
    err_msg = ''
    props = ['ich', 'chg', 'mul']
    for i, prop in enumerate(props):
        if prop in spc_dct_i:
            props[i] = spc_dct_i[prop]
        else:
            err_msg = prop
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return props
