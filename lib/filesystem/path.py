"""
Build paths for spc and thy file systems
"""

import autofile
from lib import moldr


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


def get_rxn_fs(run_prefix, save_prefix, ts):
    """ get filesystems for a reaction
    """
    rxn_ichs = ts['rxn_ichs']
    rxn_chgs = ts['rxn_chgs']
    rxn_muls = ts['rxn_muls']
    ts_mul = ts['mul']

    rxn_run_fs = autofile.fs.reaction(run_prefix)
    rxn_run_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_run_path = rxn_run_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    rxn_save_fs = autofile.fs.reaction(save_prefix)
    rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
    rxn_save_path = rxn_save_fs.leaf.path(
        [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    return rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path
