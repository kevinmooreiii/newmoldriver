"""
Functions to check the filesystem
"""

from routines import util
from lib import msg


def check_save(save_fs, tsk, obj):
    """ check if information exists in a save directory for
        conformers, tau, scan
    """
    assert obj in ('conf', 'tau', 'scan')
    avail = True
    if obj == 'conf':
        save_locs = util.min_energy_conformer_locators(save_fs)
    else:
        save_locs = [save_fs.leaf.existing()[0]]
    if not save_locs:
        msg.ini_info_noavail_msg(tsk)
        avail = False
    elif not save_fs.leaf.file.geometry.exists(save_locs):
        msg.ini_info_noavail_msg(tsk)
        avail = False
    return avail
