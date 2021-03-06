"""
Functions to check the filesystem
"""

from lib.filesystem import minc as fsmin
from lib import printmsg


def check_save(save_fs, tsk, obj):
    """ check if information exists in a save directory for
        conformers, tau, scan
    """
    print('save fs leaf exist')
    print(save_fs.leaf.existing())
    print(save_fs.leaf.existing()[0])
    print([save_fs.leaf.existing()[0]])
    assert obj in ('conf', 'tau', 'scan')
    avail = True
    if obj == 'conf':
        save_locs = fsmin.min_energy_conformer_locators(save_fs)
    # Doesn't work well for tau since there are multipole tau we wish to comp
    elif obj == 'tau':
        save_locs = save_fs.leaf.existing()[0]
    else:
        save_locs = [save_fs.leaf.existing()[0]]
    if not save_locs:
        printmsg.ini_info_noavail_msg(tsk)
        avail = False
    elif not save_fs.leaf.file.geometry.exists(save_locs):
        printmsg.ini_info_noavail_msg(tsk)
        avail = False
    return avail
