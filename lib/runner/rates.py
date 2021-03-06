"""
  Run rates with MESS
"""

import os
import autofile

# New Libs
from lib.runner import script
from lib.filesystem import orb as fsorb


def run_rates(
        header_str, energy_trans_str, well_str, bim_str, ts_str,
        tsdct, rxn_save_path,
        model_dct, thy_dct, run_model):
    """ Generate k(T,P) by first compiling all the MESS strings
        and then running MESS
    """
    pf_levels = lmech.set_es_model_info(model_dct['es'], thy_dct)[0]

    ts_info = (tsdct['ich'], tsdct['chg'], tsdct['mul'])
    orb_restr = fsorb.orbital_restriction(ts_info, thy_info)
    ref_level = thy_info[1:3]
    ref_level.append(orb_restr)
    print('ref level test:', ref_level)
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(ref_level)
    thy_save_path = thy_save_fs.leaf.path(ref_level)

    mess_inp_str = '\n'.join(
        [header_str, energy_trans_str, well_str, bim_str, ts_str])

    print('mess input file')
    print(mess_inp_str)

    bld_locs = ['MESS', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs.leaf.create(bld_locs)
    mess_path = bld_save_fs.leaf.path(bld_locs)
    print('Build Path for MESS rate files:')
    print(mess_path)
    with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
        mess_file.write(mess_inp_str)
    script.run_script(script.MESSRATE, mess_path)
    return mess_path
