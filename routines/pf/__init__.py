""" Libraries of physical data used by the drivers
"""

from routines.pf import rates
from routines.pf import fit
from routines.pf import messf
from routines.pf.messf import get_zero_point_energy
from routines.pf.messf import get_high_level_energy


__all__ = [
    'rates',
    'fit',
    'messf',
    'get_zero_point_energy',
    'get_high_level_energy'
]


def write_channel_mess_strs(spc_dct, rxn_lst, pes_formula,
                            ts_model, pf_levels, multi_info, pst_params,
                            spc_save_fs, save_prefix, idx_dct):
    """ Write all the MESS input file strings for the reaction channels
    """
    mess_strs = ['', '', '']

    for spc in full_queue:
        spc_info = (spc_dct[spc]['ich'],
                    spc_dct[spc]['chg'],
                    spc_dct[spc]['mul'])
        print('SPC CHECK')
        print(spc)
        print(spc_dct)
        if 'ts_' in spc:
            spc_dct[spc] = lmech.set_sadpt_info(
                ts_tsk_lst, spc_dct, spc, thy_dct,
                run_prefix, save_prefix,
                kickoff=(0.1, False))
            spc_save_path = spc_dct[spc]['rxn_fs'][3]
            saddle = True
            save_path = spc_save_path
        else:
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            saddle = False
            save_path = save_prefix

        # Cacluate the ZPVE and set in the dict
        spc_dct[spc]['zpe'], _ = calczpe.get_zpe(
            spc, spc_dct[spc], spc_save_path, pf_levels, pf_model)

        # Set ene and string
        spc_dct[spc]['ene'] = routines.pf.get_high_level_energy(
            spc_info=spc_info,
            thy_low_level=model_dct['geo'],
            thy_low_level=finf.get_thy_info(model_dct['geo'], thy_dct)
            thy_high_level=finf.get_thy_info(model_dct['ene'], thy_dct)
            save_prefix=save_path,
            saddle=saddle)

    first_ground_ene = 0.
    species = messrates.make_all_species_data(
        rxn_lst, spc_dct, save_prefix, ts_model, pf_levels)
    for idx, rxn in enumerate(rxn_lst):
        tsname = 'ts_{:g}'.format(idx)
        tsform = automol.geom.formula(
            automol.zmatrix.geometry(spc_dct[tsname]['original_zma']))
        if tsform != pes_formula:
            print('Reaction ist contains reactions on different potential',
                  'energy surfaces: {} and {}'.format(tsform, pes_formula))
            print('Will proceed to construct only {}'.format(pes_formula))
            continue
        mess_strs, first_ground_ene = messrates.make_channel_pfs(
            tsname, rxn, species, spc_dct, idx_dct, mess_strs,
            first_ground_ene, spc_save_fs, ts_model, pf_levels,
            multi_info,
            pst_params=pst_params)
    well_str, bim_str, ts_str = mess_strs
    ts_str += '\nEnd\n'
    print(well_str)
    print(bim_str)
    print(ts_str)

    return well_str, bim_str, ts_str
