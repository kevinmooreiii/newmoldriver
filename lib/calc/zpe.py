"""
  ZPE Calculations
"""


# New
from lib.submission import substr
from routines.pf import get_zero_point_energy


def get_zpe(spc, spc_info, spc_save_path, pf_levels, spc_model):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    spc_zpe = {}
    is_atom = {}
    zero_energy_str = {}

    spc_zpe, is_atom = get_zero_point_energy(
        spc, spc_info, pf_levels, spc_model,
        pf_script_str=substr.MESSPF,
        elec_levels=[[0., 1]],
        sym_factor=1.0,
        save_prefix=spc_save_path)
    zpe_str = '{0:<8.2f}\n'.format(spc_zpe)
    if is_atom:
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
        zero_energy_str += 'End'

    return spc_zpe, zero_energy_str
