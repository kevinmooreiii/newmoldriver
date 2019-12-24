"""
Build paths and file systesm given species and theory
'info' objects
"""

import automol
import autofile

from lib.phydat import phycon
from lib.submission import theolvls
from lib import moldr


def get_thy_info(method):
    """ convert theory level dictionary to theory information array
    """
    lvldic = theolvls.ES_DCT[method]
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


def get_es_info(method):
    """
    Turn es dictionary in theory info array
    """
    if method == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = get_thy_info(method)
    return ret


def get_spc_info(spc_dct):
    """ convert species dictionary to species_info array
    """
    err_msg = ''
    props = ['ich', 'chg', 'mul']
    for i, prop in enumerate(props):
        if prop in spc_dct:
            props[i] = spc_dct[prop]
        else:
            err_msg = prop
    if err_msg:
        print('ERROR: No {} found'.format(err_msg))
    return props


def rxn_info(run_prefix, save_prefix, ts, spc_dct,
             thy_info, ini_thy_info=None):
    """ prepare rxn info and reverse the reactants and products
        if reaction is endothermic
    """
    rxn_ichs = [[], []]
    rxn_chgs = [[], []]
    rxn_muls = [[], []]
    print('\n TS for {}: {} = {}'.format(
        ts, '+'.join(spc_dct[ts]['reacs']), '+'.join(spc_dct[ts]['prods'])))
    reacs = spc_dct[ts]['reacs']
    prods = spc_dct[ts]['prods']
    print('ts dct', spc_dct[ts])
    for spc in reacs:
        rxn_ichs[0].append(spc_dct[spc]['ich'])
        rxn_chgs[0].append(spc_dct[spc]['chg'])
        rxn_muls[0].append(spc_dct[spc]['mul'])
    for spc in prods:
        rxn_ichs[1].append(spc_dct[spc]['ich'])
        rxn_chgs[1].append(spc_dct[spc]['chg'])
        rxn_muls[1].append(spc_dct[spc]['mul'])
    # check direction of reaction
    try:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, thy_info)
    except:
        rxn_exo = moldr.util.reaction_energy(
            save_prefix, rxn_ichs, rxn_chgs, rxn_muls, ini_thy_info)
    print('reaction is {:.2f} endothermic'.format(rxn_exo*phycon.EH2KCAL))
    if rxn_exo > 0 and not spc_dct[ts]['given_class']:
        rxn_ichs = rxn_ichs[::-1]
        rxn_chgs = rxn_chgs[::-1]
        rxn_muls = rxn_muls[::-1]
        spc_dct[ts]['reacs'] = prods
        spc_dct[ts]['prods'] = reacs
        print('Reaction will proceed as {}: {} = {}'.format(
            ts, '+'.join(spc_dct[ts]['reacs']),
            '+'.join(spc_dct[ts]['prods'])))
        # print('ts search will be performed in reverse direction')

    # set up the filesystem
    rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
        rxn_ichs, rxn_chgs, rxn_muls)
    low_mul = min(
        automol.mult.ts._low(rxn_muls[0]),
        automol.mult.ts._low(rxn_muls[1]))
    high_mul = max(
        automol.mult.ts._high(rxn_muls[0]),
        automol.mult.ts._high(rxn_muls[1]))

    return rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul
