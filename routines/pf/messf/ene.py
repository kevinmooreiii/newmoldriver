""" calculates certain quantities of interest using MESS+filesytem
"""

import automol
import elstruct
import autofile

# New Libs
from lib.phydat import phycon
from lib.filesystem import minc as fsmin
from lib.filesystem import orb as fsorb
from lib.filesystem import inf as finf
import routines.pf.messf.models as pfmodels
from routines.pf.messf.blocks import set_model_filesys


def get_zpe_str(spc_dct, zpe):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    if automol.geom.is_atom(automol.inchi.geom(spc_dct['ich'])):
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + str(zpe)
        zero_energy_str += 'End'

    return zero_energy_str


def get_zero_point_energy(spc, spc_dct_i, pf_levels, spc_model, save_prefix):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    # Set up input
    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])

    # Prepare the sets of file systems
    [_, _, harm_levels, _, _, tors_levels] = pf_levels
    print('spc_model')
    print(spc_model)
    tors_model, vib_model, _ = spc_model

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_prefix)

    # Set boolean to account for rad-rad reaction (not supported by vtst)
    rad_rad_ts = False
    saddle = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True
        saddle = True

    # Set the filesystem objects for various species models
    harmfs = set_model_filesys(
        thy_save_fs, spc_info, harm_levels, saddle=saddle)
    [harm_cnf_save_fs, _,
     harm_min_cnf_locs, _] = harmfs
    if tors_levels and not rad_rad_ts:
        torsfs = set_model_filesys(
            thy_save_fs, spc_info, tors_levels[0], saddle=saddle)
        [tors_cnf_save_fs, tors_cnf_save_path,
         tors_min_cnf_locs, tors_save_path] = torsfs

    # if saddle:
    if 'ts_' in spc:
        frm_bnd_key = spc_dct_i['frm_bnd_key']
        brk_bnd_key = spc_dct_i['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []

    # Get reference harmonic
    harm_zpe = 0.0
    is_atom = False
    if not harm_min_cnf_locs:
        print('ERROR: No harmonic reference geometry for this species ',
              '{}'.format(spc_info[0]))
        return harm_zpe, is_atom
    harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(harm_min_cnf_locs)
    if automol.geom.is_atom(harm_geo):
        zpe = 0.0
        is_atom = True
    else:
        hess = harm_cnf_save_fs.leaf.file.hessian.read(harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        saddle = False
        mode_start = 6
        if 'ts_' in spc:
            mode_start = mode_start + 1
            saddle = True
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]

        harm_zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

        # Determine the ZPVE based on the model
        print(vib_model)
        print(tors_model)
        print(rad_rad_ts)
        if (vib_model == 'harm' and tors_model == 'rigid') or rad_rad_ts:
            print('HARM_RIGID')
            zpe = harm_zpe
        elif vib_model == 'harm' and tors_model == '1dhr':
            print('HARM_1DHR')
            _, _, _, _, zpe = pfmodels.vib_harm_tors_1dhr(
                harm_min_cnf_locs, harm_cnf_save_fs,
                tors_min_cnf_locs, tors_cnf_save_fs,
                tors_save_path, tors_cnf_save_path,
                spc_dct_i, spc_info,
                frm_bnd_key, brk_bnd_key,
                sym_factor, elec_levels,
                saddle=saddle)
        elif vib_model == 'harm' and tors_model == 'mdhr':
            print('HARM and MDHR combination is not yet implemented')
        elif vib_model == 'harm' and tors_model == 'tau':
            print('HARM and TAU combination is not yet implemented')
        elif vib_model == 'tau' and tors_model == 'tau':
            zpe = 0.0
        elif vib_model == 'vpt2' and tors_model == 'rigid':
            print('VPT2 and RIGID combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == '1dhr':
            print('VPT2 and 1DHR combination is not yet implemented')
        elif vib_model == 'vpt2' and tors_model == 'tau':
            print('VPT2 and TAU combination is not yet implemented')

    return zpe


def get_high_level_energy(
        spc_info, thy_low_level, thy_high_level, save_prefix, saddle=False):
    """ get high level energy at low level optimized geometry
    """
    if saddle:
        spc_save_path = save_prefix
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_low_level)
    thy_low_level = thy_low_level[1:3]
    thy_low_level.append(orb_restr)

    ll_save_fs = autofile.fs.theory(spc_save_path)
    ll_save_path = ll_save_fs.leaf.path(thy_low_level)

    if saddle:
        ll_save_fs = autofile.fs.ts(ll_save_path)
        ll_save_fs.trunk.create()
        ll_save_path = ll_save_fs.trunk.path()

    cnf_save_fs = autofile.fs.conformer(ll_save_path)
    min_cnf_locs = fsmin.min_energy_conformer_locators(
        cnf_save_fs)
    if not min_cnf_locs:
        print('ERROR: No minimum conformer geometry for ',
              'this species {}'.format(spc_info[0]))
        return 0.0
    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)

    orb_restr = fsorb.orbital_restriction(
        spc_info, thy_high_level)
    thy_high_level = thy_high_level[1:3]
    thy_high_level.append(orb_restr)

    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs.leaf.create(thy_high_level)

    print('test path')
    print(thy_high_level)
    print(sp_save_fs.leaf.path(thy_high_level))
    min_ene = sp_save_fs.leaf.file.energy.read(thy_high_level)

    return min_ene


def calc_shift_ene(spc_dct, thy_dct, model_dct, rxn,
                   spc_tgt, spc_info_tgt, spc_save_fs_tgt,
                   pf_levels1, spc_model1, pf_levels2, spc_model2,
                   save_prefix_tgt, saddle_tgt):
    """ Function to shift the energy of a species to allow mixing of
        channels calculated with two different methods.
        We consider two levels of theory:
          (1) method for PES reference species and
          (2) method the user requested for this point on the PES.
        Shifted energy will be calculated as:
          E_t[1] = {E_s[1] + (E_t[2] - E_s[2])} - E_r[1]
        where E_t = ene of target species missing info at lvl 1,
              E_s = ene of other species on channel, and
              E_r = ene of species that is reference for whole PES.
        Note that this function only returns the part in curly braces and the
        final E_r[1] term is subtracted in some other part of the code.
    """

    # Read the energy for the target species from the filesystem
    tgt_ene_lvl2 = get_fs_ene_zpe(spc_dct, spc_tgt,
                                  spc_info_tgt, spc_save_fs_tgt,
                                  thy_dct, model_dct, pf_levels2, spc_model2,
                                  save_prefix_tgt, saddle=saddle_tgt)

    # Loop over the species in the channel and find one species
    # where the energy and ZPVE has been calculated at levels 1 and 2.
    chn_enes = {}
    for spc in rxn:
        # Maybe get the spc and other stuff...
        spc_info = 1
        spc_save_fs = 1
        save_prefix = 1
        saddle = bool('ts_' in spc)
        # Try and read the energies from the filesystem
        chn_ene1 = get_fs_ene_zpe(spc_dct, spc, spc_info, spc_save_fs,
                                  thy_dct, model_dct, pf_levels1, spc_model1,
                                  save_prefix, saddle=saddle)
        chn_ene2 = get_fs_ene_zpe(spc_dct, spc, spc_info, spc_save_fs,
                                  thy_dct, model_dct, pf_levels2, spc_model2,
                                  save_prefix, saddle=saddle)
        # Only add the energies to both dcts if ene1 and ene2 were found
        if chn_ene1 and chn_ene2:
            chn_enes[spc] = (chn_ene1, chn_ene2)

    # Calculate the shifted energy for the species at level 1, if possible
    if chn_enes:
        # Get lvl1 and lvl2 enes for 1st spc in dcts (shouldn't matter which)
        for spc, enes in chn_enes.values():
            chn_ene_lvl1, chn_ene_lvl2 = enes
            break
        # Calculate energy
        tgt_ene_lvl1 = chn_ene_lvl1 + (chn_ene_lvl2 - tgt_ene_lvl2)
    else:
        print('No species on the channel with energies at methods 1 and 2')
        tgt_ene_lvl1 = None

    return tgt_ene_lvl1


def get_fs_ene_zpe(spc_dct, spc, spc_info, spc_save_fs,
                   thy_dct, model_dct, pf_levels, spc_model,
                   save_prefix, saddle=False):
    """ Get the energy for a species on a channel
    """
    # Set paths
    if saddle:
        spc_save_path = spc_dct[spc]['rxn_fs'][3]
        save_path = spc_save_path
    else:
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)
        save_path = save_prefix

    # Read the electronic energy and ZPVE
    e_elec = get_high_level_energy(
        spc_info=spc_info,
        thy_low_level=finf.get_thy_info(model_dct['geo'], thy_dct),
        thy_high_level=finf.get_thy_info(model_dct['ene'], thy_dct),
        save_prefix=save_path,
        saddle=saddle)
    e_zpe = get_zero_point_energy(
        spc, spc_dct, pf_levels, spc_model,
        save_prefix=spc_save_path)

    return (e_elec + e_zpe) * phycon.EH2KCAL
