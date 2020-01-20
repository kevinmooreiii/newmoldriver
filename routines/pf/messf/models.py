"""
  Build the species string based on the model
"""

import elstruct
import automol
from routines.pf.messf import _tors as tors
from routines.pf.messf import _vib as vib


def vib_harm_tors_rigid(spc_info, harm_min_cnf_locs, harm_cnf_save_fs,
                        saddle=False):
    """ Build the species string for a model of
        harmonic vibrational frequencies and rigid torsions
    """
    # Do the freqs obtain for two species for fake and pst
    if harm_min_cnf_locs is not None:
        # Obtain geom and freqs from filesys
        harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(harm_min_cnf_locs)
        hess = harm_cnf_save_fs.leaf.file.hessian.read(harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)
        # Modify freqs lst and get imaginary frequencies
        mode_start = 6
        if saddle:
            mode_start = mode_start + 1
            imag_freq = freqs[0]
        if automol.geom.is_linear(harm_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies ',
              'for species {}'.format(spc_info[0]))
        raise ValueError

    return harm_geo, freqs, imag_freq


def vib_harm_tors_1dhr(harm_min_cnf_locs, harm_cnf_save_fs,
                       tors_min_cnf_locs, tors_cnf_save_fs,
                       tors_save_path, tors_cnf_save_path,
                       spc_dct_i, spc_info,
                       frm_bnd_key, brk_bnd_key,
                       sym_factor, elec_levels,
                       projrot_script_str,
                       hind_rot_geo=False,
                       saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs.leaf.file.energy.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs.leaf.file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        if tors_min_cnf_locs is not None:

            # Get geometry for the torsional minimum
            zma = tors_cnf_save_fs.leaf.file.zmatrix.read(
                tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(
                tors_min_cnf_locs)

            # Set torsional stuff
            tors_names = tors.get_tors_names(
                spc_dct_i, tors_cnf_save_fs, saddle=saddle)
            tors_grids = tors.get_tors_grids(
                spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key)

            # Set ts bond
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

            hind_rot_str, proj_rotors_str = tors.write_1dhr_tors_mess_strings(
                harm_geo, spc_info, sym_num, spc_dct_i,
                tors_names, tors_grids, ts_bnd, zma,
                tors_cnf_save_path, min_ene,
                saddle=False, hind_rot_geo=None)

            # Calculate ZPVES of the hindered rotors
            if saddle and tors_names is not None:
                tors_zpe = tors.calc_tors_freqs_zpe(
                    tors_geo, sym_factor, elec_levels,
                    hind_rot_str, tors_save_path)
            else:
                tors_zpe = 0.0

            # Run one vers ProjRot to proj freqs for that version
            freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
                tors_geo, hess, pot,
                proj_rotors_str, projrot_script_str,
                tors_save_path, saddle=False)

            # Now run the other version of ProjRot
            pfreqs2 = vib.projrot_freqs_2(
                tors_save_path, pot, saddle=saddle)
            [freqs2, imag_freq2,
             zpe_harm_no_tors_2, harm_zpe] = pfreqs2

            # Determine freqs and imag_freqs
            freqs, imag_freq, zpe = vib.determine_freqs_zpe(
                freqs1, freqs2, imag_freq1, imag_freq2,
                zpe_harm_no_tors, zpe_harm_no_tors_2,
                harm_zpe, tors_zpe)
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        tors_geo, freqs, imag_freq, hind_rot_str = (), (), 0.0, ''
        raise ValueError

    return tors_geo, freqs, imag_freq, hind_rot_str, zpe, sym_factor


def vib_harm_tors_mdhr_tau(harm_min_cnf_locs, harm_cnf_save_fs,
                           tors_min_cnf_locs, tors_cnf_save_fs,
                           tors_save_path, tors_cnf_save_path,
                           spc_dct_i, spc_info,
                           frm_bnd_key, brk_bnd_key,
                           sym_factor, elec_levels,
                           projrot_script_str,
                           hind_rot_geo=False,
                           saddle=False):
    """ Build the species string for a model: Harm, 1DHR
    """
    if harm_min_cnf_locs is not None:
        harm_geo = harm_cnf_save_fs.leaf.file.geometry.read(
            harm_min_cnf_locs)
        min_ene = harm_cnf_save_fs.leaf.file.energy.read(
            harm_min_cnf_locs)
        hess = harm_cnf_save_fs.leaf.file.hessian.read(
            harm_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(
            harm_geo, hess, project=False)

        if tors_min_cnf_locs is not None:

            # Get geometry for the torsional minimum
            zma = tors_cnf_save_fs.leaf.file.zmatrix.read(
                tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(
                tors_min_cnf_locs)

            # Set torsional stuff
            tors_names = tors.get_tors_names(
                spc_dct_i, tors_cnf_save_fs, saddle=saddle)
            tors_grids = tors.get_tors_grids(
                spc_dct_i, zma, tors_names, frm_bnd_key, brk_bnd_key)

            # Set ts bond
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)

            # Get the tau info string
            tau_inf_str = writer_tau_pf()

            # Calculate ZPVES of the hindered rotors
            if saddle and tors_names is not None:
                tors_zpe = tors.calc_tors_freqs_zpe(
                    tors_geo, sym_factor, elec_levels,
                    tau_inf_str, tors_save_path)
            else:
                tors_zpe = 0.0

            # Run one vers ProjRot to proj freqs for that version
            freqs1, imag_freq1, zpe_harm_no_tors = vib.projrot_freqs_1(
                tors_geo, hess, pot,
                proj_rotors_str, projrot_script_str,
                tors_save_path, saddle=False)

            # Now run the other version of ProjRot
            pfreqs2 = vib.projrot_freqs_2(
                tors_save_path, pot, saddle=saddle)
            [freqs2, imag_freq2,
             zpe_harm_no_tors_2, harm_zpe] = pfreqs2

            # Determine freqs and imag_freqs
            freqs, imag_freq, zpe = vib.determine_freqs_zpe(
                freqs1, freqs2, imag_freq1, imag_freq2,
                zpe_harm_no_tors, zpe_harm_no_tors_2,
                harm_zpe, tors_zpe)
    else:
        print('ERROR: Reference geometry is missing for harmonic frequencies',
              ' for species {}'.format(spc_info[0]))
        tors_geo, freqs, imag_freq, hind_rot_str = (), (), 0.0, ''
        raise ValueError

    return tors_geo, freqs, imag_freq, hind_rot_str, zpe, sym_factor


def get_stoich(harm_min_cnf_locs_i, harm_min_cnf_locs_j,
               harm_cnf_save_fs_i, harm_cnf_save_fs_j):
    """ get the overall combined stoichiometry
    """
    if harm_min_cnf_locs_i is not None:
        harm_geo_i = harm_cnf_save_fs_i.leaf.file.geometry.read(
            harm_min_cnf_locs_i)
        if harm_min_cnf_locs_j is not None:
            harm_geo_j = harm_cnf_save_fs_j.leaf.file.geometry.read(
                harm_min_cnf_locs_j)

    form_i = automol.geom.formula(harm_geo_i)
    form_j = automol.geom.formula(harm_geo_j)
    form = automol.formula.join(form_i, form_j)
    stoich = ''
    for key, val in form.items():
        stoich += key + str(val)

    return stoich
