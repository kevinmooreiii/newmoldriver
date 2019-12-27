"""
  Functions to read the filesystem and pull objects from it
"""

import automol
import autofile
from routines import util


def get_zmas(
        reacs, prods, spc_dct, ini_thy_info, save_prefix, run_prefix,
        kickoff_size, kickoff_backward, projrot_script_str):
    """get the zmats for reactants and products using the initial level of theory
    """
    if len(reacs) > 2:
        ich = spc_dct[reacs[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        reacs = reacs[:-1]
    elif len(prods) > 2:
        ich = spc_dct[prods[-1]]['ich']
        ichgeo = automol.inchi.geometry(ich)
        ichzma = automol.geom.zmatrix(ichgeo)
        prods = prods[:-1]
    rct_geos, rct_cnf_save_fs_lst = get_geos(
        reacs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    prd_geos, prd_cnf_save_fs_lst = get_geos(
        prods, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str)
    rct_zmas = list(map(automol.geom.zmatrix, rct_geos))
    prd_zmas = list(map(automol.geom.zmatrix, prd_geos))
    # for geo in prd_geos:
    #     xyzs = automol.geom.coordinates(geo)
    if len(rct_zmas) > 2:
        rct_zmas.append(ichzma)
    if len(prd_zmas) > 2:
        prd_zmas.append(ichzma)
    return rct_zmas, prd_zmas, rct_cnf_save_fs_lst, prd_cnf_save_fs_lst


def get_geos(
        spcs, spc_dct, ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward, projrot_script_str):
    """get geos for reactants and products using the initial level of theory
    """
    spc_geos = []
    cnf_save_fs_lst = []
    for spc in spcs:
        spc_info = [spc_dct[spc]['ich'],
                    spc_dct[spc]['chg'],
                    spc_dct[spc]['mul']]
        orb_restr = util.orbital_restriction(spc_info, ini_thy_info)
        ini_thy_level = ini_thy_info[0:3]
        ini_thy_level.append(orb_restr)
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs.leaf.create(spc_info)
        spc_run_path = spc_run_fs.leaf.path(spc_info)
        ini_thy_save_fs = autofile.fs.theory(spc_save_path)
        ini_thy_save_path = ini_thy_save_fs.leaf.path(ini_thy_level[1:4])
        ini_thy_run_fs = autofile.fs.theory(spc_run_path)
        ini_thy_run_path = ini_thy_run_fs.leaf.path(ini_thy_level[1:4])
        cnf_save_fs = autofile.fs.conformer(ini_thy_save_path)
        cnf_save_fs_lst.append(cnf_save_fs)
        cnf_run_fs = autofile.fs.conformer(ini_thy_run_path)
        min_cnf_locs = util.min_energy_conformer_locators(cnf_save_fs)
        if min_cnf_locs:
            geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        else:
            run_fs = autofile.fs.run(ini_thy_run_path)
            run_fs.trunk.create()
            tmp_ini_fs = [None, ini_thy_save_fs]
            tmp_fs = [spc_save_fs, spc_run_fs, ini_thy_save_fs, ini_thy_run_fs,
                      cnf_save_fs, cnf_run_fs, run_fs]
            geo = geom.reference_geometry(
                spc_dct[spc], ini_thy_level, ini_thy_level, tmp_fs, tmp_ini_fs,
                kickoff_size, kickoff_backward, projrot_script_str,
                overwrite=False)
        spc_geos.append(geo)
    return spc_geos, cnf_save_fs_lst
