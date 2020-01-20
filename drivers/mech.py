""" estoktp driver functions with no
better home
"""

from drivers import esdriver
from drivers import ktpdriver
from lib.filesystem import inf as finf
from lib.filesystem import path as fpath
from lib.filesystem import read as fread
from lib.reaction import rxnid
from lib.load import species as loadspc


def run_driver(pes_dct, conn_chnls_dct,
               spc_dct, cla_dct,
               thy_dct,
               model_dct,
               es_tsk_lst,
               run_jobs_lst,
               run_inp_dct,
               run_options_dct,
               driver='es_spc'):
    """ Run the ktp driver for the PESs
    """
    for formula in pes_dct:

        # Build the names list
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Select names from the names list corresponding to chnls to run
        # conn_chnls_dct[formula] = {sub_pes_idx: [channel_idxs]}
        for cvals in conn_chnls_dct[formula].values():
            run_pes = False
            rct_names_lst = []
            prd_names_lst = []
            rxn_name_lst = []
            for chn_idx, _ in enumerate(pes_rxn_name_lst):
                if chn_idx in cvals:
                    run_pes = True
                    rct_names_lst.append(pes_rct_names_lst[chn_idx])
                    prd_names_lst.append(pes_prd_names_lst[chn_idx])
                    rxn_name_lst.append(pes_rxn_name_lst[chn_idx])
                    print('running channel {}: {} = {}'.format(
                        str(chn_idx+1),
                        ' + '.join(pes_rct_names_lst[chn_idx]),
                        ' + '.join(pes_prd_names_lst[chn_idx])))

            # Call the Driver for the PES
            if run_pes:

                # Form the reaction list
                rxn_lst = format_run_rxn_lst(rct_names_lst, prd_names_lst)

                # Format the task info list
                spc_tsk_lst, ts_tsk_lst = format_tsk_lst(es_tsk_lst)

                # Add the stationary points to the spc dcts
                print('\nBegin transition state prep')
                ts_dct = loadspc.build_spc_dct_for_sadpts(
                    spc_dct, rxn_lst, rxn_name_lst,
                    rct_names_lst, prd_names_lst, cla_dct)
                for sadpt in ts_dct:
                    ts_dct[sadpt] = set_sadpt_info(
                        ts_tsk_lst, ts_dct, spc_dct, sadpt,
                        run_inp_dct['run_prefix'],
                        run_inp_dct['save_prefix'],
                        run_options_dct['kickoff'])
                    print('loop dct\n', ts_dct[sadpt]['dist_info'])
                spc_dct.update(ts_dct)

                # Run the appropriate driver
                if driver == 'es':
                    esdriver.run(
                        rxn_lst, spc_dct,
                        spc_tsk_lst,
                        model_dct,
                        run_options_dct, run_inp_dct)
                elif driver == 'ktp':
                    ktpdriver.run(
                        formula,
                        spc_dct,
                        thy_dct,
                        rxn_lst,
                        model_dct,
                        run_inp_dct,
                        run_rates=bool('rates' in run_jobs_lst),
                        run_fits=bool('params' in run_jobs_lst))


def format_tsk_lst(tsk_info_lst):
    """ Format the input task list appropriate for spc or rxns
    """
    spc_tsk_lst = []
    ts_tsk_lst = []
    ts_tsk = False
    for tsk in tsk_info_lst:
        if 'find_ts' in tsk[0]:
            ts_tsk = True
        if ts_tsk:
            ts_tsk_lst.append(tsk)
        else:
            spc_tsk_lst.append(tsk)

    return spc_tsk_lst, ts_tsk_lst


def format_run_rxn_lst(rct_names_lst, prd_names_lst):
    """ Get the lst of reactions to be run
    """
    # spc_queue
    spc_queue = []
    for rxn, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)
    run_lst = []
    for rxn, _ in enumerate(rct_names_lst):
        spc_queue = []
        rxn_spc = list(rct_names_lst[rxn])
        rxn_spc.extend(list(prd_names_lst[rxn]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)
        run_lst.append(
            {'species': spc_queue, 'reacs': list(rct_names_lst[rxn]), 'prods':
             list(prd_names_lst[rxn])})

    return run_lst


def set_sadpt_info(ts_tsk_lst, ts_dct, spc_dct, sadpt,
                   run_prefix, save_prefix, kickoff):
    """ set the saddle point dct with info
    """
    ini_thy_info = ts_tsk_lst[1][2]
    thy_info = ts_tsk_lst[1][1]
    print('ini_thy_info', ini_thy_info)
    print('thy_info', thy_info)
    # Generate rxn data, reorder if necessary, and put in spc_dct for given ts
    rxn_ichs, rxn_chgs, rxn_muls, low_mul, high_mul = finf.rxn_info(
        save_prefix, sadpt, ts_dct, spc_dct, thy_info, ini_thy_info)
    ts_dct[sadpt]['rxn_ichs'] = rxn_ichs
    ts_dct[sadpt]['rxn_chgs'] = rxn_chgs
    ts_dct[sadpt]['rxn_muls'] = rxn_muls
    ts_dct[sadpt]['low_mul'] = low_mul
    ts_dct[sadpt]['high_mul'] = high_mul

    # Generate rxn_fs from rxn_info stored in spc_dct
    [kickoff_size, kickoff_backward] = kickoff
    rxn_run_fs, rxn_save_fs, rxn_run_path, rxn_save_path = fpath.get_rxn_fs(
        run_prefix, save_prefix, ts_dct[sadpt])
    ts_dct[sadpt]['rxn_fs'] = [
        rxn_run_fs,
        rxn_save_fs,
        rxn_run_path,
        rxn_save_path]
    rct_zmas, prd_zmas, rct_cnf_save_fs, prd_cnf_save_fs = fread.get_zmas(
        ts_dct[sadpt]['reacs'], ts_dct[sadpt]['prods'], spc_dct,
        ini_thy_info, save_prefix, run_prefix, kickoff_size,
        kickoff_backward)
    ret = rxnid.ts_class(
        rct_zmas, prd_zmas, ts_dct[sadpt]['rad_rad'],
        ts_dct[sadpt]['mul'], low_mul, high_mul,
        rct_cnf_save_fs, prd_cnf_save_fs, ts_dct[sadpt]['given_class'])
    ret1, ret2 = ret
    if ret1:
        [rxn_class, spc_zma,
         dist_name, brk_name, grid,
         frm_bnd_key, brk_bnd_key,
         tors_names, update_guess] = ret1
        ts_dct[sadpt]['class'] = rxn_class
        ts_dct[sadpt]['grid'] = grid
        ts_dct[sadpt]['tors_names'] = tors_names
        ts_dct[sadpt]['original_zma'] = spc_zma
        dist_info = [dist_name, 0., update_guess, brk_name]
        ts_dct[sadpt]['dist_info'] = dist_info
        ts_dct[sadpt]['frm_bnd_key'] = frm_bnd_key
        ts_dct[sadpt]['brk_bnd_key'] = brk_bnd_key
        # Adding in the rct and prd zmas for vrctst
        ts_dct[sadpt]['rct_zmas'] = rct_zmas
        ts_dct[sadpt]['prd_zmas'] = prd_zmas
        if ret2:
            ts_dct[sadpt]['bkp_data'] = ret2
        else:
            ts_dct[sadpt]['bkp_data'] = None
    else:
        ts_dct[sadpt]['class'] = None
        ts_dct[sadpt]['bkp_data'] = None

    return ts_dct[sadpt]


def set_pf_model_info(pf_model):
    """ Set the PF model list based on the input
    """
    tors_model = pf_model['tors'] if pf_model['tors'] else 'RIGID'
    vib_model = pf_model['vib'] if pf_model['vib'] else 'HARM'
    sym_model = pf_model['sym'] if pf_model['sym'] else ''

    pf_models = [tors_model, vib_model, sym_model]
    return pf_models


def set_es_model_info(es_model, thy_dct):
    """ Set the model info
    """
    # Read the ES models from the model dictionary
    geo_lvl = es_model['geo'] if es_model['geo'] else None
    ene_lvl = es_model['ene'] if es_model['ene'] else None
    harm_lvl = es_model['harm'] if es_model['harm'] else None
    anharm_lvl = es_model['anharm'] if es_model['anharm'] else None
    sym_lvl = es_model['sym'] if es_model['sym'] else None

    # Torsional Scan which needs a reference for itself
    tors_lvl_sp = es_model['tors'][0] if es_model['tors'] else None
    tors_lvl_scn = es_model['tors'][1] if es_model['tors'] else None

    # Set the theory info objects
    geo_thy_info = finf.get_thy_info(geo_lvl, thy_dct)
    ene_thy_info = finf.get_thy_info(ene_lvl, thy_dct)
    harm_thy_info = finf.get_thy_info(harm_lvl, thy_dct)
    anharm_thy_info = (finf.get_thy_info(anharm_lvl, thy_dct)
                       if anharm_lvl else None)
    sym_thy_info = (finf.get_thy_info(sym_lvl, thy_dct)
                    if sym_lvl else None)
    tors_sp_thy_info = (finf.get_thy_info(tors_lvl_sp, thy_dct)
                        if tors_lvl_sp else None)
    tors_scn_thy_info = (finf.get_thy_info(tors_lvl_scn, thy_dct)
                         if tors_lvl_scn else None)

    # Combine levels into a list
    es_levels = [
        geo_thy_info,
        ene_thy_info,
        harm_thy_info,
        anharm_thy_info,
        sym_thy_info,
        [tors_sp_thy_info, tors_scn_thy_info]
    ]

    return es_levels
