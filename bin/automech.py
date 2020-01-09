"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import os
import sys
from drivers import mech as lmech
# from drivers import thermodriver
from lib.load import run as loadrun
from lib.load import theory as loadthy
from lib.load import model as loadmodel
from lib.load import mechanism as loadmech
from lib.load import species as loadspc
from lib.filesystem import build as lfs
from lib.submission import read_dat


# Set runtime options based on user input
JOB_PATH = sys.argv[1]
PARAMS = read_dat.params(os.path.join(JOB_PATH, 'inp/params.dat'))

# Parse the run input
RUN_INP_DCT = loadrun.build_run_inp_dct(JOB_PATH)
RUN_LST = loadrun.objects_lst(JOB_PATH)
print(RUN_LST)
[PESNUMS, CHANNELS, MODEL] = RUN_LST[0]

RUN_GLOB_OPTS_LST = loadrun.build_run_glob_opts_dct(JOB_PATH)
RUN_JOBS_LST = loadrun.build_run_jobs_lst(JOB_PATH)
RUN_ES_TSKS = loadrun.read_es_tsks(JOB_PATH)

# Parse the theory input
THY_DCT = loadthy.build_thy_dct(JOB_PATH)

# Parse the model input
ALL_MODEL_DCT = loadmodel.read_models_sections(JOB_PATH)
MODEL_DCT = ALL_MODEL_DCT[MODEL]
for x,y in MODEL_DCT.items():
    print(x, y)
print('\n\n')

# Parse the species input
SPC_DCT = loadspc.build_spc_dct(JOB_PATH, RUN_INP_DCT['spc'])
MOD_SPC_DCT = loadspc.modify_spc_dct(JOB_PATH, SPC_DCT, PARAMS.HIND_INC)

# Parse the mechanism input
PES_DCT = loadmech.parse_mechanism_file(
    JOB_PATH, RUN_INP_DCT['mech'], MOD_SPC_DCT,
    check_stereo=RUN_INP_DCT['check_stereo'],
    sort_rxns=RUN_INP_DCT['sort_rxns'],
    rad_rad_sort=RUN_INP_DCT['rad_rad_sort'])

# Set the channels on the PESs to run
loadmech.print_pes_channels(PES_DCT)
# PESNUMS_LST = loadmech.get_pes_nums(PES_DCT, PESNUMS)
CONN_CHNLS_LST = loadmech.determine_connected_pes_channels(
    PES_DCT, [PESNUMS])

# Format the task list
ES_TSK_LST = loadrun.build_run_es_tsks_lst(RUN_ES_TSKS, MODEL_DCT)
import sys
sys.exit()

# Prepare run-save filesystem for job
lfs.prefix_filesystem(RUN_INP_DCT['run_prefix'], RUN_INP_DCT['save_prefix'])

# Run the requested drivers: es, thermo, ktp
if PARAMS.RUN_ES_SPC:
    lmech.run_driver(
        PES_DCT, [PESNUMS], CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        ES_TSK_LST,
        RUN_INP_DCT['run_prefix'], RUN_INP_DCT['save_prefix'],
        options=PARAMS.OPTIONS_RATE,
        hind_inc=PARAMS.HIND_INC,
        mc_nsamp=PARAMS.MC_NSAMP0,
        multi_info=PARAMS.MULTI_INFO,
        kickoff=(PARAMS.KICKOFF_SIZE, PARAMS.KICKOFF_BACKWARD),
        driver='es_spc'
    )

if PARAMS.RUN_ES_RXN:
    lmech.run_driver(
        PES_DCT, PESNUMS_LST, PARAMS.CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        PARAMS.TSK_INFO_LST,
        PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
        ene_coeff=[1.],
        vdw_params=[False, False, True],
        options=PARAMS.OPTIONS_RATE,
        etrans=lmech.etrans_lst(PARAMS),
        pst_params=PARAMS.PST_PARAMS,
        rad_rad_ts=PARAMS.RAD_RAD_TS,
        hind_inc=PARAMS.HIND_INC,
        mc_nsamp=PARAMS.MC_NSAMP0,
        temps=PARAMS.TEMPS,
        pressures=PARAMS.PRESSURES,
        multi_info=PARAMS.MULTI_INFO,
        assess_pdep=PARAMS.ASSESS_PDEP,
        kickoff=(PARAMS.KICKOFF_SIZE, PARAMS.KICKOFF_BACKWARD),
        driver='es_rxn'
    )

# if PARAMS.RUN_THERMO:
#     thermodriver.run(
#         PARAMS.TSK_INFO_LST, SPC_DCT, PARAMS.REF_MOLS,
#         PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
#         ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)

if PARAMS.RUN_RATES:
    lmech.run_driver(
        PES_DCT, PESNUMS_LST, PARAMS.CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        PARAMS.TSK_INFO_LST,
        PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
        ene_coeff=[1.],
        vdw_params=[False, False, True],
        options=PARAMS.OPTIONS_RATE,
        etrans=lmech.etrans_lst(MODEL_DCT['etransfer']),
        pst_params=PARAMS.PST_PARAMS,
        rad_rad_ts=PARAMS.RAD_RAD_TS,
        hind_inc=PARAMS.HIND_INC,
        mc_nsamp=PARAMS.MC_NSAMP0,
        temps=PARAMS.TEMPS,
        pressures=PARAMS.PRESSURES,
        multi_info=PARAMS.MULTI_INFO,
        assess_pdep=PARAMS.ASSESS_PDEP,
        driver='ktp'
    )
