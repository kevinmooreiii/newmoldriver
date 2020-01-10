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

# Parse the species input to get a dct with ALL species in mechanism
SPC_DCT = loadspc.build_spc_dct(JOB_PATH, RUN_INP_DCT['spc'])

# Parse the mechanism input and get a dct with info on PESs user request to run
PES_DCT, CHNLS_DCT = loadmech.parse_mechanism_file(
    JOB_PATH, RUN_INP_DCT['mech'], SPC_DCT, PESNUMS,
    check_stereo=RUN_INP_DCT['check_stereo'],
    sort_rxns=RUN_INP_DCT['sort_rxns'],
    rad_rad_sort=RUN_INP_DCT['rad_rad_sort'])

# Run the requested drivers: es, thermo, ktp
if RUN_JOBS_DCT['es_spc']:
    lmech.run_driver(
        PES_DCT, [PESNUMS], CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_JOBS_DCT,
        ES_TSK_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='es_spc'
    )

if RUN_JOBS_DCT['es_rxn']:
    lmech.run_driver(
        PES_DCT, [PESNUMS], CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_JOBS_DCT,
        ES_TSK_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='es_rxn'
    )

# if RUN_JOBS_DCT['thermo'] or RUN_JOBS_DCT['nasa']:
#     thermodriver.run(
#         PARAMS.TSK_INFO_LST, SPC_DCT, PARAMS.REF_MOLS,
#         PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
#         ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)

if RUN_JOBS_DCT['rates'] or RUN_JOBS_DCT['params']:
    lmech.run_driver(
        PES_DCT, [PESNUMS], CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_JOBS_DCT,
        ES_TSK_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='ktp'
    )
