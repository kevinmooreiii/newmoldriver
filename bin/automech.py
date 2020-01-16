"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import sys
from drivers import mech as lmech
# from drivers import thermodriver
from lib.load import run as loadrun
from lib.load import theory as loadthy
from lib.load import model as loadmodel
from lib.load import mechanism as loadmech
from lib.load import species as loadspc


# Set runtime options based on user input
JOB_PATH = sys.argv[1]
# PARAMS = read_dat.params(os.path.join(JOB_PATH, 'inp/params.dat'))

# Parse the run input
RUN_INP_DCT = loadrun.build_run_inp_dct(JOB_PATH)
RUN_LST = loadrun.objects_lst(JOB_PATH)
print(RUN_LST)
import sys
sys.exit()
[PESNUMS, CHANNELS, MODEL] = RUN_LST[0]

RUN_OPTIONS_DCT = loadrun.build_run_glob_opts_dct(JOB_PATH)
RUN_JOBS_LST = loadrun.build_run_jobs_lst(JOB_PATH)
RUN_ES_TSK_STR = loadrun.read_es_tsks(JOB_PATH)

# Parse the theory input
THY_DCT = loadthy.build_thy_dct(JOB_PATH)

# Parse the model input
ALL_MODEL_DCT = loadmodel.read_models_sections(JOB_PATH)
MODEL_DCT = ALL_MODEL_DCT[MODEL]

# Parse the species input to get a dct with ALL species in mechanism
SPC_DCT = loadspc.build_spc_dct(JOB_PATH, RUN_INP_DCT['spc'])

# Parse the mechanism input and get a dct with info on PESs user request to run
PES_DCT, CHNLS_DCT = loadmech.parse_mechanism_file(
    JOB_PATH, RUN_INP_DCT['mech'], SPC_DCT, [PESNUMS],
    sort_rxns=RUN_INP_DCT['sort_rxns'])

# Do some extra work to prepare the info to pass to the drivers
RUN_ES_TSK_LST = loadrun.build_run_es_tsks_lst(RUN_ES_TSK_STR, MODEL_DCT, THY_DCT)
print('first RUN_ES_TSK_LST')
print(RUN_ES_TSK_LST)

# Print stuff for test
print('test auto input')
print('\nrun inp dct')
print(RUN_INP_DCT)
print('\nrun options dct')
print(RUN_OPTIONS_DCT)
print('\nrun jobs lst')
print(RUN_JOBS_LST)
print('\nrun es tsk lst')
print(RUN_ES_TSK_LST)
print('\ntheory dct')
print(THY_DCT)
print('\nall model dct')
print(ALL_MODEL_DCT)
print('\nspc dct')
print(SPC_DCT)
print('\npes dct')
print(PES_DCT)
print('\nchnls dct')
print(CHNLS_DCT)

# Run the requested drivers: es, thermo, ktp
print('\n\n')
print('RUNNING ES DRIVER FOR SPC')
if 'es_spc' in RUN_JOBS_LST:
    lmech.run_driver(
        PES_DCT, CHNLS_DCT,
        SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_ES_TSK_LST,
        RUN_JOBS_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='es_spc'
    )

print('RUNNING ES DRIVER FOR RXNS')
if 'es_rxn' in RUN_JOBS_LST:
    lmech.run_driver(
        PES_DCT, CHNLS_DCT,
        SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_ES_TSK_LST,
        RUN_JOBS_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='es_rxn'
    )

# if RUN_JOBS_LST['thermo'] or RUN_JOBS_LST['nasa']:
#     thermodriver.run(
#         PARAMS.TSK_INFO_LST, SPC_DCT, PARAMS.REF_MOLS,
#         PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
#         ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)

print('RUNNING KTP DRIVER')
if 'rates' in RUN_JOBS_LST or 'params' in RUN_JOBS_LST:
    lmech.run_driver(
        PES_DCT, CHNLS_DCT,
        SPC_DCT, {},
        THY_DCT,
        MODEL_DCT,
        RUN_ES_TSK_LST,
        RUN_JOBS_LST,
        RUN_INP_DCT,
        RUN_OPTIONS_DCT,
        driver='ktp'
    )
