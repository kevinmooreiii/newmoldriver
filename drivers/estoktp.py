"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import os

# New libs
from lib.load import mechanism as loadmech
from lib.load import species as loadspc
import lib.filesystem.build as lfs
from lib.submission import read_dat
import mech as lmech

# Set runtime options based on user input
INPUT = lmech.get_user_input()
DATA_PATH, MECH_PATH, MECH_TYPE, MECH_FILE, PESNUMS, CHANNELS = INPUT
GEOM_DCT = lmech.build_geom_dct(DATA_PATH)
PARAMS = read_dat.params(os.path.join(MECH_PATH, 'params.dat'))

# Prepare run-save filesystem for job
lfs.prefix_filesystem(PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX)

# Read the input files
SPC_STR = open(os.path.join(MECH_PATH, 'species.csv')).read()
MECH_STR = open(os.path.join(MECH_PATH, 'mechanism.txt')).read()

# Parse the species input
SPC_DCT = loadspc.build_spc_dct(SPC_STR, 'csv')
MOD_SPC_DCT = loadspc.modify_spc_dct(SPC_DCT, GEOM_DCT, PARAMS.HIND_INC)

# Parse the mechanism input
PES_DCT = loadmech.parse_mechanism_file(
    MECH_STR, MECH_TYPE, MOD_SPC_DCT,
    check_stereo=PARAMS.CHECK_STEREO,
    sort_rxns=PARAMS.SORT_RXNS,
    rad_rad_sort=PARAMS.RAD_RAD_SORT)

# # Run thermo if desired
# if RUN_THERMO:
#     SPC_QUEUE = list(SPC_NAMES)
#     thermodriver.driver.run(
#         PARAMS.TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, PARAMS.REF_MOLS,
#         PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
#         ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)
#
# Set the channels on the PESs to run
loadmech.print_pes_channels(PES_DCT)
PESNUMS_LST = loadmech.get_pes_nums(PES_DCT, PARAMS.PESNUMS)
CONN_CHNLS_LST = loadmech.determine_connected_pes_channels(
    PES_DCT, PESNUMS_LST)

# Run the es driver
if PARAMS.RUN_ES:
    lmech.run_driver(
        PES_DCT, PESNUMS_LST, PARAMS.CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
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
        driver='es_spc'
    )
    lmech.run_driver(
        PES_DCT, PESNUMS_LST, PARAMS.CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
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
        driver='es_rxn'
    )

# Run the rates driver
if PARAMS.RUN_RATES:
    lmech.run_driver(
        PES_DCT, PESNUMS_LST, PARAMS.CHANNELS, CONN_CHNLS_LST,
        MOD_SPC_DCT, {},
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
        driver='ktp'
    )
