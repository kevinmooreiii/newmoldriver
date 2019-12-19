"""
   Driver to parse and sort the mechanism input files and
   launch the desired drivers
"""

import lib.drivers.mech as lmechdriver
import lib.filesystem.build as lfs


# Set runtime options based on user input
INPUT = lmechdriver.get_user_input()
DATA_PATH, MECH_PATH, MECH_TYPE, MECH_FILE, PESNUMS, CHANNELS = INPUT
GEOM_DCT = lmechdriver.build_geom_dct(DATA_PATH)
PARAMS = lmechdriver.read_params_file(MECH_PATH)

# Prepare run-save filesystem for job
lfs.build_filesystem(PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX)

# Read the input files
RXN_INFO = lmechdriver.parse_mechanism_file(
    MECH_TYPE, MECH_PATH, MECH_FILE,
    PARAMS.CHECK_STEREO, PARAMS.SORT_RXNS, PARAMS.RAD_RAD_SORT)
SPC_DCT, RCT_NAMES, PRD_NAMES, RXN_NAME, FORM_STRS = RXN_INFO

# Update the species information
lmechdriver.update_spc_dct(SPC_DCT, GEOM_DCT, PARAMS.HIND_INC)

# # Run thermo if desired
# if RUN_THERMO:
#     SPC_QUEUE = list(SPC_NAMES)
#     thermodriver.driver.run(
#         PARAMS.TSK_INFO_LST, ES_DCT, SPC_DCT, SPC_QUEUE, PARAMS.REF_MOLS,
#         PARAMS.RUN_PREFIX, PARAMS.SAVE_PREFIX,
#         ene_coeff=PARAMS.ENE_COEFF, options=PARAMS.OPTIONS_THERMO)
#
# Set the channels on the PESs to run
PES_DCT = lmechdriver.build_pes_dct(
    FORM_STRS, RCT_NAMES, PRD_NAMES, RXN_NAME)
lmechdriver.print_pes_channels(PES_DCT)
PESNUMS = lmechdriver.set_pes_nums(
    PES_DCT, PARAMS.PESNUMS)
PES = lmechdriver.determine_connected_channels(
    PES_DCT, PARAMS.PESNUMS, PARAMS.CHANNELS)
PESCHNS, CONNCHNS, PES_RCT_NAMES_LST, PES_PRD_NAMES_LST, PES_RXN_NAME_LST = PES

# Run the reactions driver (es, ktp, etc)
if PARAMS.RUN_RATES:
    lmechdriver.run_reactions_driver(
        SPC_DCT, PES_DCT, CLA_DCT,
        PESNUMS, CONNCHNS, PESCHNS,
        PES_RCT_NAMES_LST,
        PES_PRD_NAMES_LST,
        PES_RXN_NAME_LST,
        PARAMS.HIND_INC,
        PARAMS.EXP_FACTOR,
        PARAMS.EXP_POWER,
        PARAMS.EXP_CUTOFF,
        PARAMS.EPS1,
        PARAMS.EPS2,
        PARAMS.SIG1,
        PARAMS.SIG2,
        PARAMS.MASS1,
        PARAMS.TSK_INFO_LST,
        RCT_NAMES,
        PRD_NAMES,
        PARAMS.RUN_PREFIX,
        PARAMS.SAVE_PREFIX,
        PARAMS.ENE_COEFF,
        PARAMS.OPTIONS_RATE,
        PARAMS.PST_PARAMS,
        PARAMS.RAD_RAD_TS)
print('EXITING')
