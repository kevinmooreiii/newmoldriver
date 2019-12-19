""" Electronic structure parameters; code, method, basis, convergence control
"""

# new levels
import lib.filesystem.inf.get_thy_info

ES_DCT = {
    'lvl_wbs': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_wbm': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31+g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_wbt': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06m': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31+g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_m06t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2d': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b2q': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvqz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b3s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'lvl_b3t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_d': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_t': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_q': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvqz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_df': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvdz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_tf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvtz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'cc_lvl_qf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvqz-f12',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'mlvl_cas_dz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvdz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        },
    'mlvl_cas_tz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvtz',
        'mc_nsamp': PARAMS.MC_NSAMP0
        }
}


def get_es_info(method, nsamp):
    """
    Turn es dictionary in theory info array
    """
    if key == 'input':
        ret = ['input_geom', None, None, None]
    else:
        ret = scripts.es.get_thy_info(es_dct[key])
    return ret
