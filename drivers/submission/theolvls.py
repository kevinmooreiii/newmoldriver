""" Electronic structure parameters; code, method, basis, convergence control
"""

ES_DCT = {
    'lvl_wbs': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31g*',
        },
    'lvl_wbm': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': '6-31+g*',
        },
    'lvl_wbt': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'wb97xd',
        'basis': 'cc-pvtz',
        },
    'lvl_m06s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31g*',
        },
    'lvl_m06m': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': '6-31+g*',
        },
    'lvl_m06t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'm062x',
        'basis': 'cc-pvtz',
        },
    'lvl_b2d': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvdz',
        },
    'lvl_b2t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvtz',
        },
    'lvl_b2q': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b2plypd3',
        'basis': 'cc-pvqz',
        },
    'lvl_b3s': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        },
    'lvl_b3t': {
        'orb_res': 'RU', 'program': 'gaussian09', 'method': 'b3lyp',
        'basis': '6-31g*',
        },
    'cc_lvl_d': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvdz',
        },
    'cc_lvl_t': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvtz',
        },
    'cc_lvl_q': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)',
        'basis': 'cc-pvqz',
        },
    'cc_lvl_df': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvdz-f12',
        },
    'cc_lvl_tf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvtz-f12',
        },
    'cc_lvl_qf': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'ccsd(t)-f12',
        'basis': 'cc-pvqz-f12',
        },
    'mlvl_cas_dz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvdz',
        },
    'mlvl_cas_tz': {
        'orb_res': 'RR', 'program': 'molpro2015', 'method': 'caspt2',
        'basis': 'cc-pvtz',
        }
}
