"""
read write configurations
"""

import json


def read_configuration(f_n):
    """
    read configuration
    """
    with open(f_n, 'r') as fp1:
        s_d = json.load(fp1)
    return s_d


def write_configuration(setting_data, f_n):
    """
    write configuration
    """
    with open(f_n, 'w') as fp2:
        json.dump(setting_data, fp2, sort_keys=False, indent=4)
