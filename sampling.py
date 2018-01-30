"""
dealing with sampling things
"""

import os
import sys
import copy
import random
import numpy as np
import global_settings
import update_settings as us
import job_drivers


def get_uniform_uncertainties(n_size=3, value=3.0):
    """
    return uniform uncertainties
    """
    result = np.ones(n_size)
    for i, _ in enumerate(result):
        result[i] = float(value)
    return result


def get_random_coef(uniform_uncertainties=None):
    """
    generate a vector of random numbers in range of 0-1
    here we constraint the size >= 1
    """
    if uniform_uncertainties is None:
        return
    n_size = len(uniform_uncertainties)
    if n_size <= 1:
        return
    result = np.ones(n_size)
    for i, _ in enumerate(result):
        result[i] = random.uniform(
            1.0 / uniform_uncertainties[i], uniform_uncertainties[i])
    return result


def run_a_sample(file_dir):
    """
    run a monte carlo sample with some size
    """
    # sensitivity analysis settings
    s_a_s = global_settings.get_s_a_setting(file_dir)

    # global k file name, all together
    g_k_f_n = os.path.join(file_dir, "output", "k_global.csv")
    if os.path.isfile(g_k_f_n):
        os.remove(g_k_f_n)
    # global target file name, all together, target is ignition delay time (ign) here
    g_t_f_n = os.path.join(file_dir, "output", "ign_global.csv")
    if os.path.isfile(g_t_f_n):
        os.remove(g_t_f_n)
    # local target file name
    l_t_f_n = os.path.join(file_dir, "output", "ign_local.csv")

    u_u = get_uniform_uncertainties(
        s_a_s['n_dim'], s_a_s['default_uncertainty'])
    # save constant uncertainty to file
    f_n_u_const = os.path.join(file_dir, "output", "uncertainties_const.csv")
    np.savetxt(f_n_u_const, u_u, fmt='%.18e',  delimiter=',', newline='\n')

    for _ in range(s_a_s['n_run']):
        r_c = get_random_coef(uniform_uncertainties=u_u)

        spe_idx_conc = copy.deepcopy(s_a_s['spe_idx_conc'])
        print(spe_idx_conc)
        for s_i in spe_idx_conc:
            if int(s_i) >= 0 and int(s_i) < len(r_c):
                spe_idx_conc[s_i] *= r_c[int(s_i)]

        us.update_s_a_setting(file_dir,
                              init_temp=s_a_s['init_temp'],
                              critical_temp=s_a_s['critical_temp'],
                              target_temp=s_a_s['target_temp'],
                              end_temp=s_a_s['end_temp'],
                              spe_idx_conc=spe_idx_conc)

        flag = job_drivers.make_run_timeout(file_dir, timeout=30)

        # local target time
        local_t_t = np.loadtxt(l_t_f_n, dtype=float, delimiter=',')
        local_t_t = [local_t_t]

        # is successfully run a sample, save to file
        if flag is True:
            r_c = r_c.reshape((1, len(r_c)))
            with open(g_k_f_n, 'ab') as f_handler:
                np.savetxt(f_handler, r_c, fmt='%.18e',
                           delimiter=',', newline='\n')
            with open(g_t_f_n, 'ab') as f_handler:
                np.savetxt(f_handler, local_t_t, fmt='%.18e',
                           delimiter=',', newline='\n')


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    run_a_sample(FILE_DIR)
