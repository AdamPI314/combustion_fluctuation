"""
main driver
"""

import os
import sys
import time

import update_settings as us
import job_drivers
import global_settings

if __name__ == '__main__':
    TIME_I = time.time()

    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)

    G_S = global_settings.get_setting(DATA_DIR)
    us.update_basic_setting(DATA_DIR, G_S)

    # run dlosde
    job_drivers.run_dlsode(DATA_DIR, G_S['traj_max_t'], G_S['traj_critical_t'])

    # update terminal species
    job_drivers.update_terminal_species_setting(DATA_DIR, G_S['terminal_spe'])

    # update chattering species and fast reactions
    job_drivers.update_chattering_species_setting(
        DATA_DIR, G_S['atom_f'])

    # write specie concentration at a time to file
    # job_drivers.spe_concentration_at_time_w2f(
    #     DATA_DIR, tau=G_S['tau'], end_t=G_S['end_t'])
    job_drivers.spe_concentration_at_time_w2f(
        DATA_DIR, tau=G_S['tau'], end_t=0.00228542952)

    # # run monte carlo trajectory
    # job_drivers.run_mc_trajectory(
    #     DATA_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
    #     init_spe=G_S['init_s'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['mc_t'],
    #     species_path=G_S['species_path'])

    # # evaluate path integral-->pathway probability
    # job_drivers.evaluate_pathway_probability(
    #     DATA_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # # convert symbolic pathway to real pathway
    # # with real species names and real reaction expression
    # job_drivers.symbolic_path_2_real_path(DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #                                       end_s_idx=None, species_path=G_S['species_path'])

    # # copy SOHR/C++ routine files
    # job_drivers.copy_sohr_files(DATA_DIR, species_path=G_S['species_path'])

    # # send email
    # job_drivers.send_email(DATA_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
