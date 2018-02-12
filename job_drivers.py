"""
Job drivers
"""

import subprocess
import os
import sys
import update_settings as us
import parse_spe_reaction_info as psri
import prepare_path_name_time as ppnt
import naming
import trajectory
import global_settings


def update_terminal_species_setting(data_dir, terminal_spe=None):
    """
    update settings.json, primarily for terminal species
    """
    us.update_terminal_species_setting(data_dir, terminal_spe=terminal_spe)


def update_chattering_species_setting(data_dir, atom_followed="C"):
    """
    update settings.json, primarily for chattering species and fast reactions
    """
    us.update_chattering_species_setting(data_dir, atom_followed)


def copy_sohr_files(data_dir, species_path=False):
    """
    copy SOHR files from C++ routine
    """
    naming.copy_sohr_files(data_dir, species_path=species_path)


def symbolic_path_2_real_path(data_dir, top_n=50, flag="", end_s_idx=None, species_path=False, max_rows=5000):
    """
    convert symbolic pathway to real pathway with real species name and real reaction name
    flag indicates a specific job, for example, pathway end time = 1.0, the j-th run,
    any unique symbol shall work
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"

    if flag == "":
        out_file_name = prefix + "pathname_prob.csv"
    else:
        out_file_name = prefix + "pathname_prob_" + str(flag) + ".csv"

    path_stat_fn = prefix + "pathway_stat.csv"

    psri.symbolic_path_2_real_path(
        data_dir,
        os.path.join(
            data_dir, "output", path_stat_fn),
        os.path.join(
            data_dir, "output", out_file_name),
        top_n, end_s_idx, max_rows=max_rows)


# path from file
def symbolic_path_2_real_path_pff(data_dir, fn):
    """
    convert symbolic pathway to real pathway with real species name and real reaction name
    """
    out_fn = fn[0:-4] + "_real_path" + ".csv"

    psri.symbolic_path_2_real_path(
        data_dir,
        os.path.join(
            data_dir, "output", fn),
        os.path.join(
            data_dir, "output", out_fn),
        10000000, None)


def delete_non_dlsode_files(data_dir):
    """
    delete none dlsode files
    """
    os.chdir(data_dir)
    cmd = ["find", "./output", "-type", "f",
           "!", "-name", "*dlsode*", "-delete"]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


def make_run(src_dir, data_dir):
    """
    make run
    """
    os.chdir(src_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


def make_run_timeout(src_dir, data_dir, timeout=None):
    """
    make run with timeout functionaliry
    https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    """
    os.chdir(src_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        subprocess.check_call(
            cmd, stdout=out_file, stderr=error_file,
            timeout=timeout)
        status = 0
    except subprocess.TimeoutExpired:
        status = -1
    except subprocess.CalledProcessError:
        status = -1

    out_file.close()
    error_file.close()

    if status == 0:
        return True
    else:
        return False


def run_dlsode(src_dir, data_dir, max_time, critical_time):
    """
    Run dlsode
    """
    os.chdir(src_dir)
    us.update_dlsode_setting(data_dir, max_time, critical_time)
    make_run(src_dir, data_dir)


def spe_concentration_at_time_w2f(src_dir, data_dir, tau=10.0, end_t=1.0):
    """
    write species concentration at a time to file
    """
    os.chdir(src_dir)
    us.update_spe_concentration_at_time_w2f(data_dir, tau=tau, end_t=end_t)
    make_run(src_dir, data_dir)


def run_mc_trajectory(src_dir, data_dir, n_traj=1000000, atom_followed="C", init_spe=114,
                      tau=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    Run mc trajectory
    """
    os.chdir(src_dir)
    us.update_mc_trajectory_setting(
        data_dir, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe,
        tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)
    make_run(src_dir, data_dir)


def evaluate_pathway_probability(src_dir, data_dir, top_n=5, num_t=1, flag="", n_traj=10000,
                                 atom_followed="C", init_spe=114, traj_max_t=100.0,
                                 tau=10.0, begin_t=0.0, end_t=1.0, top_n_s=10,
                                 spe_oriented=True, end_s_idx=None, species_path=False):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(src_dir)

    if spe_oriented is True:
        us.update_eval_path_integral(
            data_dir, top_n=top_n * top_n_s, n_traj=n_traj,
            atom_followed=atom_followed, init_spe=init_spe,
            tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)

        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                data_dir, exclude=None, top_n=top_n_s, traj_max_t=traj_max_t,
                tau=tau, end_t=end_t, tag="M", atoms=[atom_followed])
        ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path)
        ppnt.prepare_pathway_time(
            data_dir, top_n=top_n * top_n_s, num=num_t, flag=flag,
            begin_t=begin_t, end_t=end_t, species_path=species_path)
    else:
        us.update_eval_path_integral(
            data_dir, top_n=top_n, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe,
            tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)
        ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path)
        ppnt.prepare_pathway_time(
            data_dir, top_n=top_n, num=num_t, flag=flag, begin_t=begin_t, end_t=end_t, species_path=species_path)

    make_run(src_dir, data_dir)


def send_email(data_dir):
    """
    send email to elliot.srbai@gmail.com
    """
    os.chdir(data_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
           "-u", "RUNNING JOB", "-m", "JOB FINISHED." + "\n" + data_dir,
           "-a", os.path.join(data_dir, "output", "output_all.txt")]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


if __name__ == '__main__':
    SRC_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    G_S = global_settings.get_setting(DATA_DIR)
    # symbolic_path_2_real_path_pff(DATA_DIR, "heuristic_pathname_O_10_10_3.csv")
    spe_concentration_at_time_w2f(
        SRC_DIR, DATA_DIR, tau=G_S['tau'], end_t=0.00022854295)
