"""
update settings.json
"""

import os
from shutil import copy2
import read_write_configuration as rwc
import global_settings


def update_basic_setting(file_dir, g_s):
    """
    update settings.json, the basic information that's will not change for this system
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['system']['condition'] = g_s['system']['condition']
    setting['system']['initializer'] = g_s['system']['initializer']
    setting['network']['merge_chatterings'] = g_s['network']['merge_chatterings']
    setting['propagator']['primary_type'] = g_s['propagator']['primary_type']
    setting['propagator']['type'] = g_s['propagator']['type']
    setting['propagator']['sub_type'] = g_s['propagator']['sub_type']
    setting['propagator']['convert_molar_concentration_to_mole_fraction'] = g_s['propagator']['convert_molar_concentration_to_mole_fraction']
    setting['propagator']['normalize_initial_concentration'] = g_s['propagator']['normalize_initial_concentration']

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_dlsode_setting(file_dir, max_time=1.0, critical_time=0.9):
    """
    update settings.json, primarily for dlsode run and pathway generating
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['time']['critical_time'] = critical_time
    setting['time']['max_time'] = max_time

    setting['job']['job_type'] = "solve_ODEs_for_concentration_using_LSODE"
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_terminal_species_setting(file_dir, terminal_spe=None):
    """
    update settings.json, primarily for terminal species
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    t_s = []
    if terminal_spe is not None and terminal_spe is not []:
        for _, val in enumerate(terminal_spe):
            t_s.append(val)

    setting['pathway']['terminal_species'] = t_s

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
    return


def update_chattering_species_setting(file_dir, atom_followed="C"):
    """
    update settings.json, primarily for chattering species and fast reactions
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    chattering_spe = global_settings.get_chattering_species(
        file_dir, atom_followed)
    setting['pathway']['chattering_species'] = chattering_spe

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
    return


def update_spe_concentration_at_time_w2f(file_dir, tau=10.0, end_t=1.0):
    """
    update settings.json, primarily for update_spe_concentration_at_time_w2f
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['job']['job_type'] = "write_concentration_at_time_to_file"
    setting['time']['tau'] = tau
    setting['pathway']['end_t'] = end_t

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_mc_trajectory_setting(file_dir, n_traj=1000000, atom_followed="C", init_spe=114,
                                 tau=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    update settings.json, primarily for generate_pathway_running_Monte_carlo_trajectory
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    chattering_spe = global_settings.get_chattering_species(
        file_dir, atom_followed)
    setting['pathway']['chattering_species'] = chattering_spe

    setting['time']['tau'] = tau

    setting['pathway']['trajectoryNumber'] = n_traj
    setting['pathway']['atom_followed'] = atom_followed
    setting['pathway']['init_spe'] = init_spe
    setting['pathway']['begin_t'] = begin_t
    setting['pathway']['end_t'] = end_t

    if species_path is True:
        setting['job']['job_type'] = "generate_species_pathway_running_Monte_carlo_trajectory"
    else:
        setting['job']['job_type'] = "generate_pathway_running_Monte_carlo_trajectory"

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
    return


def update_eval_path_integral(file_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114,
                              tau=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    update settings.json, primarily for evaluate path integral
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    chattering_spe = global_settings.get_chattering_species(
        file_dir, atom_followed)
    setting['pathway']['chattering_species'] = chattering_spe

    if species_path is True:
        setting['job']['job_type'] = "evaluate_species_path_integral_over_time"
    else:
        setting['job']['job_type'] = "evaluate_path_integral_over_time"
    setting['pathway']['topN'] = [top_n]
    setting['pathway']['trajectoryNumber'] = n_traj
    setting['pathway']['atom_followed'] = atom_followed
    setting['pathway']['init_spe'] = init_spe

    setting['time']['tau'] = tau
    setting['pathway']['begin_t'] = begin_t
    setting['pathway']['end_t'] = end_t

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
