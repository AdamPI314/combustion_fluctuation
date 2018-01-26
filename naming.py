"""
namespace, naming stuff
"""

import os
import sys
from shutil import copyfile
import read_write_configuration as rwc


def get_suffix(file_dir, init_spe=None, atom_followed=None, begin_t=None, end_t=None):
    """
    get suffix
    """
    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))
    suffix = ""
    if init_spe is None:
        suffix += "_S" + str(setting['pathway']['init_spe'])
    else:
        suffix += "_S" + str(init_spe)
    if atom_followed is None:
        suffix += "_" + str(setting['pathway']['atom_followed'])
    else:
        suffix += "_" + str(atom_followed)
    if begin_t is None:
        suffix += "_" + str(setting['pathway']['begin_t'])
    else:
        suffix += "_" + str(begin_t)
    if end_t is None:
        suffix += "_" + str(setting['pathway']['end_t'])
    else:
        suffix += "_" + str(end_t)

    return suffix


def copy_sohr_files(file_dir, species_path=False):
    """
    make a copy of SOHR files
    1. output/pathway_stat.csv
    2. output/pathway_name_candidate.csv
    3. output/pathway_time_candidate.csv
    4. output/pathway_name_selected.csv
    5. output/pathway_prob.csv
    6. output/pathway_AT.csv
    7. output/pathname_prob.csv
    8. output/chattering_group_info.json
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = get_suffix(file_dir)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_stat.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_stat" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix +
                         "pathway_name_selected.csv")
    f_n_2 = os.path.join(file_dir, "output",
                         prefix + "pathway_name_selected" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_prob.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_AT.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_AT" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_AT_no_IT.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_AT_no_IT" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_AT_with_SP.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_AT_with_SP" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathway_SP.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathway_SP" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix + "pathname_prob.csv")
    f_n_2 = os.path.join(file_dir, "output", prefix +
                         "pathname_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", "chattering_group_info.json")
    f_n_2 = os.path.join(file_dir, "output",
                         "chattering_group_info" + suffix + ".json")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix +
                         "pathway_name_candidate.csv")
    f_n_2 = os.path.join(file_dir, "output",
                         prefix + "pathway_name_candidate" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", prefix +
                         "pathway_time_candidate.csv")
    f_n_2 = os.path.join(file_dir, "output",
                         prefix + "pathway_time_candidate" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    return


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    copy_sohr_files(FILE_DIR)
