"""
followed atom scheme, either natural atoms, or hypothesized atoms
"""

import os
import sys
import time
from shutil import copy2
import read_write_configuration as rwc


def spe_composition_2_atom_scheme(data_dir):
    """
    convert species grouped atom scheme, which refers to file named
    "spe_composition.json" generated from cantera to a new file named
    "atom_scheme_base.json"
    """
    spe_comp = rwc.read_configuration(os.path.join(
        data_dir, "input", "spe_composition.json"))

    atom_scheme = {}
    for _, s_1 in enumerate(spe_comp):
        # print(s_1)
        for atom_1 in spe_comp[s_1]:
            if atom_1 not in atom_scheme:
                atom_scheme[atom_1] = {str(s_1): spe_comp[s_1][atom_1]}
            else:
                atom_scheme[atom_1].update({str(s_1): spe_comp[s_1][atom_1]})

    fn0 = os.path.join(data_dir, "input", "atom_scheme_base_backup.json")
    fn1 = os.path.join(data_dir, "input", "atom_scheme_base.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


def spe_information_2_atom_scheme(data_dir):
    """
    convert species information
    "species_information.json" to a new file named
    "atom_scheme_base.json"
    """
    spe_comp = rwc.read_configuration(os.path.join(
        data_dir, "input", "species_information.json"))

    atom_scheme = {}
    for _, spe_idx in enumerate(spe_comp):
        s_1 = spe_comp[spe_idx]["name"]
        for atom_1 in spe_comp[spe_idx]["spe_composition"]:
            atom_number = spe_comp[spe_idx]["spe_composition"][atom_1]
            if atom_number == "0":
                continue
            if atom_1 not in atom_scheme:
                atom_scheme[atom_1] = {str(s_1): float(atom_number)}
            else:
                atom_scheme[atom_1].update({str(s_1): float(atom_number)})

    fn0 = os.path.join(data_dir, "input", "atom_scheme_base_backup.json")
    fn1 = os.path.join(data_dir, "input", "atom_scheme_base.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


def update_a_atom_entry(data_dir, source_atoms=None, entry_name="HA4", number=1.0):
    """
    update a atom entry based on atoms list
    """
    if source_atoms is None or source_atoms is []:
        return

    f_n_as = os.path.join(data_dir, "input", "atom_scheme.json")
    atom_scheme = rwc.read_configuration(f_n_as)

    new_entry = {}
    for atom in source_atoms:
        if atom not in atom_scheme:
            continue
        for key in atom_scheme[atom]:
            if key not in new_entry:
                new_entry.update({key: number})

    atom_scheme.update({entry_name: new_entry})

    fn0 = os.path.join(data_dir, "input", "atom_scheme_base_backup.json")
    fn1 = os.path.join(data_dir, "input", "atom_scheme_base.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


def atom_scheme_set_atom_number(data_dir, followed_atom="C", number=1.0):
    """
    modify "atom_scheme_base.json", change atom number to number
    """
    fn1 = os.path.join(data_dir, "input", "atom_scheme_base.json")

    atom_scheme = rwc.read_configuration(fn1)

    for key in atom_scheme[followed_atom]:
        atom_scheme[followed_atom][key] = number

    fn0 = os.path.join(data_dir, "input", "atom_scheme_base_backup.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


def get_atom_scheme(data_dir):
    """
    read "atom_scheme.json" and return a dictionary
    """

    f_n_as = os.path.join(data_dir, "input", "atom_scheme.json")
    atom_scheme = rwc.read_configuration(f_n_as)

    return atom_scheme


if __name__ == '__main__':
    INIT_TIME = time.time()

    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)

    # spe_composition_2_atom_scheme(DATA_DIR)
    # spe_information_2_atom_scheme(DATA_DIR)
    # atom_scheme_set_atom_number(DATA_DIR, followed_atom="HA3", number=1.0)
    # update_a_atom_entry(DATA_DIR, source_atoms=["C", "O", "H"], entry_name="HA4", number=1.0)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
