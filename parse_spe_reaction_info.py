"""
parse species and reaction information
"""

import re
import os
import sys
import json
import pandas as pd
import numpy as np
import read_write_configuration as rwc
# import time


def parse_spe_info(file_dir):
    """
    parse species info from file= "os.path.join(file_dir, "output", "species_labelling.csv")"
    """
    f_n = os.path.join(file_dir, "input", "species_labelling.csv")
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')

    matched_str = [re.findall(r"(\d+)\t-->\t([\w|\-|(|)]+)", line)[0]
                   for line in line_content]
    matched_str_reverse = [(x[1], x[0]) for x in matched_str]

    spe_ind_name_dict = dict(matched_str)
    spe_name_ind_dict = dict(matched_str_reverse)

    return spe_ind_name_dict, spe_name_ind_dict


def parse_species_pair_reaction(file_dir):
    """
    parse species pairs and associated reactions, coefficients
    """
    f_n = os.path.join(file_dir, "input", "species_pairs_reactions_coefs.json")

    s_p_r_c = rwc.read_configuration(f_n)
    return s_p_r_c


def read_spe_composition(f_n):
    """
    read species composition
    """
    with open(f_n, 'r') as f_h:
        data = json.load(f_h)
    return data


def parse_reaction_and_its_index(file_dir):
    """
    parse reaction info from file= "os.path.join(file_dir, "input", "reaction_labelling.csv")"
    """
    f_n = os.path.join(file_dir, "input", "reaction_labelling.csv")
    # load data
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')
    matched_tmp = [re.findall(r"([\d]+)\s+([\-\d]+)\s+([\w\(\)\-\_,\+]+\<?={1}\>?[\w\(\)\-\_,\+]+)", line)
                   for line in line_content]
    matched_ind1_ind2_str = [x[0] for x in matched_tmp if len(x) != 0]
    # map the new old reaction index
    new_old_index_dict = dict()
    for _, val in enumerate(matched_ind1_ind2_str):
        new_old_index_dict.update(
            {val[0]: str(val[1])})

    # reactant arrow product
    reactant_product = [re.findall(r"([\w|+|(|)|\-|\_|,]+)[=|>|<]+([\w|+|(|)|\-|\_|,]+)",
                                   ind1_ind2_reaction[2])[0]
                        for ind1_ind2_reaction in matched_ind1_ind2_str]
    reactant = [x[0] for x in reactant_product]
    product = [x[1] for x in reactant_product]
    # map reaction new reaction label and the exact reaction
    new_ind_reaction_dict = dict()
    for i in range(len(reactant_product)):
        if int(matched_ind1_ind2_str[i][1]) > 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: reactant[i] + '=>' + product[i]})
        elif int(matched_ind1_ind2_str[i][1]) < 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: product[i] + '=>' + reactant[i]})
    return new_old_index_dict, new_ind_reaction_dict


def reaction_name_to_real_reaction(new_ind_reaction_dict, pathway_name):
    """
    converted reaction name to their reaction format instead of index
    """
    matched_reaction = re.findall(r"R([-]?\d+)", pathway_name)
    # only reactions
    str_t = '['
    for _, val in enumerate(matched_reaction):
        if '-' not in val:
            str_t += new_ind_reaction_dict[val]
        else:
            str_t += "<-chattering->"
    str_t += ']'
    return str_t


def pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, pathway_name):
    """
    converted path to their real species name and reaction format instead of index
    """
    # always starts from species
    str_t = ""
    matched_s_r = re.findall(r"S\d+(?:R[-]?\d+)?", pathway_name)
    for idx, val in enumerate(matched_s_r):
        m_s = re.findall(r"S(\d+)", val)
        m_r = re.findall(r"R([-]?\d+)", val)

        m_s_idx = m_s[0]
        str_t += '[' + spe_ind_name_dict[m_s_idx] + ']'
        if (len(m_r) == 0):
            if idx != len(matched_s_r) - 1:
                str_t += "-->"
        elif len(m_r) > 0:
            m_r_idx = m_r[0]
            if '-' not in m_r_idx:
                str_t += new_ind_reaction_dict[m_r_idx]
                str_t += "-->"
            else:
                str_t += "<-chattering->"

    return str_t


def symbolic_path_2_real_path(file_dir, f_n_p, f_n_p_out, top_n=50, end_s_idx=None, max_rows=5000):
    """
    read species and reaction info,
    convert path info into real species and reaction instead of index and write to file
    """

    # load path data
    if end_s_idx is None or end_s_idx is []:
        path_data = pd.read_csv(f_n_p, names=['path', 'prob'], nrows=top_n + 1)
    elif end_s_idx is not None and end_s_idx is not []:
        n_spe = len(end_s_idx)
        path_data = pd.read_csv(
            f_n_p, names=['path', 'prob'], nrows=top_n * n_spe + max_rows)

    total_prob = sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    # renormalize
    path_data['prob'] /= total_prob

    # filter
    if end_s_idx is not None and end_s_idx is not []:
        end_spe_str = ['S' + str(x) for x in end_s_idx]
        end_spe_tuple = tuple(end_spe_str)
        path_data = path_data[path_data['path'].str.endswith(end_spe_tuple)]

    # load spe and reaction info
    spe_ind_name_dict, _ = parse_spe_info(file_dir)
    _, new_ind_reaction_dict = parse_reaction_and_its_index(file_dir)

    # convert species reaction index to real species and reactions
    path_data['path'] = path_data['path'].apply(
        lambda x: pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, x))

    # write to file
    path_data[0:top_n].to_csv(f_n_p_out, header=False,
                              index=False, sep=',', columns=['path', 'prob'])


def parse_reaction_net_reactant(file_dir):
    """
    return a dict of "species": number based on reaction reactant
    """
    f_n = os.path.join(file_dir, "input", "reaction_information.json")

    data = rwc.read_configuration(f_n)

    net_reactant = {}
    for _, r_idx in enumerate(data):
        entry = {}
        for val1 in data[r_idx]['net_reactant']:
            entry.update({data[r_idx]['net_reactant'][val1]['species_index']:
                          data[r_idx]['net_reactant'][val1]['coefficient']})
        net_reactant.update({r_idx: entry})

    return net_reactant


def parse_reaction_net_product(file_dir):
    """
    return a dict of "species": number based on reaction product
    """
    f_n = os.path.join(file_dir, "input", "reaction_information.json")

    data = rwc.read_configuration(f_n)

    net_product = {}
    for _, r_idx in enumerate(data):
        entry = {}
        for val1 in data[r_idx]['net_product']:
            entry.update({data[r_idx]['net_product'][val1]['species_index']:
                          data[r_idx]['net_product'][val1]['coefficient']})
        net_product.update({r_idx: entry})

    return net_product


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    parse_reaction_net_product(FILE_DIR)
