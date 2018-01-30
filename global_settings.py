"""
global settings
"""
import sys
import os
from collections import OrderedDict, defaultdict
import union_find


def get_chattering_species(file_dir, atom_followed="C"):
    """
    get_chattering_species
    """
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import local_settings
        return local_settings.get_chattering_species(atom_followed)
    except IOError:
        return OrderedDict()


def get_union_find_group(file_dir, atom_followed="C"):
    """
    return union_find_groups
    """
    chattering_species = get_chattering_species(file_dir, atom_followed)

    counter = 0
    spe_idx_label = dict()
    label_spe_idx = dict()

    u_set = set()

    for _, pair_label in enumerate(chattering_species):
        idx1 = chattering_species[pair_label][0]
        idx2 = chattering_species[pair_label][1]
        if int(idx1) not in u_set:
            u_set.add(int(idx1))
            spe_idx_label[int(idx1)] = counter
            label_spe_idx[counter] = int(idx1)
            counter += 1
        if int(idx2) not in u_set:
            u_set.add(int(idx2))
            spe_idx_label[int(idx2)] = counter
            label_spe_idx[counter] = int(idx2)
            counter += 1
    print(spe_idx_label, label_spe_idx)
    wqnpc = union_find.WeightedQuickUnionWithPathCompression(len(u_set))
    for _, pair_label in enumerate(chattering_species):
        idx1 = chattering_species[pair_label][0]
        idx2 = chattering_species[pair_label][1]

        idx1 = int(spe_idx_label[int(idx1)])
        idx2 = int(spe_idx_label[int(idx2)])
        wqnpc.unite(idx1, idx2)

    # unique labels
    unique_labels = set()
    for idx, _ in enumerate(spe_idx_label):
        l_tmp = wqnpc.root(idx)
        unique_labels.add(l_tmp)

    # unique labels and their group
    unique_labels_group = defaultdict(set)
    for idx, _ in enumerate(spe_idx_label):
        l_tmp = wqnpc.root(idx)
        if l_tmp not in unique_labels_group:
            unique_labels_group[l_tmp] = set()
            unique_labels_group[l_tmp].add(label_spe_idx[idx])
        else:
            unique_labels_group[l_tmp].add(label_spe_idx[idx])
    # print(unique_labels_group)

    # species index and the big group it belongs to
    idx_group = defaultdict(set)
    for idx, _ in enumerate(spe_idx_label):
        if idx in unique_labels_group:
            idx_group[str(label_spe_idx[idx])] = unique_labels_group[idx]
        else:
            idx_group[str(label_spe_idx[idx])
                      ] = unique_labels_group[int(wqnpc.root(int(idx)))]
    print(idx_group)

    return idx_group


def get_setting(file_dir):
    """
    return global settings
    """
    setting = {}
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import local_settings
        return local_settings.get_local_settings()
    except IOError:
        return setting


def get_s_a_setting(file_dir):
    """
    return sensitivity analysis setting
    """
    setting = {}
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import local_settings
        return local_settings.get_s_a_setting()
    except IOError:
        return setting


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    get_union_find_group(FILE_DIR, atom_followed="C")
    print("test")
