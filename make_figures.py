"""
plot jobs
"""
import os
import sys
from copy import deepcopy
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.ticker import FormatStrFormatter

import parse_spe_reaction_info as psri
import global_settings
from tools import get_colors_markers_linestyles
import trajectory


def plot_concentrations(file_dir, spe_idx=None, tau=10.0, end_t=1.0, tag="fraction", exclude_names=None,
                        renormalization=True, semilogy=False, hasTemp=True):
    """
    plot concentrations give species index list, if exclude is not None, means we are going
    to renormalize the molelar fraction
    """
    if exclude_names is None:
        exclude_names = []

    spe_idx_tmp = deepcopy(spe_idx)
    if spe_idx_tmp is None:
        spe_idx_tmp = [0]

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)

    if hasTemp is True:
        s_idx_n["-1"] = "Temp"
        spe_idx_tmp.append(-1)

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    conc = trajectory.get_normalized_concentration(
        file_dir, tag=tag, exclude_names=exclude_names, renormalization=renormalization)

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx_tmp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            if semilogy is True:
                a_x_left.semilogy(time[0:end_point], conc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            else:
                a_x_left.plot(time[0:end_point], conc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                              color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.set_xlim([0, tau * end_t])
    a_x_left.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    a_x_left.set_xlabel("Time/sec")
    a_x_left.set_ylabel("[X]")

    if hasTemp is True:
        leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
        leg_right.get_frame().set_alpha(0.7)
        a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx_tmp)
    # plt.title(s_n_str)

    fig.savefig(os.path.join(file_dir, "output",
                             "trajectory_" + s_n_str + ".jpg"), dpi=500)
    plt.close()


def plot_spe_drc(file_dir, spe_idx=None, tau=10.0, end_t=1.0, tag="fraction", reciprocal=False):
    """
    plot species destruction rate constant, give species index list
    """
    spe_idx_tmp = deepcopy(spe_idx)
    if spe_idx_tmp is None:
        spe_idx_tmp = [0]

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)
    s_idx_n["-1"] = "Temp"

    spe_idx_tmp.append(-1)

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    spe_drc = np.loadtxt(os.path.join(file_dir, "output",
                                      "drc_dlsode_" + str(tag) + ".csv"), delimiter=",")
    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx_tmp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point], markevery=delta_n,
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            if reciprocal is False:
                a_x_left.semilogy(time[0:end_point], spe_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            else:
                a_x_left.semilogy(time[0:end_point], 1.0 / spe_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    if reciprocal is False:
        leg_left = a_x_left.legend(loc=9, fancybox=True, prop={'size': 10.0})
    else:
        leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})

    leg_right = a_x_right.legend(loc=4, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.set_xlim([0.05, time[end_point]])

    a_x_left.set_xlabel("time/s")
    if reciprocal is False:
        a_x_left.set_ylabel("k/s$^{-1}$")
    else:
        a_x_left.set_ylabel("k$^{-1}/s$")

    a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx_tmp)
    # plt.title(s_n_str)

    if reciprocal is False:
        fig.savefig(os.path.join(file_dir, "output",
                                 "spe_drc_" + s_n_str + ".jpg"), dpi=500)
    else:
        fig.savefig(os.path.join(file_dir, "output",
                                 "spe_drc_reciprocal_" + s_n_str + ".jpg"), dpi=500)

    plt.close()


def plot_reaction_rates(file_dir, reaction_idx=None, tau=10.0, end_t=1.0, tag="fraction"):
    """
    plot reaction rates give reaction index list
    """

    colors, markers, _ = get_colors_markers_linestyles()

    _, rxn_idx_n = psri.parse_reaction_and_its_index(file_dir)
    rxn_idx_n["-1"] = "Temp"
    reaction_idx.append(-1)

    if reaction_idx is None:
        reaction_idx = [0]
    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    rxn_rates = np.loadtxt(os.path.join(file_dir, "output",
                                        "reaction_rate_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in reaction_idx:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=rxn_idx_n[str(s_idx)],
                           markevery=delta_n)
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            a_x_left.semilogy(time[0:end_point], rxn_rates[0:end_point, s_idx], marker=m_k,
                              color=colors[counter % (len(colors) - 1)],
                              label=rxn_idx_n[str(s_idx)],
                              markevery=delta_n)
            counter += 1
    leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})
    leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    a_x_left.set_xlabel("Time/sec")

    a_x_left.set_ylabel("R")
    a_x_right.set_ylabel("T/K")

    rxn_idx_str = "_".join(str(x) for x in reaction_idx)
    plt.title("reaction rates and Temp")

    fig.savefig(os.path.join(file_dir, "output",
                             "reaction_rate_" + rxn_idx_str + "_" + str(end_t) + ".jpg"), dpi=500)
    plt.close()


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    G_S = global_settings.get_setting(FILE_DIR)

    SPE_LIST = [3, 4, 5, 6, 7]
    # SPE_LIST, _, _ = trajectory.get_species_with_top_n_concentration(
    #     FILE_DIR, exclude=None, top_n=10,
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], end_t=G_S['end_t'],
    #     tag="M", atoms=["C", "O"])
    plot_concentrations(FILE_DIR, spe_idx=SPE_LIST,
                        tau=G_S['tau'], end_t=0.001, tag="M", exclude_names=None,
                        renormalization=False, semilogy=True, hasTemp=True)
