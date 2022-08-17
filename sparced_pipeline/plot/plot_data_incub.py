#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Post-Simulation Plot Sparced Data with Incubation"""

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control',      default="control.txt",          help="name of the control file")
    parser.add_argument('-e', '--experiment',   default="experiment.txt",       help="name of the experiment file")
    parser.add_argument('-i', '--incubcontrol', default="control_incub.txt",    help="name of the control incubation file")
    parser.add_argument('-j', '--incubexp',     default="experiment_incub.txt", help="name of the experimnet incubation file")
    parser.add_argument('-l', '--clabel',       default="Control",              help="label for control")
    parser.add_argument('-n', '--elabel',       default="Experiment",           help="label for experiment")
    parser.add_argument('-s', '--incubtime',    default="time_inucb.txt",       help="name of the incubation timepoints file")
    parser.add_argument('-t', '--time',         default="time.txt",             help="name of the timepoints file")
    parser.add_argument('-v', '--verbose',      action='store_true',            help="display additional details during execution")
    return(parser.parse_args())

def read_file(f):
    file = pd.read_csv(f, header=0, sep="\t")
    del file[file.columns[0]] # Drop first column (index)
    return(file)

def read_time(f):
    with open(f, 'r') as f:
        t = f.read()
        time = t.split("\t")
        time.pop(-1) # Remove final ''
        time = [float(ele) / 3600.0 for ele in time]
    return(time)

def set_filenames(args):
    return(args.control, args.incubcontrol, args.experiment, args.incubexp, args.time, args.incubtime)

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # LOAD DATA
    f_cont, f_cont_incub, f_exp, f_exp_incub, f_time, f_time_incub = set_filenames(args)
    # Control
    control = read_file(f_cont_incub)
    control_simulation = read_file(f_cont)
    control = control.append(control_simulation, ignore_index=True)
    # Experiment
    experiment = read_file(f_exp_incub)
    experiment_simulation = read_file(f_exp)
    experiment = experiment.append(experiment_simulation, ignore_index=True)
    # Time
    t_incubation = read_time(f_time_incub)
    t_simul = read_time(f_time)
    for t in range(0,len(t_simul)):
        t_simul[t] += t_incubation[-1]
    tt = t_incubation + t_simul

    # Plot
    if args.verbose: print("SPARCED: Ready to plot")
    plot_loop = True
    while plot_loop:
        compound = input("Enter compound (see species list) : ")
        if compound != "STOP":
            fig = plt.figure(figsize=(8, 2), dpi=300, facecolor='w', edgecolor='k')
            plt.plot(tt, control[compound], 'b', linewidth=2, markersize=1, label=args.clabel)
            plt.plot(tt, experiment[compound], 'r', linewidth=2, markersize=1, label=args.elabel)
            plt.axvline(t_incubation[-1], color='g')
            plt.xlim([-0.5, tt[-1]+1])
            plt.legend()
            plt.ylabel('Concentration (nM)')
            plt.xlabel('Time (hours)')
            plt.title(compound)
            plt.xticks(np.arange(0, tt[-1]+1, step=0.5))
            plt.savefig(compound, format='jpeg')
            plt.show()
        else:
            plot_loop = False
    if args.verbose: print("SPARCED: End of plot loop")
