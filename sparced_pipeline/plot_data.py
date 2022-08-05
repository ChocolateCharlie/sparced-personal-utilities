#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Post-Simulation Plot Sparced Data"""

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control',    default="control.txt",      help="name of the control file")
    parser.add_argument('-e', '--experiment', default="experiment.txt",   help="name of the experiment file")
    parser.add_argument('-s', '--species',    default="species_list.txt", help="name of the list of the species")
    parser.add_argument('-t', '--time',       default="time.txt",         help="name of the timepoints file")
    parser.add_argument('-v', '--verbose',    action='store_true',        help="display additional details during execution")
    return(parser.parse_args())

def read_species(f):
    with open(f, 'r') as f:
        s = f.read()
        species = s.split("\t")
    return(species)

def set_filenames(args):
    return(args.control, args.experiment, args.species, args.time)

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Load data
    f_cont, f_exp, f_spec, f_time = set_filenames(args)
    control = pd.read_csv(f_cont, header=0, sep="\t")
    del control[control.columns[0]] # Drop first column (index)
    experiment = pd.read_csv(f_exp, header=0, sep="\t")
    del experiment[experiment.columns[0]] # Drop first column (index)
    species = read_species(f_spec)
    time = np.array(line.strip().split("\t") for line in open(f_time))
    tt = time / 3600.0

    # Plot
    if args.verbose: print("SPARCED: Ready to plot")
    plot_loop = True
    while plot_loop:
        compound = input("Enter compound (see species list) : ")
        if compound != "STOP":
            fig = plt.figure(figsize=(8, 2), dpi=300, facecolor='w', edgecolor='k')
            plt.plot(tt, control[:,species.index(compound)], 'b', linewidth=2, markersize=1, label='control')
            plt.plot(tt, experiment[:,species.index(compound)], 'r', linewidth=2, markersize=1, label='experiment')
            plt.xlim([-0.5, tt[-1]+1])
            plt.legend()
            plt.ylabel('Concentration (nM)')
            plt.xlabel('Time (hours)')
            plt.title('Name')
            plt.xticks(np.arange(0, tt[-1]+1, step=1))
            plt.savefig('name', format='jpeg')
            plt.show()
        else:
            plot_loop = False
    if args.verbose: print("SPARCED: End of plot loop")
