#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: AUTOMATIC Post-Simulation Plot Sparced Data"""

import os
import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', default="plot",      help="name of the directory containing EXCLUSIVELY the data files")
    parser.add_argument('-t', '--time',      default="time.txt",  help="name of the timepoints file")
    parser.add_argument('-v', '--verbose',   action='store_true', help="display additional details during execution")
    return(parser.parse_args())

def read_time(f):
    with open(f, 'r') as f:
        t = f.read()
        time = t.split("\t")
        time.pop(-1) # Remove final ''
        time = [float(ele) / 3600.0 for ele in time]
    return(time)

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Load time
    f_time = args.time
    tt = read_time(f_time)

    # Load directory
    path = os.getcwd() + "/" + args.directory
    if args.verbose: print(path)
    files = os.listdir(path)
    if args.verbose: print(files)

    # Plot
    if args.verbose: print("SPARCED: Ready to plot")
    plot_loop = True
    while plot_loop:
        compound = input("Enter compound (see species list) : ")
        if compound != "STOP":
            fig = plt.figure(figsize=(8, 2), dpi=300, facecolor='w', edgecolor='k')
            # Each file in the directory will be loaded, which might cause a performance issue since this is performed on each loop
            for f in files:
                if os.path.isfile(os.path.join(path, f)):
                    # Load file
                    data = pd.read_csv(path + "/" + f, header=0, sep="\t")
                    del data[data.columns[0]] # Drop first column (index)
                    if args.verbose: print("SPARCED: Successfully loaded " + f)
                    # Plot data
                    plt.plot(tt, data[compound], color='b', label=f[:-4])
            # Figure settings
            plt.xlim([-0.5, tt[-1]+1])
            plt.legend()
            plt.ylabel('Concentration (nM)')
            plt.xlabel('Time (hours)')
            plt.title(compound)
            plt.xticks(np.arange(0, tt[-1]+1, step=1))
            plt.savefig(compound, format='jpeg')
            plt.show()
        else:
            plot_loop = False
    if args.verbose: print("SPARCED: End of plot loop")
