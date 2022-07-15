#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Creation"""

import os
from antimony import *
import argparse
import numpy as np
from bin.copydir import copy_directory

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--antimony',     default="SPARCED",          help="name of the antimony model's name")
    parser.add_argument('-c', '--compartments', default="Compartments.txt", help="name of the compartments' file")
    parser.add_argument('-i', '--inputdir',     default="/input_files",     help="relative path to input files directory")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt",     help="name of the stoichiometric matrix' file")
    parser.add_argument('-o', '--outputparams', default="ParamsAll.txt",    help="name of the output parameters' file")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt",     help="name of the rate laws' file")
    parser.add_argument('-s', '--species',      default="Species.txt",      help="name of the species' file")
    return(parser.parse_args())

def set_io_filenames(args):
    return(args.compartments, args.stoichmatrix, args.outputparams, args.ratelaws, args.species)

def antimony_init(f):
    comp = []
    vol = []
    sheet = np.array([np.array(line.strip().split("\t")) for line in open(f)])
    for row in sheet[1:]: # Skip header row
        comp.append(row[0])
        vol.append(row[1])
    return((comp, vol))

def antimony_write_compartments(f, comp):
    f.write("# Compartments and Species:\n")
    for i in range(len(comp)):
        f.write("Compartment {comp_name}; ".format(comp_name=comp[i]))
    f.write("\n\n")
    

if __name__ == '__main__':
    args = parse_args()

    # Initialize filenames
    f_comp, f_stoi, f_outp, f_rate, f_spec = set_io_filenames(args)
    antimony_model_name = args.antimony

    # Copy intput files in working directory
    current_dir = os.getcwd()
    copy_directory(current_dir + args.inputdir, current_dir)

    # ANTIMONY
    with open(antimony_model_name + ".txt", "w") as antimony_model:
        antimony_model.write("# PanCancer Model by Birtwistle Lab\n")
        antimony_model.write("model {antimony}()\n\n".format(antimony=args.antimony))
        # Initialize compartment and volume lists
        compartments, volumes = antimony_init(f_comp)
        # Write antimony compartments
        antimony_write_compartments(antimony_model, compartments)
