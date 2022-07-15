#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Creation"""

import os
import argparse
from bin.copydir import copy_directory

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument()
    parser.add_argument('-a', '--antimony',     default="SPARCED.txt",      help="name of the antimony model's name")
    parser.add_argument('-c', '--compartments', default="Compartments.txt", help="name of the compartments' file")
    parser.add_argument('-i', '--inputdir',     default="/input_files",     help="relative path to input files directory")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt",     help="name of the stoichiometric matrix' file")
    parser.add_argument('-o', '--outputparams', default="ParamsAll.txt",    help="name of the output parameters' file")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt",     help="name of the rate laws' file")
    parser.add_argument('-s', '--species',      default="Species.txt",      help="name of the species' file")
    return(parser.parse_args())

def set_io_filenames(args):
    return(args.compartments, args.stoichmatrix, args.outputparams, args.ratelaws, args.species)

if __name__ == '__main__':
    args = parse_args()

    # Initialize filenames
    f_comp, f_stoi, f_outp, f_rate, f_spec = set_io_filenames(args)
    antimony_model_name = args.antimony

    # Copy intput files in working directory
    current_dir = os.getcwd()
    copy_directory(current_dir + args.inputdir, current_dir)

    # Antimony
    antimony_model = open(antimony_model_name, "w")
# To be continued...
    antimony_model.write("model SPARCED()\n\n")