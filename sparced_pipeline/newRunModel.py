#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Simulation"""

import os
import sys

import argparse
import importlib
import libsbml

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--sbml',     default="SPARCED", help="name of the SBML model")
    parser.add_argument('-c', '--compound', default=None,      help="compound's name")
    parser.add_argument('-d', '--dose',     default=0.0,       help="compound's concentration in nM")
    parser.add_argument('-e', '--egf',      default=1.0,       help="EGF concentration in nM")
    parser.add_argument('-n', '--name',     default="sparced", help="name of the simulation")
    parser.add_argument('-p', '--pop',      default=3,         help="number of cells in the population")
    parser.add_argument('-t', '--time',     default=72.0,      help="duration of the simulation in virtual time")
    parser.add_argument('-v', '--verbose',  default=True,      help="display additional details during execution")
    return(parser.parse_args())

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()
# Import model
sbml_model_name = model_output_dir = args.sbml
sbml_file_name = sbml_model_name + ".xml"
sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(sbml_model_name)
if args.verbose: print(model_module)

