#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Simulation

This is the fast run model code adjusted to perform incubation after adding
the compound into the cell.
After the compound is added, the cell is incubated during the specified
duration. Please note that the simulation will start only after this step,
including the addition of growth factors.
"""

import os
import sys

import argparse
import importlib
import libsbml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bin.run_sparced_fast import run_sparced_fast

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--sbml',          default="SPARCED",          help="name of the SBML model")
    parser.add_argument('-c', '--compound',      default=None,               help="compound's name")
    parser.add_argument('-d', '--dose',          default=0.0,                help="compound's concentration in nM")
    parser.add_argument('-D', '--deterministic', action='store_const',       help="flag D", const=1, default=0)
    parser.add_argument('-e', '--egf',           default=1.0,                help="EGF concentration in nM")
    parser.add_argument('-i', '--ins',           default=17.21,              help="EGF concentration in nM")
		parser.add_argument('-k', '--incubation',    default=24.0,               help="duration of the incubation")
    parser.add_argument('-l', '--list',          default="species_list.txt", help="name of the list of species file")
    parser.add_argument('-n', '--name',          default="GrowthStim",       help="name of the simulation")
    parser.add_argument('-p', '--pop',           default=1,                  help="number of cells in the population")
    parser.add_argument('-s', '--species',       default="Species.txt",      help="name of the species' file")
    parser.add_argument('-t', '--time',          default=1.0,                help="duration of the simulation in virtual time")
    parser.add_argument('-v', '--verbose',       action='store_true',        help="display additional details during execution")
    parser.add_argument('-x', '--exchange',      default=30,                 help="information exchange between modules timeframe")
    return(parser.parse_args())

def save_output(model, file_prefix, cell_number, xoutS_all, xoutG_all, tout_all):
    # xoutS
    columnsS = [ele for ele in model.getStateIds()]
    condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS)
    condsSDF.to_csv(file_prefix+'_S_'+str(cell_number)+'.txt',sep="\t")  
    condsSDF = None
    #xoutG
    columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
    columnsG = columnsG[1:] # Skip header
    resa = [sub.replace('m_', 'ag_') for sub in columnsG]
    resi = [sub.replace('m_', 'ig_') for sub in columnsG]
    columnsG2 = np.concatenate((resa, resi), axis=None)
    condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
    condsGDF.to_csv(file_prefix+'_G_'+str(cell_number)+'.txt',sep="\t") 
    condsGDF = None
    #tout
    np.savetxt(file_prefix+'_T_'+str(cell_number)+'.txt', tout_all, newline="\t", fmt="%s")

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Import model
    sbml_model_name = model_output_dir = args.sbml
    sys.path.insert(0, os.path.abspath(model_output_dir))
    model_module = importlib.import_module(sbml_model_name)
    if args.verbose: print(model_module)
    # Set model
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, args.exchange, 2))
    # Save list of species
    species_all = list(model.getStateIds())
    np.savetxt(args.list, species_all, newline="\t", fmt="%s")

    # Run simulations
    cell_number = 0
    while cell_number < int(args.pop):
        # INITIALIZATION
        # Set initial conditions
        species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(args.species, encoding='latin-1')])
        species_initializations = []
        for row in species_sheet[1:]:
            species_initializations = np.append(species_initializations, float(row[2]))
            species_initializations = np.array(species_initializations)
            species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
				# INCUBATION
				# Add compound
				if args.compound is no None: species_initializations[species_all.index(args.compound)] = args.dose
				model.setInitialStates(species_initializations)
				if args.verbose: print("SPARCED: Now ready to start incubation")
				xoutS_incub, xoutG_incub, tout_incub = run_sparced_fact(args.deterministic, float(args.incubation), species_initializations, sbml_model_name + ".xml", model)
				save_output(model, args.name + "_incub", cell_number, xoutS_incub, xoutG_incub, tout_incub)
				if args.verbose: print("SPARCED: Incubation is now over")
				# SIMULATION
				for idx in range(len(species_initializations)):
					species_initializations[idx] = xoutS_incub[-1:,idx]
				model.setInitialStates(species_initializations)
			  # Input ligand concentrations (in order): EGF, Her, HGF, PDGF, FGF, IGF, INS
				STIMligs = [float(args.egf), 0.0, 0.0, 0.0, 0.0, 0.0, float(args.ins)] # in nM, in extracellular volume
        species_initializations[155:162] = STIMligs
        if args.verbose: print("SPARCED: Now ready to run a simulation")
        xoutS_all, xoutG_all, tout_all = run_sparced_fast(args.deterministic, float(args.time), species_initializations, sbml_model_name + ".xml", model)
        # SAVE OUTPUT
        save_output(model, args.name, cell_number, xoutS_all, xoutG_all, tout_all)
        if args.verbose: print("SPARCED: Simulation is now over")
        cell_number += 1
