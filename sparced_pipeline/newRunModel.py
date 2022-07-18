#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Simulation"""

import os
import sys

import argparse
import importlib
import libsbml
import matplotlib.pyplot as plt
import numpy as np
import panda as pd

from bin.modules.RunSPARCED import RunSPARCED

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--sbml',          default="SPARCED",     help="name of the SBML model")
    parser.add_argument('-c', '--compound',      default=None,          help="compound's name")
    parser.add_argument('-d', '--dose',          default=0.0,           help="compound's concentration in nM")
    parser.add_argument('-D', '--deterministic', default=1,             help="flag D")
    parser.add_argument('-e', '--egf',           default=1.0,           help="EGF concentration in nM")
    parser.add_argument('-n', '--name',          default="sparced",     help="name of the simulation")
    parser.add_argument('-p', '--pop',           default=1,             help="number of cells in the population")
    parser.add_argument('-s', '--species',       default="Species.txt", help="name of the species' file")
    parser.add_argument('-t', '--time',          default=1.0,           help="duration of the simulation in virtual time")
    parser.add_argument('-v', '--verbose',       default=True,          help="display additional details during execution")
    parser.add_argument('-x', '--exchange',      default=30,            help="information exchange between modules timeframe")
    return(parser.parse_args())

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Import model
    sbml_model_name = model_output_dir = args.sbml
    sys.path.insert(0, os.path.abspath(model_output_dir))
    model_module = importlib.import_module(sbml_model_name)
    if args.verbose: print(model_module)

    # Set parameters
    cell_number = 0
    nmoutfile = 'GrowthStim_' # File prior for saving the simulation trajectory and gene states

    # Create solver instance
    model = model_module.getModel()
    solver = model.getSolver()
    solver.setMaxSteps = 1e10    
    model.setTimepoints(np.linspace(0, args.exchange, 2))
    # model.setTimepoints(np.linspace(0, th*3600, int(3600 * th + 1)))

    # Run simulations
    while cell_number < args.pop:
        # INITIALIZATION
        # Set initial conditions
        species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(args.species)])
        species_initializations = []
        for row in species_sheet[1:]: # Skip header
            species_initializations.append(float(row[2]))
        species_initializations = np.array(species_initializations)
        species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
        # Input ligand concentrations (in order): EGF, Her, HGF, PDGF, FGF, IGF, INS
        STIMligs = [args.egf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # in nM, in extracellular volume
        species_initializations[155:162] = STIMligs
        # Compound
        species_all = list(model.getStateIds())
        if args.compound is not None: species_initializations[species_all.index(args.compound)] = args.dose
        # model.setInitialStates(species_initializations)
        # SIMULATION
        xoutS_all, xoutG_all, tout_all = RunSPARCED(args.deterministic, args.time, species_initializations,[], sbml_model_name + ".xml", model)
        # SAVE OUTPUT
        columnsS = [ele for ele in model.getStateIds()]
        columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
        columnsG = columnsG[1:]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_csv(nmoutfile+'S_'+str(cell_number)+'.txt',sep="\t")  
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_csv(nmoutfile+'G_'+str(cell_number)+'.txt',sep="\t") 
        condsGDF = None
        cell_number += 1
