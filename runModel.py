#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:37:30 2022
Sparced run model utility
@author: charlie
"""

import libsbml
import importlib
import amici
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import argparse
import sys
import os
from bin.modules.RunSPARCED import RunSPARCED

"""
Import model
"""
#%%
# Model import
sbml_file = "SPARCEDtslap1.xml"
model_name = sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, "/home/charlie/Documents/troubleshoot-lapatinib/SPARCEDtslap1")  # TODO: FIX PATH
ts = 30 # The time frame at which stochastic gene module and deterministic SBML module exhange/update information
cellNumber = 0 # used for filename to save
model_module = importlib.import_module(model_name)
print(model_module)
model = model_module.getModel()
solver = model.getSolver() # Create solver instance
solver.setMaxSteps = 1e10

#%%

"""
Custom inputs
"""
#%%
# Inputs
flagD = 1
th = 1
# Input ligand concentrations (in order): EGF, Her, HGF, PDGF, FGF, IGF, INS
# TIMligs = [100.0,100.0,100.0,100.0,100.0,100.0,1721.0] # in nM, in extracellular volume
# File name prior for saving the simulation trajectory and gene states
nmoutfile = 'GrowthStim_'

# Species initial conditions
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
# species_initializations[155:162] = STIMligs
#%%

"""
Simulation
"""
"""
# Input ligand concentrations (in order): EGF, Her, HGF, PDGF, FGF, IGF, INS
STIMligs = [100.0,0.0,0.0,0.0,0.0,0.0,0.0] # in nM, in extracellular volume
species_initializations[155:162] = STIMligs
species_initializations[167] = 3.8
"""
th = 1
# Control
"""
species_initializations[155] = 100 # EGF
species_initializations[773] = 0 # Lapatinib
species_initializations[162] = 50 # HER1
species_initializations[164] = 3.8 # HER2
species_initializations[179] = 0.0144 # E2E2
"""
species_initializations[0] = 1000 # Lapatinib
model.setInitialStates(species_initializations) 
model.setTimepoints(np.linspace(0,th,3600*th+1)) # set timepoints
rdata = amici.runAmiciSimulation(model, solver)
"""
# xoutS_all_control, xoutG_all_control, tout_all_control = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)
# saves output
columnsS = [ele for ele in model.getStateIds()]
columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
columnsG = columnsG[1:]
resa = [sub.replace('m_', 'ag_') for sub in columnsG]
resi = [sub.replace('m_', 'ig_') for sub in columnsG]
columnsG2 = np.concatenate((resa, resi), axis=None)
condsSDF = pd.DataFrame(data=xoutS_all_control,columns=columnsS)
condsSDF.to_csv(nmoutfile+'S_'+str(cellNumber)+'.txt',sep="\t")  
condsSDF = None
condsGDF = pd.DataFrame(data=xoutG_all_control,columns=columnsG2)
condsGDF.to_csv(nmoutfile+'G_'+str(cellNumber)+'.txt',sep="\t") 
condsGDF = None
cellNumber+=1
fig = plt.figure(figsize=(8, 2), dpi=300, facecolor='w', edgecolor='k')
tt = tout_all_control/3600.0
"""
# Lapatinib
"""
species_initializations[155] = 100 # EGF
species_initializations[773] = 1000 # Lapatinib
species_initializations[162] = 50 # HER1
species_initializations[164] = 3.8 # HER2
species_initializations[179] = 0.0144 # E2E2
xoutS_all_lapatinib, xoutG_all_lapatinib, tout_all_lapatinib = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)
# saves output
columnsS = [ele for ele in model.getStateIds()]
columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
columnsG = columnsG[1:]
resa = [sub.replace('m_', 'ag_') for sub in columnsG]
resi = [sub.replace('m_', 'ig_') for sub in columnsG]
columnsG2 = np.concatenate((resa, resi), axis=None)
condsSDF = pd.DataFrame(data=xoutS_all_lapatinib,columns=columnsS)
condsSDF.to_csv(nmoutfile+'S_'+str(cellNumber)+'.txt',sep="\t")  
condsSDF = None
condsGDF = pd.DataFrame(data=xoutG_all_lapatinib,columns=columnsG2)
condsGDF.to_csv(nmoutfile+'G_'+str(cellNumber)+'.txt',sep="\t") 
condsGDF = None
cellNumber+=1
"""
#%%
"""
# Cheatsheet
species_all = list(model.getStateIds())
print(species_all.index'Specie')
print(species_initializations[index])
"""
#%%
# Plots
fig = plt.figure(figsize=(8, 2), dpi=300, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(1,1,1)
#plt.plot(tt, xoutS_all_normal[:,717], 'blue', linewidth=2, markersize=12, label='Normal')
#plt.plot(tt, xoutS_all_normal_resp[:,717], 'magenta', linewidth=2, markersize=12, label='HER2-')
plt.plot(tt, xoutS_all_control[:,777], 'blue', linewidth=2, markersize=12, label='Control')
plt.plot(tt, xoutS_all_lapatinib[:,777], 'red', linewidth=2, markersize=12, label='Lapatinib')
plt.xlim([-1,th+1])
plt.legend()
ax1.set_ylabel('Concentration of lap_E4 (nM)')
ax1.set_xlabel('Time (hr)')
plt.xticks(np.arange(0, th+1, step=4))  # Set label locations.
