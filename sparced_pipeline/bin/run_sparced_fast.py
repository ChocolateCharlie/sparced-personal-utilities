#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Arnab's fast simulation run"""

import os
import libsbml
import importlib
import amici
import numpy as np
import pandas as pd

from bin.SGEmodule import SGEmodule
from bin.run_prep import run_prep

def run_sparced_fast(flagD, th, spdata, sbml_file, omics_input='OmicsData.txt',genereg_input='GeneReg.txt'):
    wd = str(os.getcwd())

    ts = 30 # time-step to update mRNA numbers
    NSteps = int(th*3600/ts)
    tout_all = np.arange(0,th*3600+1,30)   
    
    # Read-in the model SBML to get compartmental volumes (used to convert nM to mpc and vice versa)
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file)
    model = sbml_doc.getModel()
    Vc = model.getCompartment(0).getVolume() # Should be the index for Cytoplasm
    Vn = model.getCompartment(2).getVolume() # Should be the index for Nuclues
    mpc2nM_Vc = (1E9/(Vc*6.023E+23))
    splist = list(model.getStateIds())
    
    PARPind = [ind for ind,ele in enumerate(splist) if ele in {'PARP'}] # find the index for PARP
    cPARPind = [ind for ind,ele in enumerate(splist) if ele in {'cPARP'}]
    
    genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs, mrna_idx = run_prep(flagD,Vn,model,wd,omics_input,genereg_input)
    
    if len(spdata)==0:
        spdata0 = pd.read_csv(os.path.join(wd,'input_files','Species.txt'),header=0,index_col=0,sep="\t")
        spdata = np.float(spdata0.values[:,1])
    xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
    xoutS_all[0,:] = spdata # 24hr time point     
    
    if len(genedata)==0:
        print("genedata0")
        exit()
        # genedata = genedata0
    xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
    xoutG_all[0,:] = genedata
    
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    
    flagA = 0
    
    n_sp = len(splist)
    
    for qq in range(NSteps):
        genedata,xmN,AllGenesVec = SGEmodule(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs,spIDs,mrna_idx)
        xoutS_all[qq,mrna_idx] = np.dot(xmN,mpc2nM_Vc)
        model.setInitialStates(xoutS_all[qq,:])
        rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
        xoutS_all[qq+1,:] = rdata._swigptr.x[-n_sp:]
        xoutG_all[qq+1,:] = genedata
        rdata = None
        
        if xoutS_all[qq+1,PARPind] < xoutS_all[-1,cPARPind]: 
            print('Apoptosis happened')
            flagA = 1
            break
        

    xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
    xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
    tout_all = tout_all[0:len(xoutS_all)]
    
    return xoutS_all, xoutG_all, tout_all, flagA