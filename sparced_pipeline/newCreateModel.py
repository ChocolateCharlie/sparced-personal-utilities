#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Creation"""

import os
from antimony import *
import argparse
import numpy as np
import pandas as pd
import re
from bin.copydir import copy_directory

# ---------------------------------------------------------------------------------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def antimony_init(f_cv, f_s):
    # Compartments and Volumes
    comp = []
    vol = []
    sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_cv)])
    for row in sheet[1:]: # Skip header row
        comp.append(row[0])
        vol.append(row[1])
    # Species
    spec = np.array([np.array(line.strip().split("\t")) for line in open(f_s)], dtype="object")
    return((comp, vol, spec))

def antimony_terminal(f):
    # Constant variables
    f.write("# Other declarations:\nconst ")
    constants = ['Cytoplasm', 'Extracellular', 'Nucleus', 'Mitochondrion']
    for c in constants[:-1]:
        f.write("{name} ".format(name=c))
    f.write("{last_name}\n\n".format(last_name=constants[-1]))
    # Unit definitions
    f.write("# Unit definitions:\n")
    f.write("  unit time_unit = second;\n")
    f.write("  unit volume = litre;\n")
    f.write("  unit substance = 1e-9 mole;\n")
    f.write("  unit nM = 1e-9 mole / litre;\n\n")

def antimony_write_compartments(f, comp):
    f.write("# Compartments and Species:\n")
    for i in range(len(comp)):
        f.write("Compartment {comp_name}; ".format(comp_name=comp[i]))
    f.write("\n\n")

def antimony_write_init_compartments(f, comp, vol):
    f.write("# Compartments initialization:\n")
    for i in range(len(comp)):
        f.write("{name} = {volume}.6e;\n{name} has volume;\n".format(name=comp[i], volume=np.double(vol[i])))
    f.write("\n")

def antimony_write_init_reactions(f, p_names, p_vals):
    f.write("# Reactions' parameters initializations:\n ")
    for i, val in enumerate(p_names):
        f.write("{name} = {value}.6e;\n".format(name=val, value=np.double(p_vals[i])))
    f.write("\n")

def antimony_write_init_species(f, spec):
    f.write("# Species initialization:\n")
    for i, val in enumerate(species[1:]): # Skip header
        f.write("{name} = {concentration}.6e;\n".format(name=val[0], concentration=np.double(val[2])))
    f.write("\n")

def antimony_write_reactions(f, f_rl, f_sm, f_outp):
    f.write("# Reactions:\n")
    stoic_sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_sm)], dtype="object")
    ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_rl)], dtype="object")
    ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]], dtype="object")
    # ========== COPY/PASTE ==========
    #gets first column minus blank space at the beginning, adds to stoic data list
    stoic_columnnames = stoic_sheet[0]
    stoic_rownames = [line[0] for line in stoic_sheet[1:]]
    stoic_data = np.array([line[1:] for line in stoic_sheet[1:]])
    # builds the important ratelaw+stoic lines into the txt file 
    paramnames = []
    paramvals = []
    paramrxns = []
    paramidxs = []
    for rowNum, ratelaw in enumerate(ratelaw_data):
        reactants = []
        products = []
        formula="k"+str(rowNum+1)+"*"
    
        for i, stoic_rowname in enumerate(stoic_rownames):
            stoic_value = int(stoic_data[i][rowNum])
            if stoic_value < 0:
                for j in range(0,stoic_value*-1):
                    reactants.append(stoic_rowname)
                    formula=formula+stoic_rowname+"*"
            elif stoic_value > 0:
                for j in range(0,stoic_value):
                    products.append(stoic_rowname)
    
        if "k" not in ratelaw[1]:
            # the mass-action formula
            formula=formula[:-1]
            #the parameter
            paramnames.append("k"+str(rowNum+1))
            paramvals.append(np.double(ratelaw[1]))
            paramrxns.append(ratelaw_sheet[rowNum+1][0])
            paramidxs.append(int(0))
        else:
            # specific formula (non-mass-action)
            formula = ratelaw[1]
            j = 1
            params = np.genfromtxt(ratelaw[2:], float) # parameters
            params = params[~np.isnan(params)]
            if len(params) == 1:
                paramnames.append("k"+str(rowNum+1)+"_"+str(j))
                paramvals.append(float(ratelaw[j+1]))
                paramrxns.append(ratelaw_sheet[rowNum+1][0])
                paramidxs.append(int(0))
                pattern = 'k\D*\d*'
                compiled = re.compile(pattern)
                matches = compiled.finditer(formula)
                for ematch in matches:
                    formula = formula.replace(ematch.group(),paramnames[-1])
            else:
                for q,p in enumerate(params):
                    paramnames.append("k"+str(rowNum+1)+"_"+str(j))
                    paramvals.append(float(ratelaw[j+1]))
                    paramrxns.append(ratelaw_sheet[rowNum+1][0])
                    paramidxs.append(q)
                    pattern1 = 'k(\D*)\d*'+'_'+str(j)
                    compiled1 = re.compile(pattern1)
                    matches1 = compiled1.finditer(formula)
                    for ematch in matches1:
                        formula = formula.replace(ematch.group(),paramnames[-1])
                    j +=1
        if ratelaw[0] == 'Cytoplasm':
            valcomp = 5.25e-12
        elif ratelaw[0] == 'Extracellular':
            valcomp = 5.00e-5
        elif ratelaw[0] == 'Nucleus':
            valcomp = 1.75e-12
        elif ratelaw[0] == 'Mitochondrion':
            valcomp = 3.675e-13
        #don't include reactions without products or reactants
        if products == [] and reactants == []:
            pass
        else:
            f.write("  %s: %s => %s; (%s)*%.6e;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, valcomp))
    
    # Export parameters for each reaction, with corresponding order within the ratelaw and its value
    params_all = pd.DataFrame({'value':paramvals,'rxn':paramrxns,'idx':paramidxs},index=paramnames)
    params_all.to_csv(f_outp,sep='\t',header=True, index=True)
    # ========== END OF COPY/PASTE ==========
    f.write("\n")
    return((paramnames, paramvals))

def antimony_write_species(f, spec):
    for i, val in enumerate(spec[1:]): # Skip header row
        f.write("Species {name} in {compartment}\n".format(name=val[0], compartment=val[1]))
    f.write("\n")
    
# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Initialize filenames
    f_comp, f_stoi, f_outp, f_rate, f_spec = set_io_filenames(args)
    antimony_model_name = args.antimony

    # Copy intput files in working directory
    current_dir = os.getcwd()
    copy_directory(current_dir + args.inputdir, current_dir)

    # Antimony
    with open(antimony_model_name + ".txt", "w") as antimony_model:
        # Write antimony file's header
        antimony_model.write("# PanCancer Model by Birtwistle Lab\n")
        antimony_model.write("model {antimony}()\n\n".format(antimony=args.antimony))
        # Initialize compartment and volume lists
        compartments, volumes, species = antimony_init(f_comp, f_spec)
        # Write antimony compartments, species and reactions
        antimony_write_compartments(antimony_model, compartments)
        antimony_write_species(antimony_model, species)
        param_names, param_vals = antimony_write_reactions(antimony_model, f_rate, f_stoi, f_outp)
        # Write antimony initial conditions (compartments, volumes, species, reactions' parameters)
        antimony_write_init_compartments(antimony_model, compartments, volumes)
        antimony_write_init_species(antimony_model, species)
        antimony_write_init_reactions(antimony_model, param_names, param_vals)
        # Write other declarations and unit definitions
        antimony_terminal(antimony_model)
        antimony_model.write("\nend") # End of antimony file


