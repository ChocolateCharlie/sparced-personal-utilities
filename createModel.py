#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:04:32 2022
Sparced Model Creation
@author: charlie
"""

import sys
import os

import libsbml
import importlib
import amici
import numpy as np
import re
import pandas as pd
from antimony import *
from bin.modules.copyDir import copyDirectory

import amici.plotting
import matplotlib.pyplot as plt

# Input file name definitions
fileComps = "Compartments.txt"
fileSpecies = "Species.txt"
fileStoic = "StoicMat.txt"
fileRatelaws = "Ratelaws.txt"
# Output: lists all parameter names, rxn names and values
fileParamsOut = "ParamsAll.txt"

# Copy input files in current directory
current_dir = os.getcwd()+"/"
copyDirectory(current_dir+"input_files", current_dir)

"""
Antimony
"""

# Create Antimony model file
fileModelName = "SPARCED_tslap1.txt"
fileModel = open(fileModelName, "w")
fileModel.write("# Troubleshooting lapatinib concentration issues\nmodel SPARCED()\n\n")

# Initialize compartment and volume lists
compartments = []
volumes = []
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open(fileComps)])
for row in compartment_sheet[1:]:   # Skip header row
    compartments.append(row[0])
    volumes.append(row[1])

# Write Antimony Compartments
fileModel.write("# Compartments and Species:\n")
for idx in range(len(compartments)):
    compName = compartments[idx]
    fileModel.write("Compartment %s;" % (compName))
fileModel.write("\n\n")

# Write Antimony Species
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(fileSpecies)])
for idx,val in enumerate(species_sheet[1:]):    # Skip header row
    fileModel.write("Species %s in %s;\n" % (val[0], val[1]))
fileModel.write("\n")

# Write Antimony Reactions
fileModel.write("# Reactions:\n")
stoic_sheet = np.array([np.array(line.strip().split("\t")) for line in open(fileStoic)])
ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open(fileRatelaws)])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])

""" COPY/PASTE """
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
        fileModel.write("  %s: %s => %s; (%s)*%.6e;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, valcomp))

# Export parameters for each reaction, with corresponding order within the ratelaw and its value
params_all = pd.DataFrame({'value':paramvals,'rxn':paramrxns,'idx':paramidxs},index=paramnames)
params_all.to_csv(fileParamsOut,sep='\t',header=True, index=True)
""" END OF COPY/PASTE"""
fileModel.write("\n")

# Write Antimony compartment initial conditions
fileModel.write("# Compartment initializations:\n")
for idx in range(len(compartments)):
    fileModel.write("%s = %.6e;\n" % (compartments[idx], np.double(volumes[idx])))
    fileModel.write("%s has volume;\n" % (compartments[idx]))
fileModel.write("\n")

# Write Antimony species initial conditions
fileModel.write("# Species initializations:\n")
for idx, val in enumerate(species_sheet[1:]):
    fileModel.write("%s = %.6e;\n" % (val[0], np.double(val[2])))
fileModel.write("\n")

# Write Antimony parameter of reactions initial conditions
fileModel.write("# Parameter initializations:\n")
for idx, val in enumerate(paramnames):
    fileModel.write("%s = %.6e;\n" % (val, np.double(paramvals[idx])))
fileModel.write("\n")

# Write other Antimony declarations
constantVars = ['Cytoplasm', 'Extracellular', 'Nucleus', 'Mitochondrion']
fileModel.write("# Other declarations:\n")
fileModel.write("const ")
for constVar in constantVars[:-1]:
    fileModel.write("%s," % (constVar))
fileModel.write("%s;\n" % (constantVars[-1]))

# Write Antimony unit definition
fileModel.write("\n  # Unit definitions:")
fileModel.write("\n  unit time_unit = second;")
fileModel.write("\n  unit volume = litre;")
fileModel.write("\n  unit substance = 1e-9 mole;")
fileModel.write("\n  unit nM = 1e-9 mole / litre;")
fileModel.write("\n")

# Close the Antimony model file
fileModel.write("\nend")
fileModel.close()

"""
SBML
"""

# Create SBML model file
sbml_file = "SPARCED_tslap1.xml"
model_name = sbml_file[0:-4]
model_output_dir = model_name

# Load Antimony model
try:
    assert not loadFile(fileModelName) == -1
except:
    print("Failed to load antimony file")
    sys.exit(1)
else:
    print("Success loading antimony file")

# Convert Antimony model to SBML
try:
    assert not writeSBMLFile(sbml_file,"SPARCED") == 0
except:
    print("Failure converting antimony to SBML")
    sys.exit(1)
else:
    print("Success converting antimony to SBML")

# Import SBML file
sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

# Set species annotations using the last column of the Species input file
for idx, row in enumerate(species_sheet[1:]):
    Annot=""
    for col in range(4, (len(row))):
        aa = str(row[col].strip())
        if aa == "nan" or aa == "":
            break
        else:
            Annot = Annot + " " + row[col]
    sbml_model.getSpecies(row[0]).setAnnotation(Annot)

# Set compartment annotations using the last column of the Compartment input file
for row in compartment_sheet[1:]:
    sbml_model.getCompartment(row[0]).setAnnotation(row[2])

# Export the finalized SBML file
writer = libsbml.SBMLWriter()
writer.writeSBML(sbml_doc, sbml_file)

"""
Model Compilation
"""

# Import the SBML file
sys.path.insert(0, os.path.abspath(model_output_dir))
model_name = sbml_file[0:-4]
model_output_dir = model_name

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

sbml_importer = amici.SbmlImporter(sbml_file)

constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]

# Compile
sbml_importer.sbml2amici(model_name,
                         model_output_dir,
                         verbose=False)

print("The model is compiled!")
