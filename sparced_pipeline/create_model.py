#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Creation"""

import os
import sys

import amici
from antimony import *
import argparse
import libsbml

from bin.antimony_utils import *
from bin.copydir import copy_directory
from bin.sbml_utils import *

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--antimony',     default="SPARCED",          help="name of the antimony model")
    parser.add_argument('-b', '--sbml',         default="SPARCED",          help="name of the SBML model")
    parser.add_argument('-c', '--compartments', default="Compartments.txt", help="name of the compartments' file")
    parser.add_argument('-i', '--inputdir',     default="/input_files",     help="relative path to input files directory")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt",     help="name of the stoichiometric matrix' file")
    parser.add_argument('-o', '--outputparams', default="ParamsAll.txt",    help="name of the output parameters' file")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt",     help="name of the rate laws' file")
    parser.add_argument('-s', '--species',      default="Species.txt",      help="name of the species' file")
    parser.add_argument('-v', '--verbose',      default=True,               help="display additional details during execution")
    return(parser.parse_args())

def set_io_filenames(args):
    return(args.compartments, args.stoichmatrix, args.outputparams, args.ratelaws, args.species)

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Initialize filenames
    f_comp, f_stoi, f_outp, f_rate, f_spec = set_io_filenames(args)
    antimony_model_name = args.antimony
    antimony_file_name = antimony_model_name + ".txt"
    sbml_model_name = model_output_dir = args.sbml
    sbml_file_name = sbml_model_name + ".xml"

    # Copy intput files in working directory
    current_dir = os.getcwd()
    copy_directory(current_dir + args.inputdir, current_dir)

    # Antimony
    with open(antimony_file_name, "w") as antimony_model:
        # Write antimony file's header
        antimony_model.write("# PanCancer Model by Birtwistle Lab\n")
        antimony_model.write("model {antimony}()\n\n".format(antimony=antimony_model_name))
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

    # Load Antimony model
    try:
        assert not loadFile(antimony_file_name) == -1
    except:
        print("SPARCED: Failed to load Antimony file")
        sys.exit(0)
    else:
        if args.verbose: print("SPARCED: Success loading Antimony file")

    # Convert Antimony model into SBML
    try:
        assert not writeSBMLFile(sbml_file_name, antimony_model_name) == 0
    except:
        print("SPARCED: Failed to convert Antimony file to SBML")
        sys.exit(0)
    else:
        if args.verbose: print("SPARCED: Success converting Antimony file to SBML")

    # SBML: Annotation
    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    # Set species annotations
    write_species_annotations(sbml_model, species)
    # Set compartments annotations
    # write_compartments_annotations(sbml_model, compartments)
    # Export the annotated SBML file
    writer = libsbml.SBMLWriter()
    writer.writeSBML(sbml_doc, sbml_file_name)

    # SBML: Compilation
    # Import annotated SBML file
    sys.path.insert(0, os.path.abspath(model_output_dir))
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    sbml_importer = amici.SbmlImporter(sbml_file_name)
    const_params = [params.getId() for params in sbml_model.getListOfParameters()]
    # Compile
    sbml_importer.sbml2amici(sbml_model_name, model_output_dir, verbose=args.verbose)
    if args.verbose: print("SPARCED: Success compiling the model")