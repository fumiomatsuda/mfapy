#-------------------------------------------------------------------------------
# Name:        mfapy example 0 toymodel Check MDV calc
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import numpy as numpy
import mfapy
from matplotlib.pyplot import figure, show
import os, sys, time

if __name__ == '__main__':
    #
    # Load metabolic model from txt file to four dictionary
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_0_toymodel_model.txt", format = "text")
    #
    # Construct MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    #
    # Set constraints manusally
    #
    model.set_constrain('reaction', 'v1', "fixed", 100.0, 1)
    #
    # self.update() is required after manusal setting constraints
    #
    model.update()
    #
    # Load metabolic state from file
    #
    flux = model.load_states("Example_0_toymodel_status.csv", format = 'csv')
    #
    # Generate instances of CarbonSource class from model
    #
    cs = model.generate_carbon_source_templete()
    #
    # Set isotope labelling of carbon sources by ratio of all isotopomer (#00, #01, #10, #11)
    #
    cs.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.25, 0.25])#a mixture of 25% [2-13C]AcCoA and 25%[1,2-13C]AcCoA
    #
    # Set isotope labelling of carbon sources by specific isotopomers
    #
    cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = "no")
    #
    # Generate MDVs from metabolic state and carbon source.
    #
    mdv = model.generate_mdv(flux, cs)
    data = mdv.get_fragment_mdv("Glue")
    print("Calced MDV data of              : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Answer by Antoniewiez et al 2007: [0.3464 0.2695 0.2708 0.0807 0.0286 0.0039]")

