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
import mfapy as mfapy
from matplotlib.pyplot import figure, show
import os, sys, time

if __name__ == '__main__':
    #
    # Load metabolic model from txt file to four dictionary
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_3_C_model.txt", format = "text")
    #
    # Construct MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode="normal")
    #
    # Addition of isotope effect durint MDV calculation
    #
    model.set_configuration(add_naturalisotope_in_calmdv = "yes") #
    #
    # Requires model.reconstruct()
    #
    model.reconstruct()
    #

    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(iteration_max = 20000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 5) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    #
    # Load metabolic state from file
    #
    state_dic = model.load_states("Example_3_C_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Set isotope labelling of carbon sources 1
    #
    carbon_source1 = model.generate_carbon_source_templete()
    carbon_source1.set_each_isotopomer('GLCEX',{'#000000': 0.0, '#100000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('CO2CS',{'#0': 1.0, '#1': 0.0}, correction = 'yes')    #
    # Load measured MDV data 1
    #
    mdv_observed1 = model.load_mdv_data('Example_3_C_mdv.txt')
    mdv_observed1.set_std(0.01, method = 'absolute')
    #mdv_observed1.correct_natural_isotope(mode = 'correction')
    #mdv_observed1.set_mdvs_for_comparison(0.001)
    #
    # Set experiment 1
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)    #

    state, flux_opt1 = model.generate_initial_states(50, 1)
    print("RSS:",model.calc_rss(flux_opt1))
    results = [('template', flux_opt1)]

    model.set_configuration(callbacklevel = 8) #
    #pvalue, rss_thres = model.goodness_of_fit(flux_opt1, alpha = 0.05)

    #model.show_results(results)

    #
    # Fitting by local optimizer
    #
    method = "SLSQP"
    start = time.time()
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
    pvalue, rss_thres = model.goodness_of_fit(flux_opt1, alpha = 0.05)
    print(method, ": State", "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
    results.append((method, flux_opt1))
    #
    # Fitting by global optimizer
    #
    method = "GN_CRS2_LM"
    start = time.time()
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
    pvalue, rss_thres = model.goodness_of_fit(flux_opt1, alpha = 0.05)
    print(method, ": State", "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
    results.append((method, flux_opt1))

    #
    # Fitting by local + flobal optimizer
    #
    method = "deep"
    start = time.time()
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
    pvalue, rss_thres = model.goodness_of_fit(flux_opt1, alpha = 0.05)
    print(method, ": State", "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
    results.append((method, flux_opt1))
    model.show_results(results)





