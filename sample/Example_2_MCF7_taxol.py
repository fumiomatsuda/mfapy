#-------------------------------------------------------------------------------
# Name:        mfapy example 2 Taxol treated MCF-7 cells
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import os, sys, time
import numpy as numpy
import mfapy.mfapy as mfapy

if __name__ == '__main__':

    #
    # Load metabolic model from txt file to four dictionary
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Example_2_MCF7_taxol_model.txt")
    #
    # Construct MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0)
    model.set_configuration(iteration_max = 20000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 5) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = 2) #Number of local CPUs for Parallel python

    stdev = 0.015
    #
    # Load metabolite state from text file
    #
    state_dic = model.load_states("Example_2_MCF7_taxol_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Generate instances of CarbonSource class from model
    #
    carbon_source1 = model.generate_carbon_source_templete()
    #
    # Set isotope labelling of carbon sources 1
    #
    carbon_source1 = model.generate_carbon_source_templete()
    carbon_source1.set_all_isotopomers('SubsCO2', [0.99, 0.01])
    carbon_source1.set_each_isotopomer('SubsGlc',{'#000000': 0, '#100000': 1}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsGln',{'#00000': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsSer',{'#000': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsGly',{'#00': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsLeu',{'#000000': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsIle',{'#000000': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsVal',{'#00000': 1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsArg',{'#000000':1.}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsCys',{'#000':1.}, correction = 'yes')
    #
    # Generate instances of CarbonSource class from model
    #
    carbon_source2 = model.generate_carbon_source_templete()
    #
    # Set isotope labelling of carbon sources 2
    #
    carbon_source2 = model.generate_carbon_source_templete()
    carbon_source2.set_all_isotopomers('SubsCO2', [0.99, 0.01])
    carbon_source2.set_each_isotopomer('SubsGlc',{'#000000': 1, '#100000': 0}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsGln',{'#11111': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsSer',{'#000': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsGly',{'#00': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsLeu',{'#000000': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsIle',{'#000000': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsVal',{'#00000': 1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsArg',{'#000000':1}, correction = 'yes')
    carbon_source2.set_each_isotopomer('SubsCys',{'#000':1}, correction = 'yes')
    #
    # Load measured MDV data 1
    #
    mdv_observed1 = model.load_mdv_data('Example_2_MCF7_taxol_mdv_1_13C_Glc.txt')
    mdv_observed1.set_std(stdev, method = 'absolute')
    #
    # Set experiment 1
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)    #
    #
    # Load measured MDV data 2
    #
    mdv_observed2 = model.load_mdv_data('Example_2_MCF7_taxol_mdv_U_13C_Gln.txt')
    mdv_observed2.set_std(stdev, method = 'absolute')
    #
    # Set experiment 2
    #
    model.set_experiment('ex2', mdv_observed2, carbon_source2)
    #
    # Fitting by various solvers
    #
    model.set_configuration(callbacklevel = 0)
    #
    # Generate initial metabolic state for fitting
    #
    state, flux_opt1 = model.generate_initial_states(50, 1)
    print("RSS:",model.calc_rss(flux_opt1))
    results = [('template', flux_opt1)]
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
    # Additional fitting by local optimizer
    #
    #methods = ["SLSQP", "COBYLA", "LN_COBYLA", "LN_BOBYQA", "LN_NEWUOA", "LN_PRAXIS", "LN_SBPLX", "LN_NELDERMEAD", "GN_CRS2_LM", "deep"]
    for method in ["SLSQP"]:
        start = time.time()
        state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = method, flux = flux_opt1)
        pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)

        print(method, ": State", "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
        results.append((method, flux_opt_slsqp))
    model.show_results(results, pool_size = "off")



