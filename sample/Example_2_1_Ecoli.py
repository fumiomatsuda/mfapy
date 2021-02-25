#-------------------------------------------------------------------------------
# Name:        mfapy Example 2-1  3C-MFA of metabolically enginnered E. coli
#              Metabolic model and data used in this code is derived from Okahashi et al. Biotechnol. Bioeng. 2017, 114, 2782-2793.

#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import mfapy
import os, sys, time

if __name__ == '__main__':
    #
    # Construction of metabolic model
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Example_2_Ecoli_model.txt", format = "text", mode = "normal")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode = "normal")
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(iteration_max = 10000)
    model.set_configuration(ncpus = 4) #Number of local CPUs for joblib
    #
    # Addtion of constrains
    #
    state_dic = model.load_states("Example_2_Ecoli_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Isotope labelling of carbon sources
    #
    carbon_source1 = model.generate_carbon_source_templete()
    carbon_source1.set_each_isotopomer('SubsGlc',{'#000000': 0.02, '#100000': 0.7, '#111111': 0.28 }, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsCO2',{'#0': 1.0, '#1': 0.0}, correction = 'no')
    carbon_source1.set_each_isotopomer('SubsAla',{'#000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsAsp',{'#0000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsThr',{'#0000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsGlu',{'#00000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsPhe',{'#000000000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsTyr',{'#000000000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsIle',{'#000000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsLeu',{'#000000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsVal',{'#00000': 1.0}, correction = 'yes')
    #
    # Load measured MDV data
    #
    mdv_observed1 = model.load_mdv_data('Example_2_Ecoli_mdv.txt', )
    mdv_observed1.set_std(0.01, method = 'absolute')
    #
    # Register experiment
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)
    #
    # Generation of initial flux
    #
    print("Generating 100 random initial states using parallel processing")
    state, flux_opt1 = model.generate_initial_states(500, 100, method ="parallel")
    print("RSS of the best initial state is:",model.calc_rss(flux_opt1[0]))
    results = [('template', flux_opt1[0])]
    model.set_configuration(callbacklevel = 0) #
    #
    # Model Fitting
    #
    print("Fitting by non-liner optimizations")
    start = time.time()
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = "GN_CRS2_LM", flux = flux_opt1)
    for i in range(15):
        for method in ["SLSQP","LN_SBPLX"]:
            state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
            pvalue, rss_thres = model.goodness_of_fit(flux_opt1[0], alpha = 0.05)
            print("RSS of the best state is:",model.calc_rss(flux_opt1[0]))
        results.append((method, flux_opt1[0]))
    model.show_results(results)








