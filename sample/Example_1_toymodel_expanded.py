#-------------------------------------------------------------------------------
# Name:        mfapy example 1 toymodel
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
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_1_toymodel_model.txt", format = "text")
    #
    # Construct MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(default_reaction_lb = 0.001) #
    model.set_configuration(default_reaction_ub = 5000.0) #
    model.set_configuration(default_metabolite_lb = 0.001) #
    model.set_configuration(default_metabolite_ub = 500.0) #
    model.set_configuration(default_reversible_lb = -300) #
    model.set_configuration(default_reversible_ub = 300.0) #
    model.set_configuration(iteration_max = 10000) # Maximal iternations in optimization
    model.set_configuration(initial_search_repeats_in_grid_search = 5) #
    model.set_configuration(initial_search_iteration_max = 1000) #
    model.set_configuration(grid_search_iterations = 1) #
    model.set_configuration(number_of_repeat = 3) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = 3) #Number of local CPUs for Parallel python
    #model.set_configuration(ppservers = ("172.16.0.149","172.16.0.133")) #IP address of PPworkers for Parallel python
    #
    # Load metabolite state from text file
    #
    state_dic = model.load_states("Example_1_toymodel_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    #
    # Generate ramdom metabolic state (flux + metabolic conc)
    #
    flux_opt, state = model.generate_flux_distribution()
    flux_opt_original = flux_opt
    #
    # Save metabolic state
    #
    #model.save_states(flux_opt, "Example_1_test.csv", format = 'csv')
    #
    # Set constraints and boundaries manusally
    # "fitting (variable and used for RSS calclation)", "fixed (fixed)", "free (variable and ignored for RSS calclation)", and "pseudo (fixed and used for MDV calc but ignored in the stoichiometry matrix)"
    model.set_constrain('reaction', 'v9', "fitting", 0.0, 1)
    model.set_boundary('reaction', 'v9', 0.0, 30.0)
    #
    # self.update() is required after manusal setting constraints and boundaries
    #
    model.update()
    #
    # Generate instances of CarbonSource class from model
    #
    cs = model.generate_carbon_source_templete()
    cs2 = model.generate_carbon_source_templete()
    #
    # Set isotope labelling of carbon sources by ratio of all isotopomer (#00, #01, #10, #11)
    #
    cs.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.25, 0.25])
    cs2.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.00, 0.50])
    #
    # Set isotope labelling of carbon sources by specific isotopomers
    #
    cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = "yes")
    cs2.set_each_isotopomer('Asp', {'#0000':0.5, '#1111':0.5})
    #
    # Generate MDVs from metabolic state and carbon source.
    #
    mdv1 = model.generate_mdv(flux_opt, cs)
    mdv2 = model.generate_mdv(flux_opt, cs2)
    #
    # Add gaussian noises for MDV data. Sample (n=3) gaussian noise with STD = 0.01, then calc their average and stdev
    #
    mdv1.add_gaussian_noise(0.01, 3)
    mdv2.add_gaussian_noise(0.01, 3)
    #
    # Set standard deviations
    #
    mdv1.set_std(0.01)
    mdv2.set_std(0.01)
    #
    # Ignore MDV data less than 0.01
    #
    mdv1.set_mdvs_for_comparison(0.01, 1.0)
    mdv2.set_mdvs_for_comparison(0.01, 1.0)
    #
    # Save MDV data
    #
    mdv1.save("testmdv1.txt", format = "text")
    #mdv2.save("testmdv2.txt", format = "text")
    #
    # Load MDV data from file
    #
    #mdv3 = model.load_mdv_data("testmdv1.txt", format = "text", output = "normal")
    #
    # Set experiments (Parallel labeling experiments)
    #

    model.set_experiment('ex1', mdv1, cs)
    model.set_experiment('ex2', mdv2, cs2)
    #
    # Set v8 and FUM as free
    #
    model.set_constrain('reaction','v8','free')
    model.set_constrain('reversible','FUM','free')
    #
    # Update model
    #
    model.update()

    #
    # Stored "template" metabolic state
    #
    results = [('template', flux_opt)]
    #
    # Generate initial metabolic state for fitting
    # Initial metabolic states were generate by 50 trials (sometimes failed), then return 1 initial metabolic state with smallest RSSs.
    state, flux_initial1 = model.generate_initial_states(50, 1)
    results.append(("Initial", flux_initial1))
    #
    # Test self.calc_rss()
    #
    print("Initial metabolic states were successfully obtained. RSS:", model.calc_rss(flux_opt))
    #
    # Find best fitted metabolic state by SLSQP (scipy)
    #
    start = time.time()
    state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = 'SLSQP', flux = flux_initial1)
    print("Fitting_flux(SLSQP scipy integrate): State", "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit, time.time()-start))
    results.append(("SLSQP", flux_opt_slsqp))
    #
    # Test goodness of fit
    #
    pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)
    print("RSS", RSS_bestfit,"Results: p-value:",pvalue, "thres", rss_thres)
    #
    # Calc threshold RSS for confidence interval search
    #
    thres, number_of_measurements, degree_of_freedom  =\
    model.get_thres_confidence_interval(flux_opt_slsqp, alpha = 0.05)

    print("Results: Threshold:",thres," Number of measurements:", number_of_measurements," Degree of freedom", degree_of_freedom)

    #
    # Find best fitted metabolic state by LN_PRAXIS (nlopt)
    #
    #test other
    #method = "SLSQP"
    method = "LN_PRAXIS"
    #method = "GN_CRS2_LM"
    start = time.time()
    state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = method, flux = flux_initial1)
    print("Fitting_flux(",method, "nlopt): State", "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit, time.time()-start))
    results.append((method, flux_opt_slsqp))
    #
    # Find best fitted metabolic state by deep (SLSQP=>LN_PRAXIS)*n. This iteration was set by self.set_configuration(number_of_repeat = 3)
    #
    start = time.time()
    state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = 'deep', flux = flux_initial1)
    print("Fitting_flux(deep): State", "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit, time.time()-start))
    results.append(("deep", flux_opt_slsqp))
    #
    # Generate initial metabolic state for fitting
    # Initial metabolic states were generate by 50 trials (sometimes failed), then return 4 initial metabolic state with smallest RSSs.
    # model.set_configuration(ncpus = 3) #Number of local CPUs for Parallel python
    #
    #
    state, flux_initial2 = model.generate_initial_states(50, 4)
    start = time.time()
    #
    # Find best fitted metabolic state by SLSQP using parallel mode
    #
    state, RSS_bestfit, flux_opt_slsqp2 = model.fitting_flux(method = 'SLSQP', flux = flux_initial2)
    print("Fitting_flux(SLSQP scipy integrate, parallel python): State", "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit[0], time.time()-start))
    #
    # Show results
    #
    model.show_results(results)
    #
    # MDV comparison
    #
    mdv_initial1 = model.generate_mdv(flux_initial1, cs)
    mdv_flux_opt_slsqp = model.generate_mdv(flux_opt_slsqp, cs)

    id_array, ratio_array1, std_array, use_array, observed_fragments, data = mdv1.generate_observed_mdv()
    id_array, ratio_array_initial1, std_array, use_array, observed_fragments, data = mdv_initial1.generate_observed_mdv()
    id_array, ratio_array_flux_opt_slsqp, std_array, use_array, observed_fragments, data = mdv_flux_opt_slsqp.generate_observed_mdv()

    fig = figure()
    axL = fig.add_subplot(121)
    axL.scatter(ratio_array_initial1, ratio_array1, c = 'r')
    axL.scatter(ratio_array_flux_opt_slsqp, ratio_array1,  c = 'b')
    #
    # Target reactions to estimate 95% CI
    #
    target_reactions = [('reaction', "v6"),('reversible', "FUM")]
    ci_edge = model.generate_ci_templete(targets = target_reactions)
    #
    # Grid search
    # model.set_configuration(initial_search_repeats_in_grid_search = 5) # Number of repeat for finding initial flux in Grid search mode such as model.generate_initial_states(5, 1)
    # model.set_configuration(grid_search_iterations = 1) # Fitting interations at each grid.
    #
    print("Test: Searching confidencial interval by grid search method.")
    model.set_configuration(callbacklevel = 1)
    ci = model.search_ci(ci_edge, flux_opt_slsqp, method = 'grid')

    for rid in ci['data'].keys():
        if ci['data'][rid]['use'] != 'on':continue
        print(rid, "Lower bondary:",ci['data'][rid]['lower_boundary'], "Upper boundary:", ci['data'][rid]['upper_boundary'])

    axR = fig.add_subplot(122)
    axR.scatter(ci['data'][('reversible', "FUM")]['flux_data'], ci['data'][('reversible', "FUM")]['rss_data'], c = 'r')
    axR.scatter(ci['data'][("reaction", "v6")]['flux_data'], ci['data'][("reaction", "v6")]['rss_data'], c = 'b')
    show()
