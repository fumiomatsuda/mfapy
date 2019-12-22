#-------------------------------------------------------------------------------
# Name:        mfapy example 1 INST-13C-MFA toymodel
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
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_1_toymodel_model.txt", format = "text")
    #
    # Construct MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(iteration_max = 1000) # Maximal iternations in optimization
    model.set_configuration(initial_search_iteration_max = 1000) #
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
    #
    # Set constraints and boundaries manusally
    # "fitting (variable and used for RSS calclation)", "fixed (fixed)", "free (variable and ignored for RSS calclation)", and "pseudo (fixed and used for MDV calc but ignored in the stoichiometry matrix)"
    model.set_constrain('reaction', 'v9', "fitting", 0.0, 1)
    model.set_boundary('reaction', 'v9', 0.0, 30.0)
    model.set_constrain('reaction','v8','free')
    model.set_constrain('reversible','FUM','free')
    model.set_constrain('metabolite', 'OA', "fitting", 10.0,0.5)
    #
    # self.update() is required after manusal setting constraints and boundaries
    #
    model.update()
    #
    # Generate instances of CarbonSource class from model
    #
    cs = model.generate_carbon_source_templete()
    #
    # Set isotope labelling of carbon sources
    #
    cs.set_all_isotopomers('AcCoA', [0.0, 0.0, 0.0, 1.0])
    cs.set_each_isotopomer('Asp', {'#1111':1.0}, correction = 'no')
    cs.set_each_isotopomer('OACs', {'#1111':1.0})
    #
    # Generate starting IDV state ( IDV at time = 0)
    #
    cs2 = model.generate_carbon_source_templete()
    cs2.set_all_isotopomers('AcCoA', [1.0, 0.0, 0.0, 0.0])
    cs2.set_each_isotopomer('Asp', {'#0000':1.0})
    cs2.set_each_isotopomer('OACs', {'#0000':1.0})
    idv = model.calc_idv(flux_opt, cs2)
    #
    #
    # Generate time course MDV data
    #
    timepoints = [0.0,5.0,10.0,15.0,30.0,50.0]
    mdv_timecourse = model.generate_mdv(flux_opt, cs, timepoints, startidv = idv)

    mdv_timecourse.add_gaussian_noise(0.01, 3)
    mdv_timecourse.set_std(0.01, method = 'absolute')
    mdv_timecourse.set_mdvs_for_comparison(0.01, 1.0)
    #
    # Set experiment
    #
    model.set_experiment('tc', mdv_timecourse, cs, startidv = idv)

    rss = model.calc_rss(flux_opt)
    print("RSS of best fitted metabolic state seems to be:", rss)

    #
    # Initial metabolic state generation
    #
    model.set_configuration(callbacklevel = 7)
    print("Initial state generation.")
    start = time.time()
    results = [('start', flux_opt)]
    state, flux_start = model.generate_initial_states(50, 1)
    rss = model.calc_rss(flux_start)
    print("RSS of initial metabolic state seems to be:", rss)

    #
    # Fitting
    #

    model.set_configuration(iteration_max = 1) # For quick test purpose
    model.set_configuration(number_of_repeat = 1) # For quick test purpose

    print("Fitting is started.")
    results.append(('start', flux_start))
    state, RSS_bestfit, flux_opt_slsqp= model.fitting_flux(method = 'deep', flux = flux_start)
    results.append(('deep', flux_opt_slsqp))
    print(model.calc_rss(flux_opt_slsqp), RSS_bestfit)
    #
    # Fitting in parallel mode
    #
    print("Fitting is started in parallel mode.")
    state, flux_opt = model.generate_initial_states(50, 2)
    state, RSS_bestfit, flux_opt_parallel = model.fitting_flux(method = 'SLSQP', flux = flux_opt)
    results.extend([('deep_para', flux) for (i,flux) in enumerate(flux_opt_parallel)])
    #
    # Show results
    #
    model.show_results(results)
    #model.show_results(results, filename = "testresult.txt", format = "csv")
    #
    # chi-square test
    #
    pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)
    print("p-value:",pvalue, "thres", rss_thres)
    #
    # Calc threshold level
    #
    thres, number_of_measurements, degree_of_freedom =\
    model.get_thres_confidence_interval(flux_opt_slsqp, alpha = 0.05)
    print("Threshold:",thres," number of measurements:", number_of_measurements," degree of freedom", degree_of_freedom)
    #
    # Show MDV comparison
    #

    mdv_initial1 = model.generate_mdv(flux_start, cs, timepoints, startidv = idv)
    id_array, ratio_array1, std_array, use_array, observed_fragments, data = mdv_timecourse.generate_observed_mdv()
    id_array, ratio_array_initial1, std_array, use_array, observed_fragments, data = mdv_initial1.generate_observed_mdv()

    fig = figure()
    axL = fig.add_subplot(121)
    axL.scatter(ratio_array_initial1, ratio_array1,  c = 'b')

    #
    # Show MDV timecourse
    #
    axR = fig.add_subplot(122)
    tc = mdv_timecourse.get_timepoints()
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 0),  c = 'r')
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 1),  c = 'b')
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 2),  c = 'y')
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 3),  c = 'g')
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 4),  c = 'y')
    axR.plot(tc, mdv_timecourse.get_timecourse_data("AKGe", 5),  c = 'g')
    show()





