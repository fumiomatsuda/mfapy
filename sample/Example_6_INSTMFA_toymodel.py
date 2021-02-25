#-------------------------------------------------------------------------------
# Name:        Example_5_INSTMFA_toymodel.py
#              Example 5 INST-MFA example using a toy model.
#              Artificial Time course MDV data was genenrate and used for flux estimation
#              Toy model used in this code in modified from Antoniewicz et al Metab. Eng. 2007, 9, 68-86.
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import mfapy
from matplotlib.pyplot import figure, show

if __name__ == '__main__':
    #
    # Construction of MetabolicModel instance
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Example_4_Simulation_model.txt", format = "text")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 1) #
    model.set_configuration(iteration_max = 1000) # Maximal iternations in optimization
    model.set_configuration(initial_search_iteration_max = 1000) #
    model.set_configuration(number_of_repeat = 3) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = 3) #Number of local CPUs for Parallel python
    #
    # Loading metabolite state from text file and batch setting of constraints
    #
    state_dic = model.load_states("Example_4_Simulation_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Generation of CarbonSource instance and isotope labelling information
    #
    cs = model.generate_carbon_source_template()
    cs.set_each_isotopomer('SubsGlc', {'#000000':0.5,'#111111':0.5}, correction = 'yes')
    #
    # Generation of initial IDV state (IDV at experimenta start (time = 0)
    #
    cs2 = model.generate_carbon_source_template()
    cs2.set_each_isotopomer('SubsGlc', {'#000000':1.0}, correction = 'yes')
    idv = model.calc_idv(state_dic, cs2)
    #
    # Generation of artificial time course MDV data
    #
    timepoints = [0.0,5.0,10.0,15.0,30.0,50.0]
    mdv_timecourse = model.generate_mdv(state_dic, cs, timepoints, startidv = idv)
    #
    # Addition of Gaussian noise, standard deviation and removal of too small signals
    #
    mdv_timecourse.add_gaussian_noise(0.01, 3)
    mdv_timecourse.set_std(0.01, method = 'absolute')
    mdv_timecourse.set_mdvs_for_comparison(0.01, 1.0)
    #
    # "mdv_timecourse.save("Example_1_toymodel_mdvtimecourse") generate following files
    #
    # Example_5_INST_mdvtimecourse0.txt
    # Example_5_INST_mdvtimecourse5.txt
    # Example_5_INST_mdvtimecourse10.txt
    # Example_5_INST_mdvtimecourse15.txt
    # Example_5_INST_mdvtimecourse30.txt
    # Example_5_INST_mdvtimecourse50.txt
    #
    mdv_timecourse.save("Example_6_INST_mdvtimecourse")
    #
    # Loading of mdv data
    #
    mdv_00 = model.load_mdv_data('Example_6_INST_mdvtimecourse0.txt')
    mdv_05 = model.load_mdv_data('Example_6_INST_mdvtimecourse5.txt')
    mdv_10 = model.load_mdv_data('Example_6_INST_mdvtimecourse10.txt')
    mdv_15 = model.load_mdv_data('Example_6_INST_mdvtimecourse15.txt')
    mdv_30 = model.load_mdv_data('Example_6_INST_mdvtimecourse30.txt')
    mdv_50 = model.load_mdv_data('Example_6_INST_mdvtimecourse50.txt')
    #
    # Construction of time course mdv data from loaded mdv data using MdvTimeCourseData instance
    #
    mdv = mfapy.mdv.MdvTimeCourseData()
    mdv.add_time_point(0, mdv_00)
    mdv.add_time_point(5.0, mdv_05)
    mdv.add_time_point(10, mdv_10)
    mdv.add_time_point(15, mdv_15)
    mdv.add_time_point(30, mdv_30)
    mdv.add_time_point(50, mdv_50)
    mdv.set_std(0.01, method = 'absolute')
    #
    # Setting of labeling experiment. INST experiment requires "startidv" information
    #
    model.set_experiment('tc', mdv, cs, startidv = idv)
    #
    # Model Fitting in parallel mode
    #
    print("Fitting is started in parallel mode.")
    state, flux_opt = model.generate_initial_states(50, 2)
    results = [('initial', flux_opt[0])]
    model.set_configuration(iteration_max = 1) # For quick fitting (test purpose)
    model.set_configuration(number_of_repeat = 1) # For quick fitting (test purpose)
    state, RSS_bestfit, flux_opt_parallel = model.fitting_flux(method = 'deep', flux = flux_opt)
    results.extend([('SLSQP', flux) for (i,flux) in enumerate(flux_opt_parallel)])
    #
    # Show results
    #
    model.show_results(results)
    #
    # Show MDV comparison
    #
    mdv_predicted = model.generate_mdv(flux_opt_parallel[0], cs, timepoints, startidv = idv)
    id_array, ratio_array1, std_array, use_array, observed_fragments, data = mdv_timecourse.generate_observed_mdv()
    id_array, ratio_array_predicted, std_array, use_array, observed_fragments, data = mdv_predicted.generate_observed_mdv()

    fig = figure()
    axL = fig.add_subplot(121)
    axL.scatter(ratio_array_predicted, ratio_array1,  c = 'b')
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






