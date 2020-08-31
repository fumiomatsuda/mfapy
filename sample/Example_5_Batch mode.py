#-------------------------------------------------------------------------------
# Name:        mfapy example 5 Batch execution mode (Example2 is used for testing: Taxol treated MCF-7 cells)
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import os, sys, time
import mfapy as mfapy


if __name__ == '__main__':
    #
    # The four files are defined here
    #
    filename_metabolicmodel = "Example_2_MCF7_taxol_model.txt"
    filename_metabolicstate = "Example_2_MCF7_taxol_status.csv"
    filename_carbonsource = "Example_2_carbonsource1.txt"
    filename_mdvdata = "Example_2_MCF7_taxol_mdv_1_13C_Glc.txt"
    filename_output = "Example_5_output.csv"
    def_stdev = 0.015
    def_ncpu = 4
    #
    #
    # Construction of MetabolicModel instance
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model(filename_metabolicmodel)
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0)
    model.set_configuration(iteration_max = 20000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 5) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = def_ncpu) #Number of local CPUs for Parallel python
    stdev = def_stdev
    #
    # Loading metabolite state from text file and atch constraint setting
    #
    state_dic = model.load_states(filename_metabolicstate, format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Generation of CarbonSource instance
    #
    carbon_source1 = model.generate_carbon_source_templete()
    carbon_source1 = model.generate_carbon_source_templete()
    carbon_source1.set_set_carbonsources(filename_carbonsource, correction = 'yes')
    #
    # Loading of measured MDV data 1
    #
    mdv_observed1 = model.load_mdv_data(filename_mdvdata )
    mdv_observed1.set_std(stdev, method = 'absolute')
    #
    # Addition of labeling experiment 1
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)    #
    #
    # Generate initial metabolic state
    #
    state, flux_opt1 = model.generate_initial_states(50, 1, method ="parallel")
    print("RSS of initial state is:",model.calc_rss(flux_opt1))
    results = [('template', flux_opt1)]
    #
    # Fitting by global optimizer
    #
    method = "GN_CRS2_LM"
    start = time.time()
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
    pvalue, rss_thres = model.goodness_of_fit(flux_opt1, alpha = 0.05)
    print(method, ": State",state,  "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
    results.append((method, flux_opt1))
    #
    # Additional fitting by local optimizer
    #
    #methods = ["SLSQP", "COBYLA", "LN_COBYLA", "LN_BOBYQA", "LN_NEWUOA", "LN_PRAXIS", "LN_SBPLX", "LN_NELDERMEAD", "GN_CRS2_LM", "deep"]
    for method in ["LN_PRAXIS"]:
        start = time.time()
        state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = method, flux = flux_opt1)
        pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)

        print(method, ": State", state, "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
        results.append((method, flux_opt_slsqp))
    #
    # Output fitted metabolic flux in console
    #
    model.show_results(results, filename = filename_output, pool_size = "off",  checkrss = "off")


