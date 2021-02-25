#-------------------------------------------------------------------------------
# Name:        mfapy example 7 Monte Carlo Simulation
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import os, sys, time, copy, csv
import mfapy
import scipy.stats as stats

if __name__ == '__main__':
    #
    #
    metabolicmodel = "Example_7_CancerCell_model.txt"
    metabolicstate = "Example_7_flux_dist_CancerCell.csv" # This file is used to set model constraints
    mdvfile1 = 'Example_7_MCF7_C24_U_CIT.txt'
    mdvfile2 = 'Example_7_MCF7_T24_U_CIT.txt'
    numberofsteps = 5000
    numberofcpus = 2
    #
    # Construction of MetabolicModel instance
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model(metabolicmodel)
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0)
    model.set_configuration(iteration_max = 10000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 3) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = numberofcpus) #Number of local CPUs for Parallel computing
    #
    # Load metabolite state from text file
    #
    state_dic = model.load_states(metabolicstate, format = 'csv')
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
    carbon_source1.set_each_isotopomer('SubsGlc',{'#000000': 1.0}, correction = 'yes')
    carbon_source1.set_each_isotopomer('SubsGln',{'#11111': 1.0}, correction = 'no')
    #
    # For two mdv files
    #
    for mdvfile in [mdvfile1, mdvfile2]:

        mdv_observed1 = model.load_mdv_data(mdvfile)
        model.set_experiment('ex1', mdv_observed1, carbon_source1)

        state, flux_opt1 = model.generate_initial_states(2, 1)
        rss = model.calc_rss(flux_opt1)
        method = "GN_CRS2_LM"
        state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
        for method in ["deep", "GN_CRS2_LM", "deep"]:
            start = time.time()
            state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
        model.show_results([("test",flux_opt1)])


        rss = model.calc_rss(flux_opt1)
        df = mdv_observed1.get_number_of_measurement()
        sf = stats.chi2.sf(x = rss, df = df)
        pdf = stats.chi2.pdf(x = rss, df = df)

        print(mdvfile, "RSS:",rss, "p-value of chisq test", sf, "pdf", pdf)
        if sf < 0.95:
            print("The metabolic model failed to overfit to isotopoer data", mdvfile)
            continue
        #
        # Estimation of posterior distribution
        #
        record = model.posterior_distribution(flux_opt1, number_of_steps = numberofsteps)
        outputfilename = mdvfile+".csv"


        with open(outputfilename, 'w') as f:
            writer = csv.writer(f, lineterminator='\n') #く
            writer.writerows(record)






