#-------------------------------------------------------------------------------
# Name:        mfapy example 5-1 Monte Carlo Simulation
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import os, sys, time, copy, csv
import numpy as numpy
import mfapy
import scipy.stats as stats

if __name__ == '__main__':
    metabolicmodel = "Example_4_Simulation_model.txt"
    metabolicstate = "Example_4_Simulation_status.csv" # This file is used to set model constraints
    outputfilename = 'Example_5_output.csv'
    numberofsteps = 5000000
    #
    # Construciton of metabolic model
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model(metabolicmodel)
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0)
    model.set_configuration(iteration_max = 10000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 3) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    #
    # Loading of metabolite state from text file
    #
    state_dic = model.load_states(metabolicstate, format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Setting carbon source1
    #
    cs1 = model.generate_carbon_source_template()
    cs1.set_each_isotopomer('SubsGlc',{'#111111':0.5,'#000000':0.5}, correction = "no")
    #
    # Set experiment
    #
    mdv_observed1 = model.generate_mdv(state_dic, cs1)
    mdv_observed1.set_observed_fragments(["OAC"])
    mdv_observed1.set_mdv_for_comparison("OAC",3)
    model.set_experiment('ex1', mdv_observed1, cs1)
    #
    # Finding a initial step
    #
    state, flux_opt1 = model.generate_initial_states(2, 1)
    for method in ["SLSQP","LN_SBPLX"]:
        state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)

    rss = model.calc_rss(flux_opt1)
    if rss < 0.001:
        rss = 0.001
    df = mdv_observed1.get_number_of_measurement()
    sf = stats.chi2.sf(x = rss, df = df)
    pdf = stats.chi2.pdf(x = rss, df = df)

    print(df, "RSS:",rss, "p-value of chisq test", sf, "pdf", pdf)
    if sf < 0.95:
        print("The metabolic model failed to overfit to the isotopoer data", mdvfile)

    independent_vector, lbi, ubi, independent_list = model.get_independents(flux_opt1)
    #
    # Data store
    #
    record = []
    record_temp = [rss]
    check, flux_vector, reaction_name = model.check_independents(independent_vector)
    record_temp.extend(reaction_name)
    record.append(record_temp)


    Rm_ind = independent_vector[:]
    Rm_ind_record = [Rm_ind] * 1000

    progress = 0
    unaccepted = 0

    while progress < chain_length:
        Rm_ind_next = Rm_ind[:]
        for j in range(100):
            reac_num = numpy.random.randint(0, len(independent_vector))
            Rm_ind_next_temp = Rm_ind_next[:]
            perturbation =  (numpy.random.rand() - 0.5) * 2 * (ubi[reac_num]-lbi[reac_num])/100
            Rm_ind_next_temp[reac_num] += perturbation
            if lbi[reac_num] < Rm_ind_next_temp[reac_num] < ubi[reac_num]:
                check, flux_vector, reaction_name = model.check_independents(Rm_ind_next_temp)
                if check == True:
                    Rm_ind_next[reac_num] = Rm_ind_next[reac_num] + perturbation
                    break
        else:
            print("Can't find next step. Return to 1000 steps before.")
            Rm_ind = Rm_ind_record[0]
            continue
        rss_next = model.calc_rss(Rm_ind_next, mode = "independent")
        pdf_next = stats.chi2.pdf(x = rss_next, df = df)
        #
        # If probabirity of next point is smaller than that of present point
        #
        if pdf_next < pdf:
            if numpy.random.rand() > pdf_next/pdf:
                    progress = progress + 1
                    unaccepted = unaccepted + 1
                    Rm_ind_record.pop(0)
                    Rm_ind_record.append(Rm_ind)
                    continue
        Rm_ind[:] = Rm_ind_next[:]
        rss = rss_next
        pdf = pdf_next
        Rm_ind_record.pop(0)
        Rm_ind_record.append(Rm_ind)
        progress = progress + 1
        if progress % 1000 == 0:
            check, flux_vector, reaction_name = model.check_independents(Rm_ind)
            print(progress, rss, unaccepted, progress,Rm_ind)
            record_temp = [rss]
            record_temp.extend(flux_vector)
            record.append(record_temp)

    with open(outputfilename, 'w') as f:
        writer = csv.writer(f, lineterminator='\n') #く
        writer.writerows(record)







