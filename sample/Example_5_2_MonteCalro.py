#-------------------------------------------------------------------------------
# Name:        mfapy example 5-2  Monte Carlo Simulation
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
    model.set_configuration(ncpus = 20) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    #
    # Load metabolite state from text file
    #
    state_dic = model.load_states(metabolicstate, format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Set isotope labelling of carbon sources 1
    #
    cs1 = model.generate_carbon_source_template()
    cs1.set_each_isotopomer('SubsGlc',{'#111111':0.5,'#000000':0.5}, correction = "no")
    #
    # Estimation of
    #
    for fragment, numbers in [("PEPe", [0,1,2,3]), ("OACe", [0,1,2,3,4]), ("AKGe", [0,1,2,3,4,5])]:
        for number in numbers:
            outputfilename = 'Example_6_'+fragment+str(number)+'.csv'
            #
            # Loading of MDV data and Setting of MDV data used for RSS calculation
            #
            mdv_observed1 = model.generate_mdv(state_dic, cs1)
            mdv_observed1.set_std(0.01, method = 'absolute')
            mdv_observed1.set_observed_fragments([fragment])
            for i in numbers:
                mdv_observed1.set_unused_mdv_for_comparison(fragment, i)
            mdv_observed1.set_mdv_for_comparison(fragment, number)
            #
            # Submittion of lableling experiment
            #
            model.clear_experiment()
            model.set_experiment('ex1', mdv_observed1, cs1)
            #
            # Finding a best fitted metabolic state
            #
            state, flux_opt1 = model.generate_initial_states(2, 1)
            for method in ["SLSQP","LN_SBPLX","SLSQP","LN_SBPLX"]:
                state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)

            df = mdv_observed1.get_number_of_measurement()
            sf = stats.chi2.sf(x = RSS_bestfit, df = df)
            #
            # Check over-fitting
            #
            if sf < 0.95:
                pass
                #print("The metabolic model failed to overfit to the isotopoer data. p-value:", sf)
            #
            # Finding a good initial metabolic state with suitable rss level
            #
            rss = 0
            while rss < 0.01:
                check, perturbed_state, flux_vector, reaction_name_list = model.add_perturbation(flux_opt1)
                rss = model.calc_rss(perturbed_state)
                pdf = stats.chi2.pdf(x = rss, df = df)
            #
            # Estimation of posterior distribution
            #
            record = model.posterior_distribution(perturbed_state, number_of_steps = numberofsteps)
            start = int(numberofsteps*model.configuration['ncpus']/1000/2)
            data = []
            for row in record:
                data.append(row[4])
            data = data[start:]
            data.sort()
            p025 = int(len(data)*0.025)
            p975 = int(len(data)*0.975)
            print("v4","\t",outputfilename,"\t",numpy.mean(data),"\t",numpy.std(data),"\t",data[p025],"\t",data[p975])

            data = []
            for row in record:
                data.append(row[5])
            data = data[start:]
            data.sort()
            p025 = int(len(data)*0.025)
            p975 = int(len(data)*0.975)
            print("v5","\t",outputfilename,"\t",numpy.mean(data),"\t",numpy.std(data),"\t",data[p025],"\t",data[p975])

            data = []
            for row in record:
                data.append(row[6])
            data = data[start:]
            data.sort()
            p025 = int(len(data)*0.025)
            p975 = int(len(data)*0.975)
            print("v6","\t",outputfilename,"\t",numpy.mean(data),"\t",numpy.std(data),"\t",data[p025],"\t",data[p975])

            data = []
            for row in record:
                data.append(row[7])
            data = data[start:]
            data.sort()
            p025 = int(len(data)*0.025)
            p975 = int(len(data)*0.975)
            print("v7","\t",outputfilename,"\t",numpy.mean(data),"\t",numpy.std(data),"\t",data[p025],"\t",data[p975])


            with open(outputfilename, 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #く
                writer.writerows(record)






