#-------------------------------------------------------------------------------
# Name:        mfapy Example 2-2  3C-MFA of metabolically enginnered E. coli
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

    numofmodel = 20
    repeatn = 5
    #
    # Construction of metabolic model
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Example_2_Ecoli_model.txt", format = "text")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode = "normal")
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(iteration_max = 10000)
    model.set_configuration(ncpus = 4) #Number of local CPUs for Parallel python
    #
    # Batch setting of constrains from  file
    #
    state_dic = model.load_states("Example_2_Ecoli_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(state_dic)
    model.update()
    #
    # Set isotope labelling of carbon sources
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
    mdv_observed1.set_std(0.004, method = 'absolute')
    #
    # Set experiment
    #
    model.set_experiment('ex1', mdv_observed1, carbon_source1)
    #
    #
    print("Obtaining initial state using parallel processing using joblib")
    flux_initial_array = []
    rss_initial_array = []
    for j in range(repeatn):
        state, flux_initial = model.generate_initial_states(numofmodel*10, numofmodel, method ="parallel")
        flux_initial_array.append(flux_initial)
        rss_initial_array.extend(model.calc_rss(flux_initial))
    temparray = sorted(rss_initial_array)
    datan = len(temparray) -1
    print("initial", datan, "RSS:{0:>8.2f}  RSS:{1:>8.2f} RSS:{2:>8.2f} RSS:{3:>8.2f} RSS:{4:>8.2f}".format(temparray[0], temparray[int(datan*0.05)], temparray[int(datan*0.5)], temparray[int(datan*0.95)], temparray[datan]))

    model.set_configuration(callbacklevel = 0) #
    #
    # Fitting by deep optimizer
    #
    RSS_dict = {}

    for method in ["SLSQP","COBYLA","LN_COBYLA","LN_BOBYQA","LN_PRAXIS","LN_SBPLX","LN_NELDERMEAD","GN_DIRECT_L","GN_AGS","GN_CRS2_LM","GN_ESCH","GN_IRES","LN_NEWUOA"]:

        RSS_dict[method] = []
        for j in range(repeatn):
            start = time.time()
            flux_opt1 = flux_initial_array[j]
            for i in range(1000):
                state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
                if len(flux_opt1) == 0:
                    break
                pvalue, rss_thres = model.goodness_of_fit(flux_opt1[0], alpha = 0.05)
                now = time.time()
                if now - start > 60:
                    break
            RSS_dict[method].extend(RSS_bestfit)
        RSS_dict[method].sort()
        temparray = RSS_dict[method]
        if len(temparray) == 0:
            continue
        datan = len(temparray) -1

        print(method, datan, "RSS:{0:>8.2f}  RSS:{1:>8.2f} RSS:{2:>8.2f} RSS:{3:>8.2f} RSS:{4:>8.2f}".format(temparray[0], temparray[int(datan*0.05)], temparray[int(datan*0.5)], temparray[int(datan*0.95)], temparray[datan]))




