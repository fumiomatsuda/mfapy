#-------------------------------------------------------------------------------
# Name:        mfapy Example_4_Simulation.py
#              Simulation of 13C-MFA using artificial observed MDVs to select effective carbon source.
#              Toy model used in this code in modified from Antoniewicz et al Metab. Eng. 2007, 9, 68-86.
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import mfapy

if __name__ == '__main__':
    #
    # Construction of metabolic model
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Example_4_Simulation_model.txt", format = "text")
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    model.set_configuration(iteration_max = 20000) # Maximal iternations in optimization
    model.set_configuration(number_of_repeat = 5) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = 4)
    model.set_configuration(initial_search_repeats_in_grid_search = 5) # Number of repeat for finding initial flux in Grid search mode such as model.generate_initial_states(5, 1)
    model.set_configuration(grid_search_iterations = 1) # Fitting interations at each grid.
    #
    # Load metabolic state from file
    #
    flux = model.load_states("Example_4_Simulation_status.csv", format = 'csv')
    model.set_constraints_from_state_dict(flux)
    for method, isotope1, isotope2  in [("single",{'#100000':1.0},{}),("single",{'#111111':0.5,'#000000':0.5},{}),("single",{'#111111':0.5,'#100000':0.5},{}),("parallel",{'#111111':0.5,'#000000':0.5},{'#100000':1.0})]:
        for i in range(100):
            #
            # Generation of instances of CarbonSource class from model
            #
            cs = model.generate_carbon_source_templete()
            cs.set_each_isotopomer('SubsGlc',isotope1, correction = "no")
            #
            # Generation of MDV instances from metabolic state and carbon source.
            #
            mdv = model.generate_mdv(flux, cs)
            mdv.add_gaussian_noise(0.01)
            mdv.set_std(0.01, method = 'absolute')
            #
            # Addition of labeling experiment 1
            #
            model.clear_experiment()
            model.set_experiment('ex1', mdv, cs)
            if method == "parallel":
                #
                # Generation of instances of CarbonSource class from model
                #
                cs2 = model.generate_carbon_source_templete()
                cs2.set_each_isotopomer('SubsGlc',isotope2, correction = "no")
                #
                # Generation of MDV instances from metabolic state and carbon source.
                #
                mdv2 = model.generate_mdv(flux, cs2)
                mdv2.add_gaussian_noise(0.01)
                mdv2.set_std(0.01, method = 'absolute')
                #
                # Addition of labeling experiment 2
                #
                model.set_experiment('ex2', mdv2, cs2)
            #
            # Generate initial metabolic state for fitting
            #
            model.set_configuration(iteration_max = 20000) # Maximal iternations in optimization
            model.set_configuration(number_of_repeat = 5) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
            state, flux_initial = model.generate_initial_states(200, 8, method ="parallel")
            state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = "GN_CRS2_LM", flux = flux_initial)
            state, RSS_bestfit, flux_opt2 = model.fitting_flux(method = "deep", flux = flux_opt1)
            pvalue, rss_thres = model.goodness_of_fit(flux_opt2[0], alpha = 0.05)
            #
            # Target reactions to estimate 95% CI
            #
            target_reactions = [('reaction', "v4"),('reaction', "v5"),('reaction', "v6"),('reaction', "v7"),('reversible', "FUM")]
            ci_edge = model.generate_ci_templete(targets = target_reactions)
            #
            # Grid search
            #
            model.set_configuration(iteration_max = 10000) # Maximal iternations in optimization
            model.set_configuration(number_of_repeat = 3) #Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
            ci = model.search_ci(ci_edge, flux_opt2[0], method = 'grid')
            #
            # Output of grid search results
            #
            for rid in target_reactions:
                lb = ci['data'][rid]['lower_boundary']
                ub = ci['data'][rid]['upper_boundary']

                print(rid,isotope1,isotope2, "Interval", ub-lb, "Lower bondary:", lb, "Upper boundary:", ub, "RSS_bestfit:",RSS_bestfit[0])






