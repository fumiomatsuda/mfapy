#-------------------------------------------------------------------------------
# Name:        Tutorial_3_mfapy_functions.py
#              Python code to explain various functions of mfapy
#              Toy model used in this code in modified from Antoniewicz et al Metab. Eng. 2007, 9, 68-86.
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
import mfapy
from matplotlib.pyplot import figure, show
import os, sys, time

if __name__ == '__main__':
    #
    # Loading of metabolic model related information "reactions, reversible, metabolites, target_fragments"
    # from given model definition file
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Tutorial_3_mfapy_functions_model.txt", format = "text")
    """
    #-------------------------------------------------------------------------------
    # Name:        Tutorial_3_mfapy_functions_model.txt
    #              Model definition file for Tutorial_3_mfapy_functions.py of mfapy
    #
    # Author:      Fumio_Matsuda
    #
    # Created:     12/06/2018
    # Copyright:   (c) Fumio_Matsuda 2021
    # Licence:     MIT license
    #-------------------------------------------------------------------------------
    #
    #  Any texts after "#" are ignored as comment
    #
    # Definition of metabolic reactions
    # 1. id, text
    #	Please use alphabet (+ number) as ids. Do not use space and special characters to avoid unexpected errors.
    # 2. Stoichiometry for the construicton of stoichiometry matrix, text
    #	All id of metabolites must be defined in //Metabolites
    #	Only "+" and "-->" are accpectable separators.
    #	Here, expression of reaction using coefficient such as  "Fumx --> {2}OA" is allowed.
    #	"nd" means this reaction is an "pseudo" reaction and ignored in the flux calculation of substrate by stoichiometry matrix.
    # 3. Stoichiometry for the construicton of atom mapping, text
    #	Only "+" and "-->" are accpectable separators.
    #	Here, expression of reaction without using coefficient is needed such as  "Fumx --> OA+OA".
    #	"nd" means this reaction is an excretion reaction and ignored in the MDV calculation
    # 4. Atom mapping, text
    #	Only "+" and "-->" are accpectable separators. Please use A-Z to express carbon atom.
    #	"nd" means this reaction is an excretion reaction and ignored in the MDV calculation
    # 5. Links to external ids (Optional)
    #	Links to external ids can be described here
    # 6. Lower boundary of metabolic flux level, float (Optional)
    #	Value must be more than zero
    # 7. Upper boundary of metabolic flux level, float (Optional)
    #	Value must be more than zero
    //Reactions
    v1	AcCoA + OAC --> Cit	AcCoA + OAC --> Cit	AB + CDEF --> FEDBAC	(kegg:R00351)	0.1	300
    v2	Cit --> AKG + CO2ex	Cit --> AKG + CO2ex	ABCDEF --> ABCDE + F	(kegg:R00709)	0.1	300
    v3	AKG --> Glu	AKG --> Glu	ABCDE --> ABCDE 	(kegg:R00243)	0.1	300
    v4	AKG --> Suc + CO2ex	AKG --> Suc + CO2ex	ABCDE --> BCDE + A	(kegg:R01197)	0.1	300
    v5	Suc --> Fum	Suc --> Fum	ABCD --> ABCD	(kegg:R02164)	0.1	300
    v6	Fum --> OAC	Fum --> OAC	ABCD --> ABCD	(kegg:R01082)	0.1	300
    v7	OAC --> Fum	OAC --> Fum	ABCD --> ABCD	(kegg:R01082)	0.1	300
    v8	Asp --> OAC	Asp --> OAC	ABCD --> ABCD	(kegg:R00355)	0.1	300
    v9	Suc --> Sucex	nd	nd	noid	0.1	300
    e1	Glu --> Gluex	nd	nd 	noid	0.1	300
    v10	OAC --> OACx	OAC --> OACx	ABCD --> ABCD	noid	0.1	300
    v11	OACs --> OACx	OACs --> OACx	ABCD --> ABCD	noid	0.001	300
    v12	OACx --> OACe	nd	nd	noid	0.1	300
    v13	OAC+OAC --> Fumx	OAC+OAC --> Fumx	ABCD+EFGH --> EFAB	noid	0.1	300
    v14	OA --> Fume	nd	nd	noid	0.1	300
    v15	OAC+Fum+Suc --> Tx	OAC+Fum+Suc --> Tx	ABCD+EFGH+IJKL --> CDHL	noid	0.1	300
    v16	Tx --> Te	nd	nd	noid	0.1	300
    v17	Fumx --> {2}OA	Fumx --> OA+OA	ABCD --> AB+CD	noid	0.1	300
    #
    # Definition of metabolites
    # 1. id, text
    #	Please use alphabet and numbers. Do not use space, "_", and other special characters to avoid unexpected errors.
    #	Initial character must be alphabet. The id "6PG" was not allowed. Please convert to "m6PG" or other names
    # 2. Carbon number, int
    #   Max number of carbon in "carbonsource" metabolites are 8.
    # 3. symmetry, text
    #	Set "symmetry" metabolites such as fumarete and succinate
    # 4. corbon source, text
    #	Set "carbonsource" metabolites
    # 5. excretion, text
    #	Set "excreted" metabolites
    # 6. Links to external ids (Optional)
    #	Links to external ids can be described here
    # 7. Lower boundary of metabolite pool size, float (Optional)
    #	Value must be more than zero
    #	Required for INST-MFA mode
    # 8. Upper boundary of metabolite pool size, float (Optional)
    #	Value must be more than zero
    #	Required for INST-MFA mode
    #
    //Metabolites
    CO2ex	1	no	no	excreted	(kegg:C00011)	0.0	300
    AcCoA	2	no	carbonsource	no	(kegg:C00024)	0.0	300
    OAC	4	no	no	no	(kegg:C00036)	0.0	300
    OACs	4	no	carbonsource	no	noid	0.0	300
    OACx	4	no	no	no	noid	0.0	300
    OACe	4	no	no	excreted	noid	0.0	300
    Fumx	4	no	no	no	noid	0.0	300
    OA	2	no	no	no	noid	0.0	300
    Fume	4	no	no	excreted	noid	0.0	300
    Tx	4	no	no	no	noid	0.0	300
    Te	4	no	no	excreted	noid	0.0	300
    Cit	6	no	no	no	(kegg:C00158)	0.0	300
    AKG	5	no	no	no	(kegg:C00026)	0.0	300
    Suc	4	symmetry	no	no	(kegg:C00042)	0.0	300
    Fum	4	symmetry	no	no	(kegg:C00122)	0.0	300
    Glu	5	no	no	no	(kegg:C00025)	0.0	300
    Gluex	5	no	no	excreted	noid	0.0	300
    Sucex	4	no	no	excreted	noid	0.0	300
    Asp	4	no	carbonsource	no	(kegg:C00049)	0.0	300
    #
    # Definition of reversible reactions
    # Reversible reaction is used for
    #	(i) output of metabolic flux distribution data
    #	(ii) constrain metabolic flux
    # 	(iii) investigate confidence intervals
    # 	(iii) Composition of two or more reactions such as for the calculation of total ATP producton is supported.
    # 1. Id, text
    #	Please use alphabet and numbers. Please do not use space and special characters to avoid unexpected errors.
    # 2. Dorward reaction, text
    #	Reacton id
    #	Composition of reacton ids such as "v9+e1" is also supported.
    #	Only "+" is accpectable separator.
    # 3. Reverse reaction, text
    #	Reacton id
    #	Set "nothing" if there is no data
    # 4. Links to external ids (Optional)
    #	Links to external ids can be described here
    # 5. Lower boundary of metabolic flux level, float (Optional)
    #	Value must be more than zero
    # 6. Upper boundary of metabolic flux level, float (Optional)
    #	Value must be more than zero
    #
    //Reversible_reactions
    FUM	v6	v7	(kegg:R01082)	-300	300
    Input	v9+e1	nothing	(kegg:R01082)	-300	300
    #
    # Definition of Target fragments
    # Target fragments are metabolite out side of the metabolic model whose MDV are generarted by the calmdv function.
    # 1. id, text
    #	Please use alphabet and numbers. Please do not use space, "_" and other special characters to avoid unexpected errors.
    # 2. type, text
    #	"gcms" is supported.
    # 3. Carbon composition, text
    #	"Glu_12345" indicated that the target compound has five carbons derived from 1-5th carbons of metabolite id "Glu".
    #	In the "gcms"mode, "OAC_12+OAC_12" indicated that the target compound has for carbons derived from 1-2th carbons of metabolite id "OAC".
    # 4. Usage, text
    #	"use" means this target fragment is used in the analysis.
    # 5. Formula, text
    #	Experimental, Chemical formula of target compound
    //Target_fragments
    Glue	gcms	Glu_1:2:3:4:5	use	C5H10N2O3
    AKGe	gcms	AKG_1:2:3:4:5	use	C5H10N2O3
    OACo	gcms	OACx_1:2:3:4	use	C5H10N2O3
    OACi	gcms	OAC_1:2:3:4	use	C5H10N2O3
    OACt	gcms	Fumx_1:2:3:4	use	C5H10N2O3
    OACc	gcms	OAC_1:2+OAC_1:2	use	C5H10N2O3
    Txt	gcms	Tx_1:2:3:4	use	C5H10N2O3
    Txc	gcms	OAC_3:4+Fum_4+Suc_4	use	C5H10N2O3
    Gluf	gcms	Glu_1:2+Glu_3:4:5	no	C5H10N2O3
    AKGms	gcms	AKG_1:2:3:4:5	no
    AKGc	gcms	AKG_1:2+AKG_3:4:5	no	C5H10N2O3
    //End
    """
    #
    # Construction of MetabolicModel instance
    # mode = "debug" is available
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode="normal")
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 1) # Frequency level of callbacks from metabolic model
    model.set_configuration(iteration_max = 10000) # Maximal number of interations (steps) in each fitting task
    model.set_configuration(initial_search_repeats_in_grid_search = 5) # Number of initial metabolic flux disributions generated for each trial in grid search.
    model.set_configuration(initial_search_iteration_max = 1000) # Maximal number of interations (steps) allowed in each task to find feasible initial metabolic flux distribution.
    model.set_configuration(grid_search_iterations = 1) # Number of trials for model fitting at each point in grid search.
    model.set_configuration(number_of_repeat = 3) # Iteration in self.fitting_flux(method = 'deep') [SLSQP => LN_PRAXIS] * n
    model.set_configuration(ncpus = 3) # Number of CPUs for Parallel processing
    #
    # Loading metabolite state from text file
    #
    state_dic = model.load_states("Tutorial_3_mfapy_functions_status.csv", format = 'csv')
    """
    State	Id	type	value	stdev	lb	ub
    reaction	v1	fixed	100	1	0	300
    reaction	v2	free	100	1	0	300
    reaction	v3	free	50	1	0	300
    reaction	v4	free	50	1	0	300
    reaction	v5	free	50	1	0	300
    reaction	v6	free	125	1	0	300
    reaction	v7	free	75	1	0	300
    reaction	v8	free	50	1	0	300
    metabolite	CO2ex	fixed	100	1	0	300
    metabolite	AcCoA	fixed	100	1	0	300
    metabolite	OAC	fixed	100	1	0	300
    metabolite	Cit	fixed	100	1	0	300
    metabolite	AKG	fixed	100	1	0	300
    metabolite	Suc	fixed	100	1	0	300
    metabolite	Fum	fixed	100	1	0	300
    metabolite	Glu	fixed	100	1	0	300
    metabolite	Asp	fixed	10	1	0	300
    reversible	FUM	free	50	1	0	300
    """
    #
    # Batch setting of constraints (type, value, stdev, lower and upper boundaries) of all reactions, metabolites, reversible reactions from given metabolite state "state_dic"
    #
    model.set_constraints_from_state_dict(state_dic)
    #
    # Generation of a metabolic state with randomizaed flux level and metabolite concentration.
    #
    flux_opt, state = model.generate_state()
    flux_opt_original = flux_opt
    #
    # Save metabolic state
    #
    #model.save_states(flux_opt, "Explanation_1_metabolicstate_test.csv", format = 'csv')
    #
    # Manual setting of constraints
    # "fitting (variable and included in RSS calculation using "value" and "stdev")"
    # "fixed (fixed to "value")"
    # "free (variable and ignored for RSS calculation)"
    # "pseudo (variable and used for MDV calc but ignored in the stoichiometry matrix)"
    #
    model.set_constrain('reaction', 'v9', "fitting", 0.0, 1)
    model.set_boundary('reaction', 'v9', 0.0, 30.0)
    #
    # self.update() is required after modification of constraint "type"
    #
    model.update()
    #
    # Generation of instances of CarbonSource class
    #
    cs = model.generate_carbon_source_template()
    cs2 = model.generate_carbon_source_template()
    #
    # Isotope labelling information of carbon source is set by ratio of all isotopomer (#00, #01, #10, #11)
    #
    cs.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.25, 0.25])
    #
    # Isotope labelling information of carbon source is set by specific isotopomers
    #
    cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = "yes")
    #
    # Isotope labelling information of carbon sources are set from file
    #
    """
    Name	Isotopomer	Ratio
    Asp	#0000	0.5
    Asp	#1111	0.5
    AcCoA	#00	0.5
    AcCoA	#11	0.5
    """
    cs2.set_carbonsources('Tutorial_3_mfapy_functions_carbonsource2.txt', correction = "no")
    #
    # Artificial parallel labeling experiment to test "calcrss" function
    #
    #
    # Generation of artificial MDV data from flux_opt and two carbon source settings.
    #
    mdv1 = model.generate_mdv(flux_opt, cs)
    mdv2 = model.generate_mdv(flux_opt, cs2)
    #
    # Standard deviation of all mdv measurements are set at 0.01
    #
    mdv1.set_std(0.01)
    mdv2.set_std(0.01)
    #
    # Two labeling experiments are set to model.
    #
    model.set_experiment('ex1', mdv1, cs)
    model.set_experiment('ex2', mdv2, cs2)
    #
    # Test "calcrss" function
    #
    print("Calculatino of RSS is OK if RSS is zero, RSS:", model.calc_rss(flux_opt))
    #
    # Setting v8 as fitting
    #
    model.set_constrain('reaction','v8','fitting', value = 49, stdev = 0.1)
    #
    # Update model
    #
    model.update()
    print("Calculatino of RSS including fitting flux is OK if RSS is close to 100, RSS:", model.calc_rss(flux_opt))
    #
    # Removal of all experiments
    #
    model.clear_experiment()
    #
    # Generate two artificial MDV data from metabolic state and carbon source settings.
    #
    mdv1 = model.generate_mdv(flux_opt, cs)
    mdv2 = model.generate_mdv(flux_opt, cs2)
    #
    # Addition of gaussian noises for the artificial MDV data. Sample (n=3) gaussian noise with STD = 0.01, then calc their average and stdev
    #
    mdv1.add_gaussian_noise(0.01)
    mdv2.add_gaussian_noise(0.01)
    #
    # Standard deviation of all mdv measurements are set at 0.01
    #
    mdv1.set_std(0.01)
    mdv2.set_std(0.01)
    #
    # Ignore MDV data whose relative intensity is outside between 0.01-1.0
    #
    mdv1.set_mdvs_for_comparison(0.01, 1.0)
    mdv2.set_mdvs_for_comparison(0.01, 1.0)
    #
    # Saving MDV data in text file
    #
    mdv1.save("Tutorial_3_mfapy_functions_mdv1.txt", format = "text")
    #
    # Load MDV data from file
    #
    mdv3 = model.load_mdv_data("Tutorial_3_mfapy_functions_mdv1.txt", format = "text", output = "normal")
    #
    # Addition of two labelling experiments (Parallel labeling experiments)
    #
    model.set_experiment('ex1', mdv3, cs)
    model.set_experiment('ex2', mdv2, cs2)
    #
    # Setting v8 and FUM as free
    #
    model.set_constrain('reaction','v8','free')
    model.set_constrain('reversible','FUM','free')
    #
    # Update model
    #
    model.update()
    #
    # "template" metabolic state is stored in 'results'
    #
    results = [('template', flux_opt)]
    #
    # Generation of initial metabolic states
    # Initial metabolic states are generate by 50 trials (sometimes failed)
    # then return 1 initial metabolic state with smallest RSSs.
    # In this case 1 dictionary of metabolic state is generated.
    #
    # method = "parallel" indicates parallel processing using joblib.
    # Please change number of local CPUs for joblib my model.set_configuration(ncpus = 3)
    #
    state, flux_initial1 = model.generate_initial_states(50, 1, method = "parallel")
    results.append(("Initial", flux_initial1))
    #
    # Calculation of residual sum of square
    #
    print("Initial metabolic state was successfully obtained.")
    print("Calculation of RSS of single flux data", model.calc_rss(flux_opt))
    print("Calculation of RSS of multiple flux data", model.calc_rss([flux_opt, flux_initial1]))
    #
    # Find best fitted metabolic state by SLSQP (scipy)
    #
    start = time.time()
    state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = 'SLSQP', flux = flux_initial1)
    print("Fitting_flux(SLSQP scipy integrate): State", state, "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit, time.time()-start))
    results.append(("SLSQP", flux_opt_slsqp))
    #
    # Calculation of goodness of fit
    #
    pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)
    print("Calculation of goodness of fit. RSS", RSS_bestfit,"Results: p-value:",pvalue, "thres", rss_thres)
    #
    # Calculation of threshold RSS for confidence interval search
    #
    thres, number_of_measurements, degree_of_freedom  =\
    model.get_thres_confidence_interval(flux_opt_slsqp, alpha = 0.05)
    print("Calculation of threshold RSS for confidence interval search. Threshold:",thres," Number of measurements:", number_of_measurements," Degree of freedom", degree_of_freedom)
    #
    # Additional fitting by local and global optimizers
    #
    methods = ["COBYLA", "LN_COBYLA", "LN_BOBYQA", "LN_NEWUOA", "LN_PRAXIS", "LN_SBPLX", "LN_NELDERMEAD", "GN_CRS2_LM", "deep"]
    for method in methods:
        start = time.time()
        state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = method, flux = flux_initial1)
        pvalue, rss_thres = model.goodness_of_fit(flux_opt_slsqp, alpha = 0.05)

        print(method, ": State", state, "RSS:{0:>8.2f} Time:{1:>8.2f} Threshold:{2:>8.2f} pvalue:{3:>8.7f}".format(RSS_bestfit, time.time()-start, rss_thres, pvalue))
        results.append((method, flux_opt_slsqp))

    #
    # Initial metabolic states were generate by 50 trials (sometimes failed), then return 4 initial metabolic state with smallest RSSs.
    # In this case, a list of 4 dictionarirs of metabolic state is generated.
    #
    state, flux_initial2 = model.generate_initial_states(50, 4, method ="parallel")
    start = time.time()
    #
    # Finding best fitted metabolic state by SLSQP using parallel mode
    # Parallel mode was automatically employed when a list of dictionarirs of metabolic state is given
    #
    state, RSS_bestfit, flux_opt_slsqp2 = model.fitting_flux(method = 'SLSQP', flux = flux_initial2)
    print("Fitting_flux(SLSQP scipy integrate, parallel python): State", "RSS:{0:>8.2f} Time:{1:>8.2f}".format(RSS_bestfit[0], time.time()-start))
    #
    # Output fitted metabolic flux into CSV file
    #
    model.show_results(results, filename = "Tutorial_3_mfapy_functions_output.csv",  checkrss = "on", format = "csv")
    #
    # Output fitted metabolic flux in console
    #
    model.show_results(results,  checkrss = "on")
    #
    # Output flux balance of given metabolite
    #
    model.show_flux_balance(results, "Fum")
    model.show_flux_balance(results, "Fum", filename = "Tutorial_3_mfapy_functions_fluxbalance.csv", format = "csv")# to file
    #
    # Comparision between MDVs
    # Genaration of MDVs for given flux states
    #
    mdv_initial1 = model.generate_mdv(flux_initial1, cs)
    mdv_flux_opt_slsqp = model.generate_mdv(flux_opt_slsqp, cs)
    #
    # Array of MDV vector is generated using "generate_observed_mdc" method
    #
    id_array, ratio_array, std_array, use_array, observed_fragments, data = mdv1.generate_observed_mdv()
    id_array, ratio_array_initial, std_array, use_array, observed_fragments, data = mdv_initial1.generate_observed_mdv()
    id_array, ratio_array_flux_opt_slsqp, std_array, use_array, observed_fragments, data = mdv_flux_opt_slsqp.generate_observed_mdv()
    #
    # Scatter plot using matplotlib
    #
    fig = figure()
    axL = fig.add_subplot(121)
    axL.scatter(ratio_array_initial, ratio_array, c = 'r')
    axL.scatter(ratio_array_flux_opt_slsqp, ratio_array,  c = 'b')
    #
    # Target reactions to estimate 95% CI and generation of template disctionary to store CI results
    #
    target_reactions = [('reaction', "v6"),('reversible', "FUM")]
    ci_edge = model.generate_ci_template(targets = target_reactions)
    #
    # Grid search related configurations
    #
    print("Test: Searching confidencial interval by grid search method.")
    #model.set_configuration(callbacklevel = 0)
    model.set_configuration(initial_search_repeats_in_grid_search = 5) # Number of repeat for finding initial flux in Grid search mode such as model.generate_initial_states(5, 1)
    model.set_configuration(grid_search_iterations = 1) # Fitting interations at each grid.
    #
    # Grid search
    #
    ci = model.search_ci(ci_edge, flux_opt_slsqp, method = 'grid')
    #
    # Visualization of grid search results
    #
    for rid in ci['data'].keys():
        if ci['data'][rid]['use'] != 'on':continue
        print(rid, "Lower bondary:",ci['data'][rid]['lower_boundary'], "Upper boundary:", ci['data'][rid]['upper_boundary'])

    axR = fig.add_subplot(122)
    axR.scatter(ci['data'][('reversible', "FUM")]['flux_data'], ci['data'][('reversible', "FUM")]['rss_data'], c = 'r')
    axR.scatter(ci['data'][("reaction", "v6")]['flux_data'], ci['data'][("reaction", "v6")]['rss_data'], c = 'b')
    show()

