#-------------------------------------------------------------------------------
# Name:        optimize.py
# Purpose:     low level optimizer functions used in mfapy. These functions were separated from model instance for the parallel execution.
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------

import numpy as numpy
import scipy as scipy
import nlopt as nlopt
#from numba import jit

def initializing_Rm_fitting(numbers, vectors, matrixinv, template, initial_search_iteration_max, method = "fitting"):
    """
    Generate randomized initial flux dixtribution using scipy.optimize.minimize SLSQP

    Parameters
    ----------
    numbers: "model.numbers" including various numbers of the model.
    vectors: "model.vector" including various vectors of the model.
    matrixinv: "model.matrixinv" is a inversed matrix for the flux calculation.
    template: Dictionary of metabolic state. When template is available, metabolic state most similar to the template is generated. The function is used in the grid search.
    initial_search_iteration_max: "configure["initial_search_iteration_max"]". Number of
    method: only "fitting" is available.

    Returns
    -------
    tmp_r: List of metabolic state (tmp_r = numpy.dot(matrixinv, Rm_temp)
    Rm_temp: List of metabolic state vector
    Rm_ind: List of independent flux
    state: State of finishing condition "Failed"/"Determined"

    Examples
    --------
    >>> tmp_r, Rm_temp, Rm_ind, state = optimize.initializing_Rm_fitting(numbers, vectors, matrixinv, template ,initial_search_iteration_max)


    See Also
    --------
    calc_protrude_scipy

    """

    # number of independent flux
    independent_number = numbers['independent_number']
    total_number = numbers['total_number']
    #
    # Mas=number of MKL thread control
    #
    try:
        import mkl
        mkl.set_num_threads(1)
    except:
        print("mkl-service is not installed this python!")
    # zero independent flux
    Rm_ind = list(numpy.zeros(independent_number))
    #boundaries
    lb = list(vectors["lb"])
    ub = list(vectors["ub"])
    independent_lb = list(vectors["independent_lb"])
    independent_ub = list(vectors["independent_ub"])
    for j in range(3):
        #ランダムな反応の最低値をランダムな値にする。
        lb_modified = list(lb)
        ub_modified = list(ub)

        #適当に初期値を発生させる
        for i in range(len(Rm_ind)):
            Rm_ind[i] = (independent_ub[i] - independent_lb[i]) * numpy.random.rand() + independent_lb[i]
        # Instantiate Optimization Problem

        parameters = {}
        parameters['stoichiometric_num'] = numbers['independent_start']
        parameters['reaction_num']=numbers['independent_end']
        parameters['matrixinv']=matrixinv
        parameters['Rm_initial']=numpy.array(list(vectors["Rm_initial"]))
        parameters['lb'] = lb_modified
        parameters['ub'] = ub_modified
        parameters['template'] = template

        #
        #    Scipy
        #
        res = scipy.optimize.minimize(calc_protrude_scipy, Rm_ind, method='SLSQP', args = (parameters,))
        result_ind = res.x
        """

        #
        # nlopt
        #
        opt = nlopt.opt(nlopt.LN_COBYLA, independent_number)
        opt.set_lower_bounds(independent_lb)
        opt.set_upper_bounds(independent_ub)

        opt.set_min_objective(lambda x,grad: calc_protrude_nlopt(x,grad,parameters))
        opt.set_xtol_abs(0.0001)
        opt.set_maxeval(initial_search_iteration_max)
        result_ind = opt.optimize(Rm_ind)
        minf = opt.last_optimum_value()

        """
        result_Rm = numpy.array(list(vectors["Rm_initial"]))
        result_Rm[numbers['independent_start']: numbers['independent_end']] = result_ind[:]
        tmp_r = numpy.dot(matrixinv, result_Rm)
        check = 0;
        for i in range(len(tmp_r)):
            if tmp_r[i] < lb[i] - 0.0001:
                check = check + 1
            if tmp_r[i] > ub[i] + 0.0001:
                check = check + 1
        if check == 0:
            return(tmp_r, result_Rm, result_ind, "Determined")
        #print "mistake", tmp_r
    return(tmp_r, result_Rm, result_ind, "Failed")

def calc_protrude_scipy(independent_flux, *args):
    """
    Objective function for initializing_Rm_fitting (SLSQP)
    Calc penalty score of metabolic state out side of the feasible space.

    Parameters
    ----------
    independent_flux: vector of independent flux
    *args: list of parameters.

    Returns
    -------
    f: Penalty score


    See Also
    --------
    initializing_Rm_fitting


    """
    kwargs = args[0]
    Rm_initial = kwargs['Rm_initial']
    stoichiometric_num = kwargs['stoichiometric_num']
    reaction_num = kwargs['reaction_num']
    matrixinv = kwargs['matrixinv']
    lb = kwargs['lb']
    ub = kwargs['ub']
    template = kwargs['template']

    Rm = numpy.array(list(Rm_initial))
    Rm[stoichiometric_num: reaction_num] = independent_flux[:]
    #
    tmp_r = numpy.dot(matrixinv, Rm)
    #
    f = 0.0
    g = []
    if len(template) > 0:
        #
        # if templete flux is available
        #
        for i, flux in enumerate(tmp_r):
            #Between lower and upper boundary
            g.append(flux - ub[i])
            g.append(lb[i]- flux)
            f = f + abs(flux - template[i])

    else:
        #
        # to generate random flux
        #
        for i, flux in enumerate(tmp_r):
            #Between lower and upper boundary
            g.append(flux - ub[i])
            g.append(lb[i]- flux)

            if flux > ub[i]:
                f = f + (flux - ub[i])
            elif flux < lb[i]:
                f = f + (lb[i] - flux)
    fail = 0
    #print(f)

    return f

def calc_protrude_nlopt(independent_flux, grad, kwargs):
    """
    Objective function for initializing_Rm_fitting (nlpot)
    Calc penalty score of metabolic state out side of the feasible space.

    Parameters
    ----------
    independent_flux: vector of independent flux
    grad: not used
    *args: list of parameters.

    Returns
    -------
    f: Penalty score


    See Also
    --------
    initializing_Rm_fitting


    """
    Rm_initial = kwargs['Rm_initial']
    stoichiometric_num = kwargs['stoichiometric_num']
    reaction_num = kwargs['reaction_num']
    matrixinv = kwargs['matrixinv']
    lb = kwargs['lb']
    ub = kwargs['ub']
    template = kwargs['template']

    Rm = numpy.array(list(Rm_initial))
    Rm[stoichiometric_num: reaction_num] = independent_flux[:]
    #
    tmp_r = numpy.dot(matrixinv, Rm)
    #
    f = 0.0
    g = []
    if len(template) > 0:
        #
        # if templete flux is available
        #
        for i, flux in enumerate(tmp_r):
            #Between lower and upper boundary
            g.append(flux - ub[i])
            g.append(lb[i]- flux)
            f = f + abs(flux - template[i])

    else:
        #
        # to generate random flux
        #
        for i, flux in enumerate(tmp_r):
            #Between lower and upper boundary
            g.append(flux - ub[i])
            g.append(lb[i]- flux)

            if flux > ub[i]:
                f = f + (flux - ub[i])
            elif flux < lb[i]:
                f = f + (lb[i] - flux)
    fail = 0
    #print(f, grad)

    return f

def calc_MDV_from_flux(tmp_r, target_fragments, mdv_carbon_sources, func, timepoint = [], y0temp = []):
    """
    Low level function to calculate mdv vector and mdv hash from metabolic flux and carbon
    source MDV using calmdv. This funcition is called from mfapy.metabolicmodel.show_results.

    Parameters
    ----------
    tmp_r: list of metabolix state
    target_fragments: list of targed mdvs for MDV calclation, model.target_fragments.keys()
    mdv_carbon_sources: dict of mdv_carbon_sources in model.experiments[ex_id]['mdv_carbon_sources']
    func: Dict of functions for MDV calclation in model.func
    timepoint: For INST mode only. timepoints for MDV comparison in model.experiments[ex_id]['timepoint']
    y0temp: Start IDV state for INST mode

    Returns
    -------
    f: Penalty score

    Example
    -------
    mdv_exp, mdv_hash = optimize.calc_MDV_from_flux(tmp_r, target_fragments_temp, mdv_carbon_sources_temp, self.func)


    See Also
    --------
    mfapy.metabolicmodel.show_results


    """
    if len(timepoint)==0:
        mdv, mdv_hash = func["calmdv"](tmp_r, sorted(target_fragments), mdv_carbon_sources)
    else:
        mdv, mdv_hash = func["diffmdv"](tmp_r, [], timepoint, sorted(target_fragments), mdv_carbon_sources, y0temp)
    return mdv, mdv_hash





def fit_r_mdv_scipy(configure, experiments, numbers, vectors, matrixinv, func, flux,  method = "SLSQP"):
    """
    Low level function for model fitting using scipy.optimize.minimize

    Parameters
    ----------
    configures: "model.configures" including various configulatoins of the model.
    experiments: "model.experiments" including experiments defined in the model.
    numbers: "model.numbers" including various numbers of the model.
    vectors: "model.vector" including various vectors of the model.
    matrixinv: "model.matrixinv" is a inversed matrix for the flux calculation.
    func: Dict of functions for MDV calclation in model.func
    flux: Dictionary of initial metabolic state.
    method: "SLSQP" and "COBYLA" are available

    Returns
    -------
    state: finishing condition
    kai: Residual sum of square of fitted metabolic state
    opt_flux: list of fitted metabolix state
    Rm_ind_sol: list of fitted independent flux

    Example
    -------
    state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configure, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "SLSQP")


    See Also
    --------
    calc_MDV_residue_scipy

    """

    if isinstance(func, dict):
        calmdv = func["calmdv"]
        diffmdv = func["diffmdv"]
    else:
        locals_dic = locals()
        exec(func, globals(), locals_dic)
        calmdv = locals_dic["calmdv"]
        diffmdv = locals_dic["diffmdv"]
    #
    # Mas=number of MKL thread control
    #
    try:
        import mkl
        mkl.set_num_threads(1)
    except:
        print("mkl-service is not installed this python!")
    #
    #
    #Set max number of iteration in pyOpt
    if 'iteration_max' in configure:
        iteration_max = configure['iteration_max']
    else:
        iteration_max = 1000

    #Set callbacklevel
    if 'callbacklevel' in configure:
        callbacklevel = configure['callbacklevel']
    else:
        callbacklevel = 0
    # number of independent flux
    independent_number = numbers['independent_number']
    ind_start = numbers['independent_start']
    ind_end = numbers['independent_end']

    total_number = numbers['total_number']
    # zero independent flux
    if isinstance(flux, dict):
        Rm_ind = [flux[group][id]["value"] for (group, id) in vectors['independent_flux']]
    else:
        Rm_ind = [flux[i] for i in vectors['independent_flux_position']]
    #boundaries
    lb = list(vectors["lb"])
    ub = list(vectors["ub"])
    independent_lb = list(vectors["independent_lb"])
    independent_ub = list(vectors["independent_ub"])
    #
    # Generate MDV vector of all defined experiments
    #
    mdv_exp_original = list(vectors["value"])
    mdv_std_original = list(vectors["stdev"])
    mdv_use = list(vectors["use"])
    for experiment in sorted(experiments.keys()):
        mdv_exp_original.extend(experiments[experiment]['mdv_exp_original'])
        mdv_std_original.extend(experiments[experiment]['mdv_std_original'])
        mdv_use.extend(experiments[experiment]['mdv_use'])
    mdv_exp = numpy.array([y for x, y in enumerate(mdv_exp_original) if mdv_use[x] != 0])
    spectrum_std = numpy.array([y for x, y in enumerate(mdv_std_original) if mdv_use[x] != 0])
    #
    # Covariance matrix
    #
    covinv = numpy.zeros((len(spectrum_std),len(spectrum_std)))
    for i, std in enumerate(spectrum_std):
        covinv[i,i] = 1.0/(std**2)

    state = {'text':"Function was called", 'value': 7}
    for k in range(2):
        #try:
        ##################################################################
        if (callbacklevel >= 2):
            print(k, "Fitting Start")
        parameters = {}
        parameters['stoichiometric_num'] = ind_start
        parameters['reaction_num']=ind_end
        parameters['matrixinv']=matrixinv
        parameters['experiments']=experiments
        parameters['mdv_exp'] = mdv_exp
        parameters['mdv_use'] = mdv_use
        parameters['covinv']= covinv
        parameters['Rm_initial']=numpy.array(list(vectors["Rm_initial"]))
        parameters['lb'] = lb
        parameters['ub'] = ub
        parameters['calmdv'] = calmdv
        parameters['diffmdv'] = diffmdv
        options={'ftol': 0.000000001, 'maxiter': iteration_max}
        method_scipy = "SLSQP"
        bounds = []
        for i in range(independent_number):
            bounds.append((independent_lb[i],independent_ub[i]))
            #print(independent_lb[i],Rm_ind[i], independent_ub[i])
        if method == "SLSQP":
            options={'ftol': 0.000000001, 'maxiter': iteration_max}
            method_scipy = "SLSQP"
            res = scipy.optimize.minimize(calc_MDV_residue_scipy, Rm_ind, bounds = bounds, options = options, method=method_scipy, args = (parameters,))

        elif method == "COBYLA":
            options={'tol': 0.000000001, 'maxiter': iteration_max}
            method_scipy = "COBYLA"
            res = scipy.optimize.minimize(calc_MDV_residue_scipy, Rm_ind, options = options, method=method_scipy, args = (parameters,))
        else:
            options={'ftol': 0.000000001, 'maxiter': iteration_max}
            method_scipy = "SLSQP"
            res = scipy.optimize.minimize(calc_MDV_residue_scipy, Rm_ind, bounds = bounds, options = options, method=method_scipy, args = (parameters,))

        #Optimized flux distribution
        result_ind = res.x
        #RSS
        kai = res.fun
        #State of optimizer
        state = res.message

        if (callbacklevel >= 2):
            print(k+1, "th trial. ",state ,"RSS = ", kai)

        #Optimized flux distribution
        Rm_opt = numpy.array(list(vectors["Rm_initial"]))
        result_Rm = numpy.array(list(vectors["Rm_initial"]))
        result_Rm[numbers['independent_start']: numbers['independent_end']] = result_ind[:]

        opt_flux = numpy.dot(matrixinv, numpy.array(result_Rm))

        return(state, kai, opt_flux, result_ind)


    return(state,-1,[],[])

def fit_r_mdv_nlopt(configure, experiments, numbers, vectors, matrixinv, func, flux,  method = "LN_PRAXIS"):
    """
    Low level function for model fitting using nlopt.opt

    Parameters
    ----------
    configures: "model.configures" including various configulatoins of the model.
    experiments: "model.experiments" including experiments defined in the model.
    numbers: "model.numbers" including various numbers of the model.
    vectors: "model.vector" including various vectors of the model.
    matrixinv: "model.matrixinv" is a inversed matrix for the flux calculation.
    func: Dict of functions for MDV calclation in model.func
    flux: Dictionary of initial metabolic state.
    method: "LN_COBYLA", "LN_BOBYQA", "LN_NEWUOA", "LN_PRAXIS", "LN_SBPLX", "LN_NELDERMEAD", "GN_DIRECT_L", "GN_CRS2_LM","GN_ESCH"

    Returns
    -------
    state: finishing condition
    kai: Residual sum of square of fitted metabolic state
    opt_flux: list of fitted metabolix state
    Rm_ind_sol: list of fitted independent flux

    Example
    -------
    state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configure, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_PRAXIS")


    See Also
    --------
    calc_MDV_residue_nlopt

    """

    if isinstance(func, dict):
        calmdv = func["calmdv"]
        diffmdv = func["diffmdv"]
    else:
        locals_dic = locals()
        exec(func, globals(), locals_dic)
        calmdv = locals_dic["calmdv"]
        diffmdv = locals_dic["diffmdv"]
    #
    # Mas=number of MKL thread control
    #
    try:
        import mkl
        mkl.set_num_threads(1)
    except:
        print("mkl-service is not installed this python!")
    #
    #
    #Set max number of iteration in pyOpt
    if 'iteration_max' in configure:
        iteration_max = configure['iteration_max']
    else:
        iteration_max = 1000

    #Set callbacklevel
    if 'callbacklevel' in configure:
        callbacklevel = configure['callbacklevel']
    else:
        callbacklevel = 0
    # number of independent flux
    independent_number = numbers['independent_number']
    ind_start = numbers['independent_start']
    ind_end = numbers['independent_end']

    total_number = numbers['total_number']
    # zero independent flux
    if isinstance(flux, dict):
        Rm_ind = [flux[group][id]["value"] for (group, id) in vectors['independent_flux']]
    else:
        Rm_ind = [flux[i] for i in vectors['independent_flux_position']]
    #boundaries
    lb = list(vectors["lb"])
    ub = list(vectors["ub"])
    independent_lb = list(vectors["independent_lb"])
    independent_ub = list(vectors["independent_ub"])
    #
    # Generate MDV vector of all defined experiments
    #
    mdv_exp_original = list(vectors["value"])
    mdv_std_original = list(vectors["stdev"])
    mdv_use = list(vectors["use"])
    for experiment in sorted(experiments.keys()):
        mdv_exp_original.extend(experiments[experiment]['mdv_exp_original'])
        mdv_std_original.extend(experiments[experiment]['mdv_std_original'])
        mdv_use.extend(experiments[experiment]['mdv_use'])
    mdv_exp = numpy.array([y for x, y in enumerate(mdv_exp_original) if mdv_use[x] != 0])
    spectrum_std = numpy.array([y for x, y in enumerate(mdv_std_original) if mdv_use[x] != 0])
    #
    # Covariance matrix
    #
    covinv = numpy.zeros((len(spectrum_std),len(spectrum_std)))
    for i, std in enumerate(spectrum_std):
        covinv[i,i] = 1.0/(std**2)

    state = {'text':"Function was called", 'value': 7}
    for k in range(2):
        #try:
        ##################################################################
        if (callbacklevel >= 1):
            print(k, "Fitting Start")
        parameters = {}
        parameters['stoichiometric_num'] = ind_start
        parameters['reaction_num']=ind_end
        parameters['matrixinv']=matrixinv
        parameters['experiments']=experiments
        parameters['mdv_exp'] = mdv_exp
        parameters['mdv_use'] = mdv_use
        parameters['covinv']= covinv
        parameters['Rm_initial']=numpy.array(list(vectors["Rm_initial"]))
        parameters['lb'] = lb
        parameters['ub'] = ub
        parameters['calmdv'] = calmdv
        parameters['diffmdv'] = diffmdv

        #
        # nlopt
        #
        if method == "LN_COBYLA":
            opt = nlopt.opt(nlopt.LN_COBYLA, independent_number)
        elif method == "LN_BOBYQA":
            opt = nlopt.opt(nlopt.LN_BOBYQA, independent_number)
        elif method == "LN_NEWUOA":
            opt = nlopt.opt(nlopt.LN_NEWUOA, independent_number)
        elif method == "LN_PRAXIS":
            opt = nlopt.opt(nlopt.LN_PRAXIS, independent_number)
        elif method == "LN_SBPLX":
            opt = nlopt.opt(nlopt.LN_SBPLX, independent_number)
        elif method == "LN_NELDERMEAD":
            opt = nlopt.opt(nlopt.LN_NELDERMEAD, independent_number)
        elif method == "GN_DIRECT_L":
            opt = nlopt.opt(nlopt.GN_DIRECT_L, independent_number)
        elif method == "GN_CRS2_LM":
            opt = nlopt.opt(nlopt.GN_CRS2_LM, independent_number)
        elif method == "GN_ESCH":
            opt = nlopt.opt(nlopt.GN_ESCH, independent_number)
        else:
            opt = nlopt.opt(nlopt.LN_COBYLA, independent_number)
        opt.set_xtol_abs(0.000001)
        opt.set_maxeval(iteration_max)
        opt.set_lower_bounds(independent_lb)
        opt.set_upper_bounds(independent_ub)

        opt.set_min_objective(lambda x,grad: calc_MDV_residue_nlopt(x,grad,parameters))

        #
        # Optimizaton
        #
        result_ind = opt.optimize(Rm_ind)
        kai = opt.last_optimum_value()


        if (callbacklevel >= 1):
            print(k+1, "th trial. Fitting was successfully finished. RSS = ", kai)

        #Optimized flux distribution
        Rm_opt = numpy.array(list(vectors["Rm_initial"]))
        result_Rm = numpy.array(list(vectors["Rm_initial"]))
        result_Rm[numbers['independent_start']: numbers['independent_end']] = result_ind[:]

        opt_flux = numpy.dot(matrixinv, numpy.array(result_Rm))

        return(state, kai, opt_flux, result_ind)


    return(state,-1,[],[])

def fit_r_mdv_deep(configure, experiments, numbers, vectors, matrixinv, func, flux):
    """
    Low level function for model fitting by interations of fittings.
    2n th iteration:  SLSQP
    2n + 1 th iteration:  LN_PRAXIS
    This combination is empirically best

    Parameters
    ----------
    configures: "model.configures" including various configulatoins of the model.
    experiments: "model.experiments" including experiments defined in the model.
    numbers: "model.numbers" including various numbers of the model.
    vectors: "model.vector" including various vectors of the model.
    matrixinv: "model.matrixinv" is a inversed matrix for the flux calculation.
    func: Dict of functions for MDV calclation in model.func
    flux: Dictionary of initial metabolic state.

    Returns
    -------
    state: finishing condition
    kai: Residual sum of square of fitted metabolic state
    opt_flux: list of fitted metabolix state
    Rm_ind_sol: list of fitted independent flux

    Example
    -------
    state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_deep(configure, self.experiments, numbers, vectors, self.matrixinv, self.func, flux)


    See Also
    --------
    optimize.fit_r_nlopt
    optimize.fit_r_scipy


    """


    #Set max number of repeat
    if 'number_of_repeat' in configure:
        number_of_repeat = configure['number_of_repeat']
    else:
        number_of_repeat = 3

    #Set max number of iteration in pyOpt
    if 'iteration_max' in configure:
        iteration_max = configure['iteration_max']
    else:
        iteration_max = 1000

    #Set callbacklevel
    if 'callbacklevel' in configure:
        callbacklevel = configure['callbacklevel']
    else:
        callbacklevel = 0
    if (callbacklevel >= 2):
        print("##Start GN_CRS2_LM method#############################################################################")
    state, kai, flux, Rm_ind_sol = fit_r_mdv_nlopt(configure, experiments, numbers, vectors, matrixinv, func, flux, method = "GN_CRS2_LM")

    for k in range (number_of_repeat):
        if (callbacklevel >= 2):
            print("##Start SLSQP method#############################################################################")
        state, kai, flux, Rm_ind_sol = fit_r_mdv_scipy(configure, experiments, numbers, vectors, matrixinv, func, flux, method = "SLSQP")
        if (callbacklevel >= 2):
            print("##Start LN_PRAXIS method#########################################################################")
        state, kai, flux, Rm_ind_sol = fit_r_mdv_nlopt(configure, experiments, numbers, vectors, matrixinv, func, flux, method = "LN_PRAXIS")



    return(state, kai, flux, Rm_ind_sol)





def calc_MDV_residue_scipy(x, *args):
    """
    Low level function for residual sum of square calculation for model fitting using scipy.optimize.minimize

    Parameters
    ----------
    x: list. vector of independent flux
    *args: list of parameters.

    Returns
    -------
    f: RSS + Penalty score (When out side of the lower and upper boundaries)


    See Also
    --------
    fit_r_mdv_scipy

    """
    kwargs = args[0]
    Rm_initial = kwargs['Rm_initial']
    stoichiometric_num = kwargs['stoichiometric_num']
    reaction_num = kwargs['reaction_num']
    reac_met_num = kwargs['reaction_num']
    matrixinv = kwargs['matrixinv']
    experiments = kwargs['experiments']
    mdv_exp = numpy.array(kwargs['mdv_exp'])
    mdv_use = kwargs['mdv_use']
    covinv = kwargs['covinv']
    lb = kwargs['lb']
    ub = kwargs['ub']
    calmdv = kwargs['calmdv']
    diffmdv = kwargs['diffmdv']

    Rm = numpy.array(list(Rm_initial))
    Rm[stoichiometric_num: reaction_num] = list(x)
    tmp_r = numpy.dot(matrixinv, Rm)

    g = numpy.hstack((numpy.array(lb) - tmp_r, tmp_r - numpy.array(ub)))
    sum = 0.0
    for i in g:
        if i > 0:
            sum = sum + i * 100000
            #print(i)
    fail = 0
    #Determination of MDV
    mdv_original = list(tmp_r)
    for experiment in sorted(experiments.keys()):
        target_emu_list = experiments[experiment]['target_emu_list']
        mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']

        #
        #
        #
        if experiments[experiment]['mode'] == "ST":
            mdv_original_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)
        elif experiments[experiment]['mode'] == "INST":
            y0temp = experiments[experiment]['y0']
            timepoints = experiments[experiment]['timepoint']
            mdv_original_temp, mdv_hash = diffmdv(list(tmp_r), [], timepoints, target_emu_list, mdv_carbon_sources, y0temp)
        mdv_original.extend(mdv_original_temp)

    mdv = numpy.array([y for x, y in enumerate(mdv_original) if mdv_use[x] != 0])
    res = mdv_exp - mdv
    f = numpy.dot(res, numpy.dot(covinv, res))
    if experiments[experiment]['mode'] == "INST":
        print(f)
    return f+sum

def calc_MDV_residue_nlopt(x, grad, kwargs):
    """
    Low level function for residual sum of square calculation for model fitting using nlopt.nlopt

    Parameters
    ----------
    x: list. vector of independent flux
    *args: list of parameters.

    Returns
    -------
    f: RSS + Penalty score (When out side of the lower and upper boundaries)


    See Also
    --------
    fit_r_mdv_scipy

    """
    Rm_initial = kwargs['Rm_initial']
    stoichiometric_num = kwargs['stoichiometric_num']
    reaction_num = kwargs['reaction_num']
    reac_met_num = kwargs['reaction_num']
    matrixinv = kwargs['matrixinv']
    experiments = kwargs['experiments']
    mdv_exp = numpy.array(kwargs['mdv_exp'])
    mdv_use = kwargs['mdv_use']
    covinv = kwargs['covinv']
    lb = kwargs['lb']
    ub = kwargs['ub']
    calmdv = kwargs['calmdv']
    diffmdv = kwargs['diffmdv']

    Rm = numpy.array(list(Rm_initial))
    Rm[stoichiometric_num: reaction_num] = list(x)
    tmp_r = numpy.dot(matrixinv, Rm)

    g = numpy.hstack((numpy.array(lb) - tmp_r, tmp_r - numpy.array(ub)))
    sum = 0.0
    for i in g:
        if i > 0:
            sum = sum + i * 100000
            #print(i)

    fail = 0
    #Determination of MDV
    mdv_original = list(tmp_r)
    for experiment in sorted(experiments.keys()):
        target_emu_list = experiments[experiment]['target_emu_list']
        mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
        #
        #
        #
        if experiments[experiment]['mode'] == "ST":
            mdv_original_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)
        elif experiments[experiment]['mode'] == "INST":
            y0temp = experiments[experiment]['y0']
            timepoints = experiments[experiment]['timepoint']
            mdv_original_temp, mdv_hash = diffmdv(list(tmp_r), [], timepoints, target_emu_list, mdv_carbon_sources, y0temp)
        mdv_original.extend(mdv_original_temp)

    mdv = numpy.array([y for x, y in enumerate(mdv_original) if mdv_use[x] != 0])
    res = mdv_exp - mdv
    f = numpy.dot(res, numpy.dot(covinv, res))
    if experiments[experiment]['mode'] == "INST":
        print(f)

    return f+sum


def calc_MDV_residue(x, *args, **kwargs):
    """
    Low level function for residual sum of square calculation from mfapy.metabolicmodel.MetaboliModel.calc_rss

    Parameters
    ----------
    x: list. vector of independent flux
    *args: list of parameters.
    **kwargs: dic of parameters.
    Returns
    -------
    f: RSS + Penalty score (When out side of the lower and upper boundaries)


    See Also
    --------
    fit_r_mdv_scipy

    """
    Rm_initial = kwargs['Rm_initial']
    stoichiometric_num = kwargs['stoichiometric_num']
    reaction_num = kwargs['reaction_num']
    reac_met_num = kwargs['reaction_num']
    matrixinv = kwargs['matrixinv']
    experiments = kwargs['experiments']
    mdv_exp = numpy.array(kwargs['mdv_exp'])
    mdv_use = kwargs['mdv_use']
    covinv = kwargs['covinv']
    lb = kwargs['lb']
    ub = kwargs['ub']
    calmdv = kwargs['calmdv']
    diffmdv = kwargs['diffmdv']

    #calmdv = kwargs['calmdv']

    Rm = numpy.array(list(Rm_initial))
    Rm[stoichiometric_num: reaction_num] = list(x)
    tmp_r = numpy.dot(matrixinv, Rm)

    g = numpy.hstack((numpy.array(lb) - tmp_r, tmp_r - numpy.array(ub)))
    sum = 0.0
    for i in g:
        if i > 0:
            sum = sum + i * 100000
            #print(i)
    fail = 0
    #Determination of MDV
    mdv_original = list(tmp_r)
    for experiment in sorted(experiments.keys()):
        target_emu_list = experiments[experiment]['target_emu_list']
        mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']

        #
        #
        #
        if experiments[experiment]['mode'] == "ST":
            mdv_original_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)
        elif experiments[experiment]['mode'] == "INST":
            y0temp = experiments[experiment]['y0']
            timepoints = experiments[experiment]['timepoint']
            mdv_original_temp, mdv_hash = diffmdv(list(tmp_r), [], timepoints, target_emu_list, mdv_carbon_sources, y0temp)
        mdv_original.extend(mdv_original_temp)

    mdv = numpy.array([y for x, y in enumerate(mdv_original) if mdv_use[x] != 0])
    res = mdv_exp - mdv
    f = numpy.dot(res, numpy.dot(covinv, res))
    if experiments[experiment]['mode'] == "INST":
        print(f)
    return f + sum




def search_ci_gn(configure, experiments, numbers, vectors, matrixinv, func, flux, r_depended):
    """
    This fucntion is now obsolute.
    Low level function for finding confidence interval by Gauss-Newton method descrinbed in Antoniewicz MR, Kelleher JK, Stephanopoulos G (2006)
    Determination of confidence intervals of metabolic fluxes estimated from stable isotope measurements. Metab Eng 8(4):324-337


    """

    if isinstance(func, dict):
        calmdv = func["calmdv"]
        diffmdv = func["diffmdv"]
    else:
        exec(func)
    def null(A, eps=1e-4):
        """
        http://mail.scipy.org/pipermail/scipy-user/2005-June/004650.html
        """
        u, s, vh = numpy.linalg.svd(A)
        n = A.shape[1]   # the number of columns of A
        if len(s)<n:
            expanded_s = numpy.zeros(n, dtype = s.dtype)
            expanded_s[:len(s)] = s
            s = expanded_s
        null_mask = (s <= eps)
        null_space = numpy.compress(null_mask, vh, axis=0)
        return numpy.transpose(null_space)

    #Set max number of iteration ###############################################
    if 'iteration_max' in configure:
        iteration_max = configure['iteration_max']
    else:
        iteration_max = 1000

    #
    if 'step_for_cisearch' in configure:
        rss_increase_expected_initial = configure['step_for_cisearch']
    else:
        rss_increase_expected_initial = 0.001
    #
    # スパイクを検出する閾値
    #
    gap_thres_initial = rss_increase_expected_initial * 10.0

    #Set callbacklevel
    if configure.has_key('callbacklevel'):
        callbacklevel = configure['callbacklevel']
    else:
        callbacklevel = 0
    #
    # フラックスの下限値
    #
    thres_lower = -1.0e-8
    #
    #||delta_u||の変動上限
    #
    thres_deltau = 1.0e7
    #
    # epsの変化量
    #
    eps_increase = 3.0
    eps_decrease = 0.95
    #
    # h の初期値
    #
    h_initial = 1.0e-5

    #
    # epsの初期値
    #
    eps_initial = 1.0e-10
    delta_jacobi_initial = 0.0005
    #
    #
    # number of independent flux################################################
    independent_number = numbers['independent_number']
    stoichiometric_num = numbers['independent_start']
    reaction_num = numbers['independent_end']
    total_number = numbers['total_number']
    searching = numbers['searching']
    RSS_threshold = numbers['rss_limit']
    print("RSS_threshold", RSS_threshold)

    # zero independent flux
    Rm_ind = [flux[group][id]["value"] for (group, id) in vectors['independent_flux']]
    #boundaries
    lb = list(vectors["lb"])
    ub = list(vectors["ub"])
    independent_lb = list(vectors["independent_lb"])
    independent_ub = list(vectors["independent_ub"])
    non_fixed_reactions = list(vectors["non_fixed_reactions"])
    #
    #experiment の数だけ並べたMDVのベクトルを作成
    #
    mdv_exp_original = list(vectors["value"])
    mdv_std_original = list(vectors["stdev"])
    mdv_use = list(vectors["use"])
    for experiment in sorted(experiments.keys()):
        mdv_exp_original.extend(experiments[experiment]['mdv_exp_original'])
        mdv_std_original.extend(experiments[experiment]['mdv_std_original'])
        mdv_use.extend(experiments[experiment]['mdv_use'])
    mdv_exp = numpy.array([y for x, y in enumerate(mdv_exp_original) if mdv_use[x] != 0])
    spectrum_std = numpy.array([y for x, y in enumerate(mdv_std_original) if mdv_use[x] != 0])
    #
    # Covariance matrix
    #
    covinv = numpy.zeros((len(spectrum_std),len(spectrum_std)))
    for i, std in enumerate(spectrum_std):
        covinv[i,i] = 1.0/(std**2)
    #
    # RSS of best fit result
    #
    Rm_bestfit = list(vectors["Rm_initial"])
    Rm_bestfit[stoichiometric_num:reaction_num] = list(Rm_ind)
    tmp_r_bestfit = numpy.dot(matrixinv, Rm_bestfit)
    mdv_bestfit = list(tmp_r_bestfit)
    for experiment in sorted(experiments.keys()):
        target_emu_list = experiments[experiment]['target_emu_list']
        mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
        if experiments[experiment]['mode'] == "ST":
            mdv_temp, mdv_hash = calmdv(list(tmp_r_bestfit), target_emu_list, mdv_carbon_sources)
        elif experiments[experiment]['mode'] == "INST":
            timepoints = experiments[experiment]['timepoint']
            mdv_temp, mdv_hash = diffmdv(list(tmp_r_bestfit), [], timepoints, target_emu_list, mdv_carbon_sources)
        mdv_bestfit.extend(mdv_temp)

    mdv_bestfit = numpy.array([y for x, y in enumerate(mdv_bestfit) if mdv_use[x] != 0])
    res_bestfit = mdv_exp - mdv_bestfit
    RSSbestfit = numpy.dot(res_bestfit, numpy.dot(covinv, res_bestfit))
    #
    # Check best fit fluxes
    #
    if callbacklevel >= 3:
        print('RSS of best fit is', RSSbestfit, '. Start searching confidence interval.')
    if len([x for x in non_fixed_reactions if tmp_r_bestfit[x] < thres_lower]) > 0:
        if callbacklevel >= 3:
            print('Big problem!!!! Original flux has negative values')
            print([x for x in non_fixed_reactions if tmp_r_bestfit[x] < thres_lower] )
    else:
        if callbacklevel >= 3:
            print('Original flux has no negative value.')
    #
    # List for storing results
    #
    rss_all = []
    flux_all = []
    raw_flux_all = []
    #
    # Set Ni
    #
    Ni = numpy.zeros(independent_number)
    Ni = matrixinv[searching,stoichiometric_num:reaction_num]

    for direction in ['reverse', 'forward']:
        #
        # 5回チャレンジ
        # 毎回切り替わるのは
        eps_initial_temp = eps_initial * 1.0
        rss_increase_expected = rss_increase_expected_initial * 1.0

        for trial in range(5):
            #
            # Initializing parameters
            #
            h = h_initial
            Endpoint = {'text':'Ended at intial point', 'value' :9}
            gap_thres = gap_thres_initial
            end_counter = 0
            spike_counter = 0
            eps_thres = eps_initial_temp * 1.0
            #
            # Record initial (best fit) point
            #
            rss_one_direction = [RSSbestfit]
            flux_one_direction = [tmp_r_bestfit[searching]]
            move_record = []

            Rm_temp = list(Rm_bestfit)
            #
            # List for recording active constraints
            #
            active_constraint_array = [0] * len(tmp_r_bestfit)
            #
            #探索開始
            #
            for iteration_number in range(iteration_max):
                #
                # 1. simulare mesurements x for current fluxes
                #
                tmp_r = numpy.dot(matrixinv, Rm_temp)
                mdv_templete = list(tmp_r)
                for experiment in sorted(experiments.keys()):
                    target_emu_list = experiments[experiment]['target_emu_list']
                    mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
                    if experiments[experiment]['mode'] == "ST":
                        mdv_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)
                    elif experiments[experiment]['mode'] == "INST":
                        timepoints = experiments[experiment]['timepoint']
                        mdv_temp, mdv_hash = diffmdv(list(tmp_r), [], timepoints, target_emu_list, mdv_carbon_sources)
                    mdv_templete.extend(mdv_temp)
                mdv_templete = numpy.array([y for z,y in enumerate(mdv_templete) if mdv_use[z] != 0])
                #
                # 2. Calculate the variance-weighted sum of squared residuals Phi(u) or RSScurrent.
                #
                res = (mdv_templete - mdv_exp)
                RSScurrent = numpy.dot(numpy.array(res), numpy.dot(covinv, numpy.array(res)))
                #
                # Set temporary step h
                #
                h_temp = h * 1.0
                delta_jacobi = delta_jacobi_initial * 1.0
                #
                # Start finding feasible direction or delta_u
                #
                for modification_of_h in range(10):
                    #
                    # 2. Calculate the sensitivity matrix dx/du
                    #
                    dxdu = numpy.zeros((len(mdv_templete), independent_number))
                    #
                    # Set delta_jacobi.
                    # Avoid too large and too small

                    if h_temp > delta_jacobi_initial:
                        delta_jacobi = delta_jacobi_initial *1.0
                    elif h_temp < delta_jacobi_initial * 0.01:
                        delta_jacobi = delta_jacobi_initial * 0.01
                    else:
                        delta_jacobi = h_temp
                    #if direction == 'reverse':
                    #delta_jacobi = delta_jacobi * -1.0
                    #
                    # Calculation of dx/du
                    #
                    for i in range(independent_number):
                        #
                        # dx/dui
                        #
                        Rm_inter = numpy.copy(Rm_temp)
                        Rm_inter[stoichiometric_num + i] = Rm_temp[stoichiometric_num + i] + delta_jacobi
                        tmp_r_inter = numpy.dot(matrixinv, Rm_inter)
                        mdv_inter = list(tmp_r_inter)
                        for experiment in sorted(experiments.keys()):
                            target_emu_list = experiments[experiment]['target_emu_list']
                            mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
                            if experiments[experiment]['mode'] == "ST":
                                mdv_temp, mdv_hash = calmdv(list(tmp_r_inter), target_emu_list, mdv_carbon_sources)
                            elif experiments[experiment]['mode'] == "INST":
                                timepoints = experiments[experiment]['timepoint']
                                mdv_temp, mdv_hash = diffmdv(list(tmp_r_inter), [], timepoints, target_emu_list, mdv_carbon_sources)
                            mdv_inter.extend(mdv_temp)
                        mdv_inter = numpy.array([mdv_inter[z] for z in range(len(mdv_inter)) if mdv_use[z] != 0])

                        dxdu[:,i] = ((mdv_inter - mdv_templete)/delta_jacobi)
                    #
                    #
                    temp_r_history = []
                    #
                    # Set initial value tmp_r2
                    #
                    tmp_r2 = [-1.0] * len(tmp_r)
                    #
                    # Modification of dxdu until all flux > 0
                    #
                    # Initialize delta_u
                    #
                    delta_u = numpy.zeros(independent_number)
                    for modification_of_dxdu in range(100):
                        active_constraint = [x for x, n in enumerate(active_constraint_array) if n > 0]
                        #
                        #  Avoid extreme eps_thres level
                        #
                        if eps_thres > 0.001:
                            eps_thres = 0.001 * numpy.random.random()
                        if eps_thres < 1e-10:
                            eps_thres = 1e-10
                        #
                        if h_temp < 1e-6:
                            h_temp = h * 1.0
                        #
                        # 4. Calculate the variance-weighted sum of squared residuals Phi(u).
                        # 5. Calculate matrix A and vector b using Eq. (24).
                        #
                        mdv_gap = mdv_templete - mdv_exp
                        J = numpy.dot(numpy.transpose(dxdu),numpy.dot(covinv,mdv_gap))
                        H = numpy.dot(numpy.transpose(dxdu),numpy.dot(covinv,dxdu))
                        b = numpy.zeros(independent_number+1+len(active_constraint))
                        b[0:independent_number] = J * -2.0
                        if direction == 'forward':
                            b[independent_number] = h_temp
                        elif direction == 'reverse':
                            b[independent_number] = h_temp * -1.0
                        #A
                        #
                        #
                        #
                        try:
                            V2 = null(H, eps = eps_thres)
                        except:
                            #
                            # When H has bad condition
                            #
                            eps_thres = eps_thres * eps_increase
                            if (callbacklevel >= 3):
                                print('H has bad condition for V2. Current eps is ', eps_thres, 'at iteration(s)', modification_of_h, modification_of_dxdu)
                            continue

                        Htemp = H + eps_thres * numpy.dot(V2, numpy.transpose(V2))
                        #
                        #
                        A = numpy.zeros((independent_number + 1 + len(active_constraint), independent_number + 1 + len(active_constraint)))
                        A[0:independent_number,independent_number] = list(Ni)
                        A[independent_number, 0:independent_number] = list(Ni)
                        A[0:independent_number,0:independent_number] = Htemp * 2.0
                        #
                        # Set active constraints
                        #
                        for i, reaction in enumerate(active_constraint):
                            Niac = matrixinv[reaction,stoichiometric_num:reaction_num]
                            A[0:independent_number,independent_number + i + 1] = list(Niac)
                            A[independent_number + i + 1, 0:independent_number] = list(Niac)
                            b[independent_number + i + 1] = (tmp_r[reaction] - 1.0e-8) * -1.0
                        #
                        # 6. Calculate Dd using Eq. (26). Additional inequality constraints may be included here.
                        #
                        try:
                            #solveはなんかおかしい
                            delta_d_solve = numpy.linalg.solve(A,b)
                            #Ainv = numpy.linalg.inv(A)
                            #delta_d_inv = numpy.dot(numpy.linalg.inv(A), b)
                        except:
                            #
                            # When H has bad condition
                            #
                            eps_thres = eps_thres * eps_increase
                            if (callbacklevel >= 3):
                                print('H has bad condition. Current eps is ', eps_thres, 'at iteration(s)', modification_of_h, modification_of_dxdu)
                            continue
                        #
                        # Set delta_u
                        #
                        delta_u_solve = delta_d_solve[0:independent_number]
                        #delta_u_inv = delta_d_inv[0:independent_number]
                        #if (callbacklevel > 3):
                        #    print "{0:5d} {1:5d} {2:5d} eps:{3:10.9f} u:{4:10.3f} u_inv:{5:10.3f} diff:{6:10.3f} cond:{7:10.3f}".format(trial, modification_of_h, modification_of_dxdu, eps_thres, sum(delta_u_solve* delta_u_solve),sum(delta_u_inv * delta_u_inv), sum((delta_u_solve - delta_u_inv)*(delta_u_solve - delta_u_inv)), numpy.linalg.cond(A))
                        #
                        # chech delta_u
                        # Discard when length of delta_u is too large
                        #if sum(delta_u_solve* delta_u_solve) > sum(delta_u_inv * delta_u_inv):
                        #    delta_u = numpy.copy(delta_u_inv)
                        #else:
                        #    delta_u = numpy.copy(delta_u_solve)
                        delta_u = numpy.copy(delta_u_solve)
                        if sum(delta_u * delta_u) > thres_deltau:
                            eps_thres = eps_thres * eps_increase
                            if (callbacklevel >= 3):
                                print("Too large delta u:", sum(delta_u * delta_u),'Eps is increased to', eps_thres)
                            continue
                        #
                        # Calculation of tmp_r2
                        #
                        Rm2 = list(Rm_temp)
                        Rm2[stoichiometric_num:reaction_num] = Rm_temp[stoichiometric_num:reaction_num] + delta_u[:]
                        tmp_r2 = numpy.dot(matrixinv, Rm2)
                        #
                        # Check movement
                        #
                        diff = tmp_r2[searching]-tmp_r[searching]
                        #
                        # If movement of target flux is unequal step size h
                        #
                        if  not abs(h_temp * 0.8) < abs(diff) < abs(h_temp * 1.2):
                            if (callbacklevel >= 3):
                                print("Odd diff another diff is tested :", diff, eps_thres, h_temp)
                            continue
                        #
                        # Check direction of near zero fluxes.
                        #
                        near_zero_reactions = set([x for x in non_fixed_reactions if tmp_r[x] < 0.01])
                        if len(near_zero_reactions) > 0:
                            less_than_zero_reactions = set([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower])
                            active_constraint_temp = near_zero_reactions.intersection(less_than_zero_reactions)
                            if len(active_constraint_temp) > 0:
                                for reaction_number in active_constraint_temp:
                                    #
                                    # Ignore if reaction of interest is searching reactions
                                    #
                                    #if reaction_number == reaction_number_for_searching_forward:
                                    #    continue
                                    #if reaction_number == reaction_number_for_searching_reverse:
                                    #    continue
                                    #
                                    # Add to active constraints
                                    #
                                    active_constraint_array[reaction_number] = 20
                                    #active_constraint_array[reaction_number] = active_constraint_array[reaction_number] + 20
                                if (callbacklevel >= 3):
                                    print('Set near zero flux to active constraint(s)', active_constraint_temp, 'at iteration', modification_of_dxdu)
                                    print('Before move', [tmp_r[x] for x in active_constraint_temp])
                                    print('After move', [tmp_r2[x] for x in active_constraint_temp])
                                    print('Counter', [active_constraint_array[x] for x in active_constraint_temp])
                                continue
                        #
                        # 進んだ先を記録する。
                        #
                        flag = True
                        #
                        # Modify dxdu to avoid negative flux values
                        #
                        if len([1 for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) > 0:
                            flag = False
                            f = [x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]

                            for i in f:
                                tmp_r_d = r_depended[i]
                                for line in tmp_r_d:
                                    #
                                    # reduce dx/du related to reactions with negative values
                                    #
                                    dxdu[:,line] = dxdu[:,line]/2.0
                                    if (callbacklevel >= 3):
                                            print('dudx modified', i, f, modification_of_dxdu)
                            #
                            # If independent fluxes are a pair of reversible reactions
                            #
                            if (len(f) == 1) and (len(tmp_r_d)==2) and (len(temp_r_history) > 0):
                                if temp_r_history[f[0]] >= tmp_r2[f[0]]:#前回と差がないとき
                                    flag = True

                        #過去を記録する
                        temp_r_history = list(tmp_r2)
                        #
                        if flag == True:
                            break
                    ######################################
                    # End of Modification of dxdu:
                    ######################################
                    #
                    # If there are fluxes with negative values
                    #
                    if len([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) > 0:
                        if callbacklevel >= 3:
                            print("Enter super mode!!!")
                            print(trial, "th tiral")
                            print('modification_of_h', modification_of_h)
                            print('modification_of_dxdu', modification_of_dxdu)
                            print('eps_thres',eps_thres)
                            print('h_temp',h_temp)
                            print('sum_deltau', sum(delta_u * delta_u))
                            print('diff',abs(diff))
                            print('a',[x for x in non_fixed_reactions if tmp_r2[x] < thres_lower])
                            print('Before move', [tmp_r[x] for x in active_constraint_temp])
                            print('After move', [tmp_r2[x] for x in active_constraint_temp])
                            print('Counter', [active_constraint_array[x] for x in active_constraint_temp])
                        #
                        # chech delta_u
                        # Discard when length of delta_u is too large
                        #
                        if sum(delta_u * delta_u) > 0.0: #そもそもdelta_d が計算されない場合要
                            if sum(delta_u * delta_u) > thres_deltau:
                                eps_thres = eps_thres * eps_increase
                                delta_jacobi =  delta_jacobi * 1.2
                                if (callbacklevel >= 3):
                                    print("Too large delta u:", sum(delta_u * delta_u),'Delta Jacobi is increased to', delta_jacobi)
                        #
                        # Reduce step size
                        #　最初4回はhを減らす。
                        if  modification_of_h < 4:
                            h_temp = h / (2.0 * (modification_of_h + 1.0))
                            if (callbacklevel >= 2):
                                print("Cannot find forward direction, h is reduced to:", h_temp)
                        #

                        if 2 < modification_of_h < 4:
                            #
                            # Add fluxes in problem to active constraints
                            # Check whether tmp_r2 == [-1]*
                            #
                            if len([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) < len(non_fixed_reactions):
                                #
                                # Set active constrains
                                #
                                active_constraint_temp = [x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]

                                for reaction_number in active_constraint_temp:
                                    active_constraint_array[reaction_number] = active_constraint_array[reaction_number] + 10

                                    if (callbacklevel >= 3):
                                        print('Set NON zero flux to active constraint(s)', active_constraint_temp)
                                        print('Before move', [tmp_r[x] for x in active_constraint_temp])
                                        print('After move', [tmp_r2[x] for x in active_constraint_temp])
                                        print('Counter', [active_constraint_array[x] for x in active_constraint_temp])
                        continue
                    break
                ######################################
                # End of Modification of h:
                ######################################
                #
                # Finish iteration when forward direction can not be found,
                #
                if len([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) > 0:
                    if callbacklevel >= 1:
                        print(trial, 'th trial was failed at', iteration_number)
                    break

                #
                # Calculation of RSS after movement
                #
                #print "a"
                mdv2 = list(tmp_r2)
                for experiment in sorted(experiments.keys()):
                    target_emu_list = experiments[experiment]['target_emu_list']
                    mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
                    if experiments[experiment]['mode'] == "ST":
                        mdv_temp, mdv_hash = calmdv(list(tmp_r2), target_emu_list, mdv_carbon_sources)
                    elif experiments[experiment]['mode'] == "INST":
                        timepoints = experiments[experiment]['timepoint']
                        mdv_temp, mdv_hash = diffmdv(list(tmp_r2), [], timepoints, target_emu_list, mdv_carbon_sources)
                    mdv2.extend(mdv_temp)
                mdv2 = numpy.array([mdv2[z] for z in range(len(mdv2)) if mdv_use[z] != 0])

                res2 = (mdv2 - mdv_exp)
                rss2 = numpy.dot(numpy.array(res2), numpy.dot(covinv, numpy.array(res2)))
                #
                # Claculation of next step size
                #
                a1 = numpy.dot(numpy.transpose(delta_u),J) * 2.0
                a2 = numpy.dot(numpy.transpose(delta_u),numpy.dot(H, delta_u))
                deltaRSS_approx = a1 * h_temp + a2 * h_temp * h_temp
                a3 = abs((rss2-RSScurrent) - deltaRSS_approx)/(h_temp*h_temp*h_temp)
                h2 = ((rss_increase_expected/ a3)**(1.0/3.0))
                #
                # Dsiplay iterations
                #
                if (callbacklevel >= 2):
                    diff = tmp_r2[searching]-tmp_r_bestfit[searching]
                    print("{0:5d} h:{1:10.9f} h':{2:10.9f} h2:{3:10.9f} rss2:{4:10.9f} diff:{5:10.9f} esp:{6:12.11f} gap:{7:10.1f} move:{8:+8.7f} loop:{9:3d} end counter:{10:3d}".format(iteration_number, h, h_temp, h2, rss2, diff, eps_thres, (rss2-RSScurrent)/deltaRSS_approx, rss2-RSScurrent, modification_of_dxdu, end_counter),)
                    print("active constraints:", [i for i, x in enumerate(active_constraint_array) if x > 0])
                #
                # Check spikes
                # Check coutinous spikes
                #
                if spike_counter > 10:
                    #
                    # Reset counter
                    #
                    spike_counter = 0
                    #
                    # Relax threshold for spike detection
                    #
                    gap_thres = gap_thres * 2
                    if (callbacklevel >= 3):
                        print('Threshold for spike detection is increased to:',gap_thres)
                #
                # Spike detection. If
                # RSS is suddenly increased : rss2-RSScurrent > gap_thres
                # and a lager gap between predicted and acutual movement : abs((rss2-RSScurrent)/deltaRSS_approx) > 1000
                # This threshold is very heulistic
                if (rss2-RSScurrent) > gap_thres:
                    #
                    # Stabilize H by increase of eps_thres
                    #
                    eps_thres = eps_thres * eps_increase
                    if (callbacklevel >= 3):
                        print("May be spikes. This iteration is rejected. eps is changed to ", eps_thres, "Spike counter is ",spike_counter, gap_thres )
                    #
                    #　Add counter
                    #
                    spike_counter = spike_counter + 1
                    continue
                #
                # When movement without spike, eps_thres are decreaed
                #
                #eps_thres = eps_thres * eps_decrease
                #
                #
                #
                eps_thres = eps_initial_temp * 1.0

                if iteration_number > 11:
                    average = sum(move_record[-10:])/10.0
                    stdev = numpy.std(move_record[-10:])
                    if 0.2 < abs(stdev/average) < 0.5:
                        eps_thres = eps_thres * 5
                    elif 0.5 <= abs(stdev/average) < 1.0:
                        eps_thres = eps_thres * 25
                    elif 1.0 < abs(stdev/average):
                        eps_thres = eps_thres * 100
                #
                # Reset counter
                #
                spike_counter = 0
                #
                # Control step size at first iteration
                #
                if iteration_number == 0:
                    if abs(h2) > 0.0001:
                        h = 0.0001
                    else:
                        h = h2 * 1.0
                    continue
                #
                # Avoid too small step size
                #
                if h2 < 1e-6:
                    h2 = 1e-3
                #
                # Record rss and movement
                #
                move_record.append(rss2-RSScurrent)
                #
                rss_one_direction.append(rss2)
                flux_one_direction.append(tmp_r2[searching])
                #
                # Check end counter
                #
                if rss2 > RSS_threshold:
                    end_counter = end_counter + 1
                else:
                    end_counter = 0
                if end_counter > 100:
                    Endpoint = {'text':'Over threshold at '+str(iteration_number)+" th iteration", 'value' :0}
                    break
                #
                # Go to nest step
                #
                Rm_temp = numpy.copy(Rm2)
                h = h2 * 1.0
                #
                # Reduce active constraint counter
                #
                for reaction_number in [i for i, x in enumerate(active_constraint_array) if x > 0]:
                    active_constraint_array[reaction_number] = active_constraint_array[reaction_number]-1
                #
                # When last iteraction
                #
                if iteration_number == iteration_max - 1:
                    Endpoint = {'text':'Iteration max', 'value' :9}
                    break
                #
                # 見ている反応がゼロよりも小さくなったら終了
                #
                if tmp_r2[searching] < thres_lower:
                    Endpoint = {'text':'Target flux is less than zero', 'value' :8}
                    break
                #
                # 見ている反応がゼロよりも小さくなったら終了
                #
                if tmp_r2[searching] > 1000:
                    Endpoint = {'text':'Target flux is over 1000', 'value' :8}
                    break
                elif tmp_r2[searching] < 0.001:
                    Endpoint = {'text':'Target flux is less than 0.001', 'value' :8}
                    break
            ######################################
            # End of iteration_number:
            ######################################
            #
            # Finish iteration when forward direction can not be found,
            #
            if len([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) > 0:

                eps_initial_temp = eps_initial * (10.0 ** (trial + 1 ))
                rss_increase_expected = rss_increase_expected * 0.5
                if (callbacklevel >= 2):
                    print('Start again*********************************************************************')
                    print('eps_initial_temp and rss_increase_expected was set to', eps_initial_temp, rss_increase_expected)
                continue
            if (callbacklevel >= 2):
                print("finished at ", trial, "th trial", Endpoint)
            break
        ######################################
        # End of trial:
        ######################################
        if len([x for x in non_fixed_reactions if tmp_r2[x] < thres_lower]) > 0:
            if (callbacklevel >= 3):
                print('Cannot find forward direction')
            Endpoint = {'text':'Cannot find forward direction', 'value' :1}
        if direction == 'reverse':
            flux_reverse = [flux_one_direction[x] for x in range(len(flux_one_direction)) if rss_one_direction[x] < RSS_threshold][-1]
            rss_one_direction.reverse()
            flux_one_direction.reverse()
            rss_all.extend(rss_one_direction)
            flux_all.extend(flux_one_direction)
            state_reverse = Endpoint

            #flux_reverse = list(tmp_r)
        elif direction == 'forward':
            flux_forward = [flux_one_direction[x] for x in range(len(flux_one_direction)) if rss_one_direction[x] < RSS_threshold][-1]
            rss_all.extend(rss_one_direction[1:])
            flux_all.extend(flux_one_direction[1:])
            state_foward = Endpoint
    #####################################3
    # End of ['reverse', 'forward']:
    ######################################
    return(flux_all, rss_all, [], flux_reverse, state_reverse, flux_forward, state_foward)

