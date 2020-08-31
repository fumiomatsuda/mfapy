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
    Funcition to generate randomized initial flux dixtribution using scipy.optimize.minimize SLSQP

    Parameters
    ----------
    numbers: "model.numbers" including various number related data of the model.
    vectors: "model.vector" including various vector related data of the model.
    matrixinv: "model.matrixinv" is a inversed matrix for flux calculation.
    template: Dictionary of metabolic state. When template is available, metabolic state most similar to the template is generated. The function is used in the grid search.
    initial_search_iteration_max: "configure["initial_search_iteration_max"]". Maximal number of interations (steps) allowed in each task to find feasible initial metabolic flux distribution.
    method: "fitting" is only available.

    Returns
    -------
    tmp_r: (list) metabolic state data (tmp_r = numpy.dot(matrixinv, Rm_temp)
    Rm_temp: (list) metabolic state vector
    Rm_ind: (list) independent flux vector
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
        if callbacklevel > 1:
            print("mkl-service is not installed this python!")
    # zero independent flux
    Rm_ind = list(numpy.zeros(independent_number))
    #boundaries
    lb = list(vectors["lb"])
    ub = list(vectors["ub"])
    independent_lb = list(vectors["independent_lb"])
    independent_ub = list(vectors["independent_ub"])

    tmp_r = []
    result_Rm = []
    result_ind = []
    message = "Initial state"
    try:
        for j in range(3):
            message = "Initial state"

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
                message = "Determined"
                break
            else:
                message = "Failed"


    except Exception as e:
        message = e
    else:
        pass
    finally:
        return(tmp_r, result_Rm, result_ind, message)




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
    #
    # Mas=number of MKL thread control
    #
    try:
        import mkl
        mkl.set_num_threads(1)
    except:
        if callbacklevel > 1:
            print("mkl-service is not installed this python!")
    #
    # Initial state
    #
    state = "Initial state"
    kai = -1.0
    opt_flux = []
    result_ind = []

    try:
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
        #try:
        ##################################################################
        if (callbacklevel >= 4):
            print("Fitting Start in fit_r_mdv_scipy using ", method)
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
        parameters['callbacklevel'] = callbacklevel

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

        if (callbacklevel >= 4):
            print("Fitting was successfully finished. RSS = ", kai)

        #Optimized flux distribution
        Rm_opt = numpy.array(list(vectors["Rm_initial"]))
        result_Rm = numpy.array(list(vectors["Rm_initial"]))
        result_Rm[numbers['independent_start']: numbers['independent_end']] = result_ind[:]

        opt_flux = numpy.dot(matrixinv, numpy.array(result_Rm))
    except Exception as e:
        state = e
    else:
        pass
    finally:
        return(state, kai, opt_flux, result_ind)


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
    #
    # Mas=number of MKL thread control
    #
    try:
        import mkl
        mkl.set_num_threads(1)
    except:
        if callbacklevel > 1:
            print("mkl-service is not installed this python!")
    #
    # Initial state
    #
    state = "Initial state"
    kai = -1.0
    opt_flux = []
    result_ind = []

    try:
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
        #try:
        ##################################################################
        if (callbacklevel >= 4):
            print("Fitting Start infit_r_mdv_nlopt using ", method)
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
        parameters['callbacklevel'] = callbacklevel
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


        if (callbacklevel >= 4):
            print("Fitting was successfully finished. RSS = ", kai)

        #Optimized flux distribution
        Rm_opt = numpy.array(list(vectors["Rm_initial"]))
        result_Rm = numpy.array(list(vectors["Rm_initial"]))
        result_Rm[numbers['independent_start']: numbers['independent_end']] = result_ind[:]

        opt_flux = numpy.dot(matrixinv, numpy.array(result_Rm))

        #return(state, kai, opt_flux, result_ind)


        #return(state,-1,[],[])
    except Exception as e:
        state = e
    else:
        pass
    finally:
        return(state, kai, opt_flux, result_ind)


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
    if (callbacklevel >= 4):
        print("##Start GN_CRS2_LM method######################################################################")
    state, kai, flux, Rm_ind_sol = fit_r_mdv_nlopt(configure, experiments, numbers, vectors, matrixinv, func, flux, method = "GN_CRS2_LM")

    for k in range (number_of_repeat):
        if (callbacklevel >= 4):
            print("Deep",k,"Start SLSQP method######################################################################")
        state, kai, flux, Rm_ind_sol = fit_r_mdv_scipy(configure, experiments, numbers, vectors, matrixinv, func, flux, method = "SLSQP")
        if (callbacklevel >= 4):
            print("Deep",k,"Start LN_PRAXIS method##################################################################")
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
    callbacklevel = kwargs['callbacklevel']

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
        if callbacklevel >= 2:
            print("RSS:", f)
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
    callbacklevel = kwargs['callbacklevel']

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
        if callbacklevel >= 4:
            print("RSS:", f)

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
    callbacklevel = kwargs['callbacklevel']

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
        if callbacklevel >= 4:
            print("RSS:", f)
    return f + sum



