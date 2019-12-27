#-------------------------------------------------------------------------------
# Name:        metabolicmodel.py
# Purpose:     MetabolicModel class in mfapy
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------

import numpy as numpy
import scipy as scipy
import itertools, re, math
from . import mdv
from . import optimize
from . import carbonsource
import copy,time
#from numba import jit

#from assimulo.problem import Explicit_Problem
#from assimulo.solvers import CVode

class MetabolicModel:
    """
    Class of Metabolic model

    Instances of this class have fuctions to analyze mass spec. data derived from 13C metabolic flux analysis.
    Other classes including MDVData, CarbonSouece are constructed from a instances of this class


    Parameters
    ----------
    reactions : Dictionary describing metabolite reactions
    reversible : Dictionary for defining reversible reactions
    metabolites : Dictionary including metabolite information
    target_fragments : Dictionary of target_fragments

    Returns
    -------
    model : instance of metabolic model


    See Also
    --------


    Examples
    --------
    >>> reactions, reversible, metabolites, target_fragments = LoadMetabolicModel("filename.txt')
    >>> model = MetabolicModel(reactions, reversible, metabolites, target_fragments)

    Attributes
    --------
        self.reactions = reactions
        self.metabolites = metabolites
        self.reversible = reversible
        self.target_fragments = target_fragments
        self.symmetry = {}
        self.carbon_source = {}
        self.experiments = {}

    """

    def __init__(self, reactions, reversible, metabolites, target_fragments, mode = "normal"):
        '''
        Generator of new instances
        '''
        #
        # Keep original information
        #
        self.reactions = copy.deepcopy(reactions)
        self.metabolites = copy.deepcopy(metabolites)
        self.reversible = copy.deepcopy(reversible)
        self.target_fragments = copy.deepcopy(target_fragments)
        self.symmetry = {}
        self.carbon_source = {}
        #
        # MFA experimental sets
        #
        self.experiments = {}
        #
        # Class configuration data
        #
        self.configuration={
        'callbacklevel': 0,
        'default_reaction_lb': 0.001,
        'default_reaction_ub': 5000.0,
        'default_metabolite_lb': 0.001,
        'default_metabolite_ub': 500.0,
        'default_reversible_lb': -300,
        'default_reversible_ub': 300.0,
        'iteration_max': 1000,
        'initial_search_repeats_in_grid_search': 50,
        'initial_search_iteration_max': 1000,
        'grid_search_iterations':1, #170623 for Grid Maeda
        'add_naturalisotope_in_calmdv':"no", #180616
        'number_of_repeat':3,
        'ncpus':3,
        'ppservers': ("",),
        'odesolver': "scipy" # or "sundials"
        }
        #
        # Class configuration data
        # 0: No call back
        # 1: Simple report for grid search
        # 2: Detailed report for grid search
        # 3: Simple report for fitting_flux using parallel python
        # 4: Detailed report for fitting_flux using parallel python
        # 5: Simple report report for single analysis
        # 6: Detailed report for single analysis
        # 7: Detailed report for model constrution
        # 8: Debug
        # Symmetry metabolites
        #
        for id in self.metabolites:
            if self.metabolites[id]['symmetry'] == 'symmetry':
                self.symmetry[id] = 1
        #
        # Generate templete for carbon_sources dictionary
        #
        for id in self.metabolites:
            if self.metabolites[id]['carbonsource'] == 'carbonsource':
                self.carbon_source[id] = {
                    'size':self.metabolites[id]['C_number'],
                    'IDV':list(numpy.zeros(2**self.metabolites[id]['C_number']))
                }
        #
        # Initilize flux values
        #
        for id in self.reactions:
            self.reactions[id]['value'] = 0.0
            self.reactions[id]['stdev'] = 1.0
            self.reactions[id]['type'] = 'free'
            self.reactions[id]['reversible'] = 'no'
            if 'lb' not in self.reactions[id]:
                self.metabolites[id]['lb'] = self.configuration["default_reaction_lb"]
            if 'ub' not in self.reactions[id]:
                self.metabolites[id]['ub'] = self.configuration["default_reaction_ub"]
        #
        # Initilize metabolite conc values
        # Impremented in ver 1.03
        #
        for id in self.metabolites:
            self.metabolites[id]['value'] = 1.0
            self.metabolites[id]['stdev'] = 1.0
            self.metabolites[id]['type'] = 'fixed'
            if 'lb' not in self.metabolites[id]:
                self.metabolites[id]['lb'] = self.configuration["default_metabolite_lb"]
            if 'ub' not in self.metabolites[id]:
                self.metabolites[id]['ub'] = self.configuration["default_metabolite_ub"]
        #
        # Reversible reactions
        #
        for id in self.reversible:
            self.reversible[id]['type'] = 'free'
            self.reversible[id]['value'] = 0.0
            self.reversible[id]['stdev'] = 1.0
            #forward = self.reversible[id]['forward']
            #reverse = self.reversible[id]['reverse']
            #self.reactions[forward]['reversible'] = id
            #self.reactions[reverse]['reversible'] = id
            if 'lb' not in self.reversible[id]:
                self.reversible[id]['lb'] = self.configuration["default_reversible_lb"]
            if 'ub' not in self.reversible[id]:
                self.reversible[id]['ub'] = self.configuration["default_reversible_ub"]
        #
        # Define vector lengthes of mass spectra of oberved fragments
        #
        for id in self.target_fragments:
            atommap = self.target_fragments[id]['atommap']

            #
            # In the case of direct analysis of intermeidates, vector length is carbon number + 1
            #

            if self.target_fragments[id]['type'] == 'intermediate':
                metabolite, positions = atommap.split("_")
                self.target_fragments[id]['number'] = len(positions.split(":")) + 1
            #
            # In the case of amino acid analysis using GC-MS, vector length is (sum of carbon numbers) + 1
            #
            elif self.target_fragments[id]['type'] == 'gcms':
                number = 0
                for fragment in atommap.split("+"):
                    metabolite, positions = fragment.split("_")
                    number = number + len(positions.split(":"))
                self.target_fragments[id]['number'] = number + 1
            #
            # In the case of intermediate analysis using MS/MS, vector length is ([nl] + 1) * ([prod]+1)
            #
            elif self.target_fragments[id]['type'] == 'msms':
                precursor, nl, product = atommap.split("+")
                nl_positions = nl.split("_")[1]
                product_positions = product.split("_")[1]
                self.target_fragments[id]['number'] = (len(positions.split(":"))+1) * (len(positions.split(":"))+1)
            #
            # Set matrixes on isotope effects
            #
            formula = self.target_fragments[id]['formula']
            if not formula == "":
                num = self.target_fragments[id]['number']
                self.target_fragments[id]['natural_isotope_cancellation'] = mdv.INV_correcting(formula)[0:num,0:num]
                self.target_fragments[id]['natural_isotope_addition'] = mdv.transition_matrix(formula)[0:num,0:num]
        if mode=="debug":
            self.configuration["callbacklevel"] = 7

        #
        # Update the metabolic model to generate stoichiometry matrix
        #
        self.reconstruct()
        if mode=="debug":
            self.configuration["callbacklevel"] = 1

    def update(self):
        """
        Update the metabolic model to generate stoichiometry matrix. The metabolic model
        have to be updated when reaction types (fixed, free, fitting, etc) in the metabolic model
        are modified.
        When total number of metabolic reactions were changed, please perform self.reconstruction()


        Parameters
        ----------
        nothing


        Examples
        --------
        >>> model.update()

        See Also
        --------
        reconstruct

        """
        def rref(M):
            """
            Find the matrix in reduced row echelon form
            """
            if not M: return
            lead = 0
            rowCount = len(M)
            columnCount = len(M[0])
            for r in range(rowCount):
                if lead >= columnCount:
                    return M
                i = r
                while M[i][lead] == 0:
                    i += 1
                    if i == rowCount:
                        i = r
                        lead += 1
                        if columnCount == lead:
                            return M
                M[i],M[r] = M[r],M[i]
                lv = M[r][lead]
                M[r] = [ mrx / lv for mrx in M[r]]
                for i in range(rowCount):
                    if i != r:
                        lv = M[i][lead]
                        M[i] = [ iv - lv*rv for rv,iv in zip(M[r],M[i])]
                        for j in range(columnCount):
                            if abs(M[i][j]) < 0.000001:
                                M[i][j] = 0;
                lead += 1
            for r in range(rowCount):
                for c in range(columnCount):
                    if abs(M[r][c]) < 0.000001:
                        M[r][c] = 0.0
            return M

        #
        # List of reaction id
        #
        self.reaction_ids = [];
        self.metabolite_ids = [];
        self.reversible_ids = [];
        #
        #list of Rm_initial
        #
        Rm_initial = []
        #
        # row_names
        #
        self.matrix_row_names = [];
        #
        # dictionary for storing information for stoichiometry matrix
        #
        stoichiometry_matirx_hash = {}
        #
        # Set vectors
        #
        self.vector ={"independent_lb":[], "independent_ub":[], "lb":[], "ub":[], "use":[], "value":[], "stdev":[], "ids":[]}
        self.state_ids = []
        counter = 0

        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            # ignore unused reactions
            #
            if self.reactions[id]['use'] != 'use':
                continue
            #
            # Order of the reaction in self.reaction_ids is stored
            #
            self.reactions[id]["position_in_tmp_r"] = counter
            #
            # Reaction id is appended to self.reaction_ids
            #
            self.reaction_ids.append(id)
            #self.state_ids.append(("reaction", id))
            self.vector["ids"].append(("reaction", id))
            counter += 1
            #
            # Check fixed reactions
            #
            self.vector["value"].append(float(self.reactions[id]['value']))
            self.vector["lb"].append(self.reactions[id]['lb'])
            self.vector["ub"].append(self.reactions[id]['ub'])
            self.vector["stdev"].append(self.reactions[id]['stdev'])

            if self.reactions[id]['type'] == "fitting":
                self.vector["use"].append(float(1.0))
            else:
                self.vector["use"].append(float(0.0))
        for id in sorted(self.metabolites.keys(), key=lambda x: self.metabolites[x]['order']):
            #
            # ignore excreted and carbonsource metabolites
            #
            if self.metabolites[id]["excreted"] == "excreted":
                continue
            if self.metabolites[id]["carbonsource"] == "carbonsource":
                continue
            #
            # Order of the reaction in self.reaction_ids is stored
            #
            self.metabolites[id]["position_in_tmp_r"] = counter
            #
            # Reaction id is appended to self.metabolite_ids
            #
            self.metabolite_ids.append(id)
            #self.state_ids.append(("metabolite", id))
            self.vector["ids"].append(("metabolite", id))
            counter += 1
            #
            # Check fixed reactions
            #
            self.vector["value"].append(float(self.metabolites[id]['value']))
            self.vector["lb"].append(self.metabolites[id]['lb'])
            self.vector["ub"].append(self.metabolites[id]['ub'])
            self.vector["stdev"].append(self.metabolites[id]['stdev'])
            if self.metabolites[id]['type'] == "fitting":
                self.vector["use"].append(float(1.0))
            else:
                self.vector["use"].append(float(0.0))

        for id in sorted(self.reversible.keys(), key=lambda x: self.reversible[x]['order']):
            #
            # Order of the reaction in self.reaction_ids is stored
            #
            self.reversible[id]["position_in_tmp_r"] = counter
            #
            # Reaction id is appended to self.metabolite_ids
            #
            self.reversible_ids.append(id)
            #self.state_ids.append(("reversible", id))
            self.vector["ids"].append(("reversible", id))
            counter += 1
            #
            # Check fixed reactions
            #
            self.vector["value"].append(float(self.reversible[id]['value']))
            self.vector["lb"].append(self.reversible[id]['lb'])
            self.vector["ub"].append(self.reversible[id]['ub'])
            self.vector["stdev"].append(self.reversible[id]['stdev'])
            if self.reversible[id]['type'] == "fitting":
                self.vector["use"].append(float(1.0))
            else:
                self.vector["use"].append(float(0.0))

        #Searching metabolites
        #Set for keeping visited metabolite names
        #
        intermediates = set()
        for id in self.reaction_ids:
            #
            # get reaction formula and substrate, product names
            #
            reaction = self.reactions[id]['stoichiometry']
            reaction.replace(" ","")
            substrate, product = reaction.split("-->")
            #
            # Substrate
            #
            for metabolite in substrate.split("+"):
                #
                # Check coefficient: {0.5}Pyr
                #
                if len(metabolite.split("}")) == 1:
                    # no coefficient
                    metabolite_name = metabolite
                elif len(metabolite.split("}")) == 2:
                    # coefficient data
                    stnum, metabolite_name = metabolite.split("}")
                #
                # Ignore if carbon source metabolites
                #
                if self.metabolites[metabolite_name]['carbonsource'] == 'carbonsource':
                    continue
                #
                # Ignore if excreted metabolites
                #
                if self.metabolites[metabolite_name]['excreted'] == 'excreted':
                    continue
                #
                # metabolite_name is visited
                #
                intermediates.update([metabolite_name])
            #
            # Product
            #
            for metabolite in product.split("+"):
                if len(metabolite.split("}")) == 1:
                    metabolite_name = metabolite
                elif len(metabolite.split("}")) == 2:
                    stnum, metabolite_name = metabolite.split("}")
                if self.metabolites[metabolite_name]['carbonsource'] == 'carbonsource':
                    continue
                if self.metabolites[metabolite_name]['excreted'] == 'excreted':
                    continue
                intermediates.update([metabolite_name])
        #
        # Initialize stoichiometry_matirx_hash for storing information for stoichiometry matrix
        #
        for metabolite in intermediates:
            stoichiometry_matirx_hash[metabolite] = {}
            for idtemp in self.reaction_ids:
                stoichiometry_matirx_hash[metabolite].setdefault(idtemp, 0.0)
        #
        # Obtain information for stoichiometry matrix construction
        #
        for id in self.reaction_ids:
            #remove spaces in reaction
            reaction = self.reactions[id]['stoichiometry']
            reaction.replace(" ","")
            substrate, product = reaction.split("-->")
            #
            # Substrate
            #
            for metabolite in substrate.split("+"):
                #
                # Check coefficient: {0.5}Pyr
                #
                if len(metabolite.split("}")) == 1:
                    metabolite_name = metabolite
                    stnum = 1.0
                elif len(metabolite.split("}")) == 2:
                    stnum, metabolite_name = metabolite.split("}")
                    stnum = float(stnum.replace("{",""))
                #
                # Ignore if carbon source metabolites
                #
                if self.metabolites[metabolite_name]['carbonsource'] == 'carbonsource':
                    continue
                #
                # Ignore if excreted metabolites
                #
                if self.metabolites[metabolite_name]['excreted'] == 'excreted':
                    continue
                #
                # Ignore if substrates in pseudo reactions
                #
                if self.reactions[id]['type'] == 'pseudo':
                    continue
                #
                stoichiometry_matirx_hash[metabolite_name][id] = stoichiometry_matirx_hash[metabolite_name][id] - float(stnum)
            #
            # Product
            #
            for metabolite in product.split("+"):
                if len(metabolite.split("}")) == 1:
                    metabolite_name = metabolite
                    stnum = 1.0
                elif len(metabolite.split("}")) == 2:
                    stnum, metabolite_name = metabolite.split("}")
                    stnum = stnum.replace("{","")
                if self.metabolites[metabolite_name]['carbonsource'] == 'carbonsource':
                    continue
                if self.metabolites[metabolite_name]['excreted'] == 'excreted':
                    continue
                stoichiometry_matirx_hash[metabolite_name][id] = stoichiometry_matirx_hash[metabolite_name][id]+float(stnum)
        #
        # Initialize stoichiometry matrix
        #
        stoichiometry_matrix = numpy.zeros((len(intermediates), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # List of intermediates
        # This is a row oder in the matrix
        #
        intermediates_sorted = sorted(list(intermediates))
        #
        # Construct a stoichiometry matrix
        #
        for i, intermediate in enumerate (intermediates_sorted):
            self.matrix_row_names.append(intermediate)
            Rm_initial.append(0.0)
            for id in self.reaction_ids:
                stoichiometry_matrix[i][self.reactions[id]["position_in_tmp_r"]] = stoichiometry_matirx_hash[intermediate][id]
                # Callback stoichiometry_matrix information
                if (self.configuration['callbacklevel'] >= 8):

                    if (id in self.reactions[id]['stoichiometry']) and (stoichiometry_matirx_hash[intermediates_sorted[i]][id] == 0):
                        print("Stoichiometry_matrix information:", intermediates_sorted[i], id, self.reactions[id]['stoichiometry'], stoichiometry_matirx_hash[intermediates_sorted[i]][id])
                    if (id not in self.reactions[id]['stoichiometry']) and (stoichiometry_matirx_hash[intermediates_sorted[i]][id] != 0):
                        print("Stoichiometry_matrix information:", intermediates_sorted[i], id, self.reactions[id]['stoichiometry'], stoichiometry_matirx_hash[intermediates_sorted[i]][id])

        #
        # Initialize fixed_flux matrix
        #
        fixed_reactions = [id for id in self.reaction_ids if self.reactions[id]["type"] == "fixed"]
        #fixed_reactions.extend([id for id in self.reaction_ids if self.reactions[id]["type"] == "pseudo"])
        fixedflux_matrix = numpy.zeros((len(fixed_reactions), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # List of fixed reactions
        # This is a row oder in the matrix
        #
        fixed_sorted = list(fixed_reactions)
        for i, ids in enumerate(fixed_sorted):
            self.matrix_row_names.append(ids)
            self.reactions[ids]["position_in_Rm_initial"] = len(Rm_initial)
            Rm_initial.append(self.reactions[ids]['value'])
            fixedflux_matrix[i][self.reactions[ids]["position_in_tmp_r"]] = 1.0

        #
        # Combine Stoichiometry matrix and fixed_flux matrix
        #
        # [ Stoichiometry matrix ]
        # [ fixed_flux matrix    ]
        #
        matrix_st_ff = numpy.vstack((stoichiometry_matrix, fixedflux_matrix))
        #
        #
        # Initialize fixed_met matrix
        #
        fixed_metabolites = [id for id in self.metabolite_ids if self.metabolites[id]["type"] == "fixed"]
        fixedmetabolite_matrix = numpy.zeros((len(fixed_metabolites), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # List of fixed reactions
        # This is a row oder in the matrix
        #
        for i, ids in enumerate(list(fixed_metabolites)):
            self.matrix_row_names.append(ids)
            self.metabolites[ids]["position_in_Rm_initial"] = len(Rm_initial)
            Rm_initial.append(self.metabolites[ids]['value'])
            fixedmetabolite_matrix[i][self.metabolites[ids]["position_in_tmp_r"]] = 1.0
        #
        # Combine Stoichiometry matrix and fixed_flux matrix
        #
        # [ Stoichiometry matrix ]
        # [ fixedflux matrix    ]
        # [ fixedmetabolite matrix    ]
        matrix_st_ff_fm = numpy.vstack((matrix_st_ff, fixedmetabolite_matrix))


        #
        # Initialize fixed_reversible matrix
        #
        fixed_reversible = [id for id in self.reversible_ids if self.reversible[id]["type"] == "fixed"]
        fixedreversible_matrix = numpy.zeros((len(fixed_reversible), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # List of fixed reactions
        # This is a row oder in the matrix
        #
        for i, ids in enumerate(list(fixed_reversible)):
            self.matrix_row_names.append(ids)
            self.reversible[ids]["position_in_Rm_initial"] = len(Rm_initial)
            Rm_initial.append(self.reversible[ids]['value'])

            fixedreversible_matrix[i][self.reversible[ids]["position_in_tmp_r"]] = 1.0


        #
        # Combine Stoichiometry matrix and fixed_flux matrix
        #
        # [ Stoichiometry matrix ]
        # [ fixedflux matrix    ]
        # [ fixedmetabolite matrix    ]
        # [ fixedreversible matrix    ]
        matrix_st_ff_fm_fr = numpy.vstack((matrix_st_ff_fm, fixedreversible_matrix))

        #
        # constraints for reversible reaction
        #
        reversible_matrix = numpy.zeros((len(self.reversible_ids), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # List of fixed reactions
        # This is a row oder in the matrix
        #
        for i, ids in enumerate(list(self.reversible_ids)):
            self.matrix_row_names.append(ids)
            self.reversible[ids]["position_in_Rm_initial"] = len(Rm_initial)
            Rm_initial.append(0.0)
            forward_ids = self.reversible[ids]["forward"]
            reverse_ids = self.reversible[ids]["reverse"]
            reversible_matrix[i][self.reversible[ids]["position_in_tmp_r"]] = -1.0
            #reversible_matrix[i][self.reactions[forward_id]["position_in_tmp_r"]] = 1.0
            #reversible_matrix[i][self.reactions[reverse_id]["position_in_tmp_r"]] = -1.0
            #
            #
            # Added for
            #
            for forward_id in forward_ids.split("+"):
                if forward_id in self.reactions:
                    reversible_matrix[i][self.reactions[forward_id]["position_in_tmp_r"]] = 1.0
            for reverse_id in reverse_ids.split("+"):
                if reverse_id in self.reactions:
                    reversible_matrix[i][self.reactions[reverse_id]["position_in_tmp_r"]] = -1.0

        # [ Stoichiometry matrix ]
        # [ fixedflux matrix    ]
        # [ fixedmetabolite matrix    ]
        # [ fixedreversible matrix    ]
        # [ reversible matrix    ]
        matrix_st_ff_fm_fr_rev = numpy.vstack((matrix_st_ff_fm_fr, reversible_matrix))

        #
        # Set independent fluxes from a reduced row eschron form of matrix_before
        # Selection of independent fluxes depend on the order of the reactions
        #

        tmp_rref = numpy.array(rref(matrix_st_ff_fm_fr_rev.tolist()))
        if (self.configuration['callbacklevel'] >= 8):
            print("stoichiometry_matrix", stoichiometry_matrix)
            print("fixedflux_matrix", fixedflux_matrix)
            print("fixedmetabolite_matrix",fixedmetabolite_matrix)
            print("fixedreversible_matrix", fixedreversible_matrix)
            print("reversible_matrix", reversible_matrix)
            print("matrix_st_ff_fm_fr_rev", matrix_st_ff_fm_fr_rev)
            print("tmp_rref",tmp_rref)
        independence = numpy.ones(tmp_rref.shape[1])
        for i in range(tmp_rref.shape[0]):
            nz = sorted(numpy.nonzero(tmp_rref[i]))[0]
            independence[nz[0]] = 0.0
        independent_flux = numpy.nonzero(independence)[0]
        #
        #independent flux matrix
        #
        independentflux_matrix = numpy.zeros((len(independent_flux), len(self.reaction_ids)+len(self.metabolite_ids)+len(self.reversible_ids)))
        #
        # independent flux
        #
        self.vector["independent_flux"] = [""] * len(independent_flux)

        for i, position in enumerate(independent_flux):
            #
            # independend flux out side of flux value
            #
            for id in self.reaction_ids:
                if position == self.reactions[id]["position_in_tmp_r"]:
                    self.matrix_row_names.append(id)
                    self.reactions[id]["position_in_Rm_initial"] = len(Rm_initial)
                    self.vector["independent_flux"][i] = ("reaction", id)

            for id in self.metabolite_ids:
                if position == self.metabolites[id]["position_in_tmp_r"]:
                    self.matrix_row_names.append(id)
                    self.metabolites[id]["position_in_Rm_initial"] = len(Rm_initial)
                    self.vector["independent_flux"][i] = ("metabolite", id)

            for id in self.reversible_ids:
                if position == self.reversible[id]["position_in_tmp_r"]:
                    self.matrix_row_names.append(id)
                    self.reversible[id]["position_in_Rm_initial"] = len(Rm_initial)
                    self.vector["independent_flux"][i] = ("reversible", id)

            Rm_initial.append(0.0)
            independentflux_matrix[i][independent_flux[i]] = 1
        #
        # Construct matrix
        #
        # [ Stoichiometry matrix ]
        # [ fixedflux matrix    ]
        # [ fixedmetabolite matrix ]
        # [ fixedreversible matrix    ]
        # [ reversible matrix    ]
        # [ indepenent matrix    ]
        self.matrix = numpy.vstack((matrix_st_ff_fm_fr_rev, independentflux_matrix))
        #
        #
        # Set numbers
        self.numbers = {}
        self.numbers["stoichiometric_start"] = 0
        self.numbers["stoichiometric_end"] = stoichiometry_matrix.shape[0]
        self.numbers["fixed_reaction_start"] = self.numbers["stoichiometric_end"]
        self.numbers["fixed_reaction_end"] = self.numbers["fixed_reaction_start"] + fixedflux_matrix.shape[0]
        self.numbers["fixed_metabolite_start"] = self.numbers["fixed_reaction_end"]
        self.numbers["fixed_metabolite_end"] = self.numbers["fixed_metabolite_start"] + fixedmetabolite_matrix.shape[0]
        self.numbers["fixed_reversible_start"] = self.numbers["fixed_metabolite_end"]
        self.numbers["fixed_reversible_end"] = self.numbers["fixed_reversible_start"] + fixedreversible_matrix.shape[0]
        self.numbers["reversible_start"] = self.numbers["fixed_reversible_end"]
        self.numbers["reversible_end"] = self.numbers["reversible_start"] + reversible_matrix.shape[0]
        self.numbers["independent_start"] = self.numbers["reversible_end"]
        self.numbers["independent_end"] = self.numbers["independent_start"] + independentflux_matrix.shape[0]
        self.numbers["independent_number"] = independentflux_matrix.shape[0]
        self.numbers["total_number"] = self.numbers["independent_end"] * 1
        #
        # for diffmdv
        #
        self.numbers["reac_met_number"] = len(self.reaction_ids)+len(self.metabolite_ids)
        #
        # Rm_initial
        #
        self.vector["Rm_initial"] = numpy.array(Rm_initial)
        self.vector["independent_lb"] = [self.vector["lb"][v] for v in independent_flux]
        self.vector["independent_ub"] = [self.vector["ub"][v] for v in independent_flux]
        #
        # independent flux
        #
        self.vector["independent_flux_position"] = list(independent_flux)
        #
        # Construct inversed matrix
        #
        self.matrixinv = numpy.linalg.inv(self.matrix)
        #me
        non_fixed_reactions = []
        for i, id in enumerate(self.reaction_ids):
            if self.reactions[id]['type'] != 'fixed':
                non_fixed_reactions.append(i)
        self.vector["non_fixed_reactions"] = list(independent_flux)
        #
        # Reaction depencency
        # This information is used in the Kajihata's LM
        #
        self.r_depended = []
        for i in range(tmp_rref.shape[1]):
            temp = []
            for i, position in enumerate(independent_flux):
                #
                # independend flux out side of flux value
                #
                group, id = self.vector["independent_flux"][i]
                if group == "reaction":
                    pos_in_Rm = self.reactions[id]["position_in_Rm_initial"]
                elif group == "metabolite":
                    pos_in_Rm = self.metabolites[id]["position_in_Rm_initial"]
                elif group == "reversible":
                    pos_in_Rm = self.reversible[id]["position_in_Rm_initial"]
                if self.matrixinv[i][pos_in_Rm] != 0:
                    temp.append(pos_in_Rm - self.numbers["independent_start"])
            self.r_depended.append(list(temp))
        return True


    def reconstruct(self, mode = "no"):
        """
        Reconstruct the metabolic model to generate new stoichiometry matrix and a function to
        calculate MDV data. The metabolic model have to be reconstructed when total number of metabolic
        reactions were changed.
        function_text = model.reconstruct()


        Parameters
        ----------
        mode: experimental "no"


        Examples
        --------
        >>> function_text = model.reconstruct()
        >>> exec(function_text)
        >>> fitting_flux(calmdv, flux)
        >>> f = open('calmdv.py', 'w')
        >>> f.write(function_text)
        >>> f.close()
        >>> from calmdv import calmdv

        Returns
        ----------
        test : "calmdv" and "diffmdv" str data describing a python script to determine MDV from metabolic flux.

        See Also
        --------
        update

        """

        self.update()
        calmdv_text, ccalmdv_text = self.generate_calmdv(mode)
        self.calmdv_text = calmdv_text
        #print(calmdv_text)
        #
        # local_dic was added for Python 3
        #
        locals_dic = locals()
        exec(calmdv_text, globals(), locals_dic)
        calmdv = locals_dic["calmdv"]
        diffmdv = locals_dic["diffmdv"]

        self.func = {}

        self.func["calmdv"] = calmdv
        self.func["diffmdv"] = diffmdv
        return calmdv_text

    def generate_calmdv(self,mode = 'normal'):
        """
        Constructer of the python function to calculate MDV of target fragments from a metabolic flux distribution.
        This function is performed in reconstruct()

        Parameters
        ----------
        mode: optimize EMU network reduction (experimental)

        Returns
        --------
        thres: a threshold level of RSS for searching confidence interval

        Examples
        --------
        >>> calmdv = model.generate_calmdv()

        """
        def generate_iterative_list(list):
            result = []
            temp = []
            result = step(0, list, result,temp[:])
            return(result)
        def step(number, list, result, temp):
            for i in range(list[number] + 1):
                tempp = temp[:]
                tempp.append(i)
                if (number == len(list)-1):
                    result.append(tempp)
                else:
                    result = step(number+1, list, result, tempp[:])
            return(result)
        def decompose_EMU(string):
            compound, atom_string = string.split('_')
            #atoms = list(atom_string)
            atoms = atom_string.split(':')
            return (compound, atoms)
        def compose_EMU(compound, atoms):
            atomstring = ":".join([str(i) for i in atoms])
            return (compound+"_"+atomstring)


        def generate_parent_EMU(reaction, emu, rev, symmetry, metabolites):

            compound, numbers = decompose_EMU(emu)
            if self.configuration['callbacklevel'] == 7:
                print(emu, compound, numbers)
            reaction.replace(' ', '')#remode space
            parent_atoms = convert_reaction_to_EMU(reaction)
            # when rec == 1
            if rev == 1:
                parent_atoms = [parent_atoms[int(x)] for x in sorted(range(len(parent_atoms)), reverse = True)]
            emu_candidates = [parent_atoms[int(x) - 1] for x in numbers]
            emu_dict = {};
            emu_list = [];

            for emu_atom in emu_candidates:
                compound, number = emu_atom.split('_')
                if compound in emu_dict:
                    emu_dict[compound].append(number)
                else:
                    emu_dict[compound] = [number]
            for compound_tmp in sorted(emu_dict.keys()):#Generate EMU symbols for each substrate
                #
                # [1]DHAP => DHAP
                #
                num, compound = compound_tmp.split("]")
                if self.configuration['callbacklevel'] == 7:
                    print(compound, sorted(emu_dict[compound_tmp]))
                emu = compose_EMU(compound, sorted(emu_dict[compound_tmp]))
                #for number in sorted(emu_dict[compound_tmp]):
                #    emu += number
                if compound in symmetry:# Symmetry metabolite
                    emu_rev = compound + "_";
                    tmp = [int(metabolites[compound]['C_number']) + 1 - int(x) for x in emu_dict[compound_tmp]]
                    for number in sorted(tmp):
                        emu_rev += str(number)
                    emu_rev = compose_EMU(compound, sorted(sorted(tmp)))
                    emu = sorted([emu, emu_rev])[0]
                emu_list.append(emu)

            return (emu_list);
        def convert_reaction_to_EMU(string):
            string.replace(' ', '')
            result = []
            for emu in string.split('+'):
                result = result + convert_one_symbol_to_EMU(emu)
            return(result)
        def convert_one_symbol_to_EMU(string):
            string.replace(' ', '')
            compound, atom_strings = string.split('_')
            result = []
            for atom in atom_strings.split(":"):
                result.append(compound + "_" + atom)
            return (result)
        def size_of_EMU(string):
            string.replace(' ', '')
            numbers = re.findall("[_|:]([1-9]+)", string)
            return(len(numbers))
            #return ( int(sum ([len(x) for x in numbers])));
        def compound_of_EMU(string):
            string.replace(' ', '')
            compound, atoms = string.split('_')
            return (compound)

        def permute(A):
            #A=[2,3]
            #B=[[0,1],[0,1,2]]
            B=[range(i) for i in A]
            n=len(B)
            if n==2:
                return list(itertools.product(B[0],B[1]))

            elif n==3:
                return  list(itertools.product(B[0],B[1],B[2]))

            elif n==4:
                return  list(itertools.product(B[0],B[1],B[2],B[3]))

            elif n==5:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4]))

            elif n==6:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5]))

            elif n==7:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6]))

            elif n==8:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7]))

            elif n==9:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8]))

            elif n==10:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9]))
            elif n==11:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10]))

            elif n==12:
                return  list(itertools.product(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11]))

            else:
                print("Error! number of substrate are too much")

        #
        # String of constructed function
        #
        string = ""
        #
        # String of constructed function
        #
        cython_string = ""
        #
        # Maximal carbon number (EMU size) of metabolites in the model
        #
        emu_size_max = max([self.metabolites[x]['C_number'] for x in self.metabolites.keys()])
        #
        # Preparation of "pathway"
        #
        pathway = []
        #
        #FBP-->DHAP+GAP ABCDEF-->CBA+DEF
        #
        for id in self.reaction_ids:
            reaction = self.reactions[id]['reaction']
            atommap = self.reactions[id]['atommap']
            #
            #
            #FBP-->DHAP+GAP
            #ABCDEF-->CBA+DEF
            #
            # ignore when atom map information is unavilable
            #
            if atommap == "nd":
                continue
            #
            # ignore is atom map string does not indlude "-->"
            #
            if not '-->' in atommap:
                continue
            #
            # Discard space and "+"
            #
            reaction = reaction.replace(" ", "")
            atommap = atommap.replace(" ", "")
            atommap = atommap.replace("+", "")
            #
            #FBP-->DHAP+GAP
            #ABCDEF-->CBADEF
            #
            # reaction is splitted by '-->'
            # substrate: FBP
            # product: DHAP+GAP
            substrates,products = reaction.split("-->")
            #
            # atom map is splitted by '-->'
            # substrate: ABCDEF
            # product: CBADEF
            atom_substrate, atom_product = atommap.split("-->")
            #
            # List and dictionary for storing some information
            #
            substrate_atoms = []
            substrate_atom_hash = {}
            product_atoms = []
            product_atom_hash = {}
            #
            # FBP -> FBP
            #
            for n, substrate in enumerate(substrates.split("+")):
                #
                # ignore if extreted metabolite
                #
                if self.metabolites[substrate]['excreted'] == 'excreted':
                    continue
                #
                # separate metabolite to each carbon
                #
                for i in range(int(self.metabolites[substrate]['C_number'])):
                    substrate_atoms.append("["+str(n)+"]"+substrate + "_" + str(i + 1))
                # FBP ->
                #substrate_atoms = [[1]FBP_1,[1]FBP_2,[1]FBP_3,[1]FBP_4,[1]FBP_5,[1]FBP_6]
                #
            #
            #
            for i, substrate_atom  in enumerate (substrate_atoms):
                substrate_atom_hash[atom_substrate[i]] = substrate_atom
                #substrate_atom_hash[A]:[1]FBP_1
                #substrate_atom_hash[B]:[1]FBP_2
                #substrate_atom_hash[C]:[1]FBP_3
                #substrate_atom_hash[D]:[1]FBP_4
                #substrate_atom_hash[E]:[1]FBP_5
                #substrate_atom_hash[F]:[1]FBP_6
            #
            #DHAP+GAP
            #
            for n, product in enumerate(products.split("+")):
                #
                # This is not a good idea
                #
                if self.metabolites[product]['excreted'] == 'excreted':
                    continue
                for i in range(int(self.metabolites[product]['C_number'])):

                    product_atoms.append("["+str(n)+"]"+product + "_" + str(i + 1))
            #
            # product_atoms_atoms = [[1]DHAP_1,[1]DHAP_2,[1]DHAP_3,[2]GAP_1,[2]GAP_2,[2]GAP_3]
            #
            for i in range(len(product_atoms)):
                product_atom_hash[atom_product[i]] = product_atoms[i]
                #product_atom_hash[A]:[1]DHAP_3
                #product_atom_hash[B]:[1]DHAP_2
                #product_atom_hash[C]:[1]DHAP_1
                #product_atom_hash[D]:[2]GAP_1
                #product_atom_hash[E]:[2]GAP_2
                #product_atom_hash[F]:[2]GAP_3
            #
            # atom mapping from products to substrate
            #
            p_atom_to_s_atom = {};
            for atom in product_atom_hash:
                p_atom_to_s_atom[product_atom_hash[atom]] = substrate_atom_hash[atom]
            #p_atom_to_s_atom
            #[1]DHAP_3:[1]FBP_1
            #[1]DHAP_2:[1]FBP_2
            #[1]DHAP_1:[1]FBP_3
            #[2]GAP_1:[1]FBP_4
            #[2]GAP_2:[1]FBP_5
            #[2]GAP_3:[1]FBP_6

            # for each product
            # DHAP+GAP
            #
            for i,product in enumerate(products.split("+")):
                #
                # Ignore exceted metabolites
                #
                if self.metabolites[product]['excreted'] == 'excreted':
                    continue
                #
                # Symmetry product
                #
                if product in self.symmetry:
                    id_for_output = id + "*"
                else:
                    id_for_output = id
                current_compound = ""
                substrate_emu = ""
                #
                #product:DHAP
                #
                product_ordered = "["+str(i)+"]"+product
                #
                substrate_list = [x for x in (sorted(p_atom_to_s_atom.keys())) if product_ordered in x]
                #
                #substrate_list = [DHAP_1, DHAP_2, DHAP_3]
                #
                for emu_atom in sorted(substrate_list):
                    if current_compound == p_atom_to_s_atom[emu_atom].split("_")[0]:
                        substrate_emu = substrate_emu+":"+p_atom_to_s_atom[emu_atom].split("_")[1]
                        #2nd loop
                        #substrate_emu = p_atom_to_s_atom[[1]DHAP_2] = [1]FBP_32
                        #3rd loop
                        #substrate_emu = p_atom_to_s_atom[[1]DHAP_3] = [1]FBP_321
                    elif substrate_emu == "":
                        #1st loop
                        substrate_emu = p_atom_to_s_atom[emu_atom]
                        #substrate_emu = p_atom_to_s_atom[[1]DHAP_1] = [1]FBP_3
                        current_compound = p_atom_to_s_atom[emu_atom].split("_")[0]
                        #current_compound = [1]FBP
                    else:
                        #Compound is changed in the atom list
                        substrate_emu = substrate_emu+"+"+p_atom_to_s_atom[emu_atom]
                        current_compound = p_atom_to_s_atom[emu_atom].split("_")[0]
                # Recored as "v1 DHAP [1]FBP_321"
                pathway.append(id_for_output+" "+product+" "+substrate_emu)
                #
                # [1]DHAP: [1]FBP_321
                # [2]GAP: [1]FBP_123
                if self.configuration['callbacklevel'] == 7:
                    print("Define product-> substrate_carbon relationship", id_for_output+" "+product+" "+substrate_emu)


        #Initialize dictionaries
        matrix = {};
        hit_emuset = set();
        previous_hit_emuset = set();
        hit_emu_reactionset = set();

        emu_intermediate = {};
        emu_source = {};

        self.emu_order_in_X = [];
        self.emu_order_in_y = [];

        # Search metabolic network from target compounds
        # ex. Leu57	gcms	AcCOAmit_12 + PYRmit_23 + PYRmit_23	no
        #
        for target_fragment in self.target_fragments:
            #
            # Ignore non used fragment
            #
            if self.target_fragments[target_fragment]['use'] != 'use':
                continue
            #
            # Leu57:AcCOAmit_1:2 + PYRmit_2:3 + PYRmit_2:3
            #
            target_emus = self.target_fragments[target_fragment]['atommap']
            #
            # AcCOAmit_1:2 + PYRmit_2:3 + PYRmit_2:3 => [AcCOAmit_1:2, PYRmit_2:3, PYRmit_2:3]
            #
            emus = target_emus.replace(" ", "").split("+")
            #
            # submit EMUs to hit_emuset
            #
            hit_emuset.update(emus)
            #
            # While all emus are checked,
            #
            while True:
                #
                # for each EMU in hit_emuset
                #
                for product_emu in hit_emuset.copy():
                    #
                    # ignore if the EMUS has already been cheched
                    #
                    if product_emu in previous_hit_emuset:
                        continue
                    #
                    # get names and corbon list of EMS
                    # decompose_EMU('FBP_3:4:2') => ('FBP', ['3', '4', '2'])
                    product, numbers = decompose_EMU(product_emu)
                    if self.configuration['callbacklevel'] == 7:
                        print('Checking:', product_emu, product, numbers)
                    #
                    # for each pathway data
                    #
                    for data in pathway:
                        #
                        # ignore when incomplete data
                        #
                        if len(data.split(' ')) < 3:
                            continue
                        reaction_id, r_product, reaction = data.split(' ');
                        #
                        # if a reaction producing indentical product is found
                        #
                        if product == r_product:
                            if self.configuration['callbacklevel'] == 7:
                                print('|-Hitted:', reaction_id, r_product, reaction)
                            #
                            # get substrate EMU(s)
                            #
                            #
                            # when symmetry metabolites
                            #
                            if '*' in reaction_id:
                                # When symmetry product
                                substrate_emus = generate_parent_EMU(reaction, product_emu, 0, self.symmetry, self.metabolites)
                                substrate_emus_reaction = '+'.join(substrate_emus)
                                hit_emu_reactionset.update([substrate_emus_reaction])
                                # Add substrate to hit_emu
                                for substrate_emu in substrate_emus:
                                    hit_emuset.update([substrate_emu])
                                if self.configuration['callbacklevel'] == 7:
                                    print('  |-Found the forward relationship for symmetry metabolile', reaction, r_product, product_emu, reaction_id, substrate_emus)
                                # Reversed
                                substrate_emus = generate_parent_EMU(reaction, product_emu, 1, self.symmetry, self.metabolites)
                                substrate_emus_reaction_rev = '+'.join(substrate_emus)
                                hit_emu_reactionset.update([substrate_emus_reaction_rev])
                                # Add substrate to hit_emu
                                for substrate_emu in substrate_emus:
                                    hit_emuset.update([substrate_emu])
                                if self.configuration['callbacklevel'] == 7:
                                    print('  |-Found the reverse relationship for symmetry metabolile', reaction, r_product, product_emu, reaction_id, substrate_emus)
                            #
                            # MS/MS
                            #
                            elif r_product in self.target_fragments and self.target_fragments[r_product]['type'] == "msms":
                                #decode
                                substrate_emus = reaction.split("+")
                                substrate_emus_reaction = '+'.join(substrate_emus)
                                #
                                #
                                hit_emu_reactionset.update([substrate_emus_reaction])
                                # Add substrate to hit_emu
                                for substrate_emu in substrate_emus:
                                    hit_emuset.update([substrate_emu])
                            #
                            # non synmetry metabolite
                            #
                            else:
                                # When non-symmetry product

                                substrate_emus = generate_parent_EMU(reaction, product_emu, 0, self.symmetry, self.metabolites)
                                if self.configuration['callbacklevel'] == 7:
                                    print('  |-Found the relationship for non symmetry metabolile', reaction, product_emu, substrate_emus)
                                substrate_emus_reaction = '+'.join(substrate_emus)
                                hit_emu_reactionset.update([substrate_emus_reaction])
                                # Add substrate to hit_emu
                                for substrate_emu in substrate_emus:
                                    hit_emuset.update([substrate_emu])

                            # Add product to hit_emu_reactionset
                            hit_emu_reactionset.update([product_emu])
                            #
                            # Add data to matrix
                            #
                            if product_emu in matrix:
                                # When product_emu exists in matrix

                                if substrate_emus_reaction in matrix[product_emu]:
                                    # When substrate_emus_reaction is in product_emu
                                    # Add id
                                    matrix[product_emu][substrate_emus_reaction].append(reaction_id)
                                    if self.configuration['callbacklevel'] == 7:
                                        print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction, reaction_id)

                                else:
                                    # Prepare new data
                                    matrix[product_emu].setdefault(substrate_emus_reaction, [reaction_id])
                                    if self.configuration['callbacklevel'] == 7:
                                        print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction, reaction_id )
                            else:
                                # When product_emu is NOT in matrix
                                # Prepare new data
                                matrix[product_emu] = {substrate_emus_reaction:[reaction_id]}
                                if self.configuration['callbacklevel'] == 7:
                                    print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction, reaction_id)
                            # For orthogonal data
                            if product_emu in matrix[product_emu]:
                                matrix[product_emu][product_emu].append(reaction_id)
                                if self.configuration['callbacklevel'] == 7:
                                    print('    |-Added the diagonal relationship to matrix', product_emu, reaction_id)
                            else:
                                matrix[product_emu].setdefault(product_emu, [reaction_id])
                                if self.configuration['callbacklevel'] == 7:
                                    print('    |-Added the diagonal relationship to matrix', product_emu, reaction_id)

                            if '*' in reaction_id:
                                # Add data to matrix
                                if product_emu in matrix:
                                    if substrate_emus_reaction_rev in matrix[product_emu]:
                                        # When substrate_emus_reaction is in product_emu
                                        # Add id
                                        matrix[product_emu][substrate_emus_reaction_rev].append(reaction_id)
                                        if self.configuration['callbacklevel'] == 7:
                                            print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction_rev, reaction_id)
                                    else:
                                        # Prepare new data
                                        matrix[product_emu].setdefault(substrate_emus_reaction_rev, [reaction_id])
                                        if self.configuration['callbacklevel'] == 7:
                                            print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction_rev, reaction_id)
                                else:
                                    # When product_emu exists is NOT in matrix
                                    # Prepare new data
                                    matrix[product_emu] = {substrate_emus_reaction_rev:[reaction_id]}
                                    if self.configuration['callbacklevel'] == 7:
                                        print('    |-Added the relationship to matrix', product_emu, substrate_emus_reaction_rev, reaction_id)
                                if product_emu in matrix[product_emu]:
                                    matrix[product_emu][product_emu].append(reaction_id)
                                    if self.configuration['callbacklevel'] == 7:
                                        print('    |-Added the diagonal relationship to matrix', product_emu, reaction_id)
                                else:
                                    matrix[product_emu].setdefault(product_emu, [reaction_id])
                                    if self.configuration['callbacklevel'] == 7:
                                        print('    |-Added the diagonal relationship to matrix', product_emu, reaction_id)
                    # Update previous_hit_emuset by newly found product_emu
                    previous_hit_emuset.update([product_emu])
                # When no new product_emu is found
                if previous_hit_emuset == hit_emuset:
                    break


        #
        # Import 'deepcopy'
        #
        from copy import deepcopy
        target_fragments = deepcopy(self.target_fragments)

        if mode == "optimize":
            #
            # This is an experimental function to reduce matrix
            # There are several bugs in the function
            #
            emus_in_fragment = set()
            for target_fragment in target_fragments:
                # Ignore
                if target_fragments[target_fragment]['use'] != 'use':
                    continue
                # Leu57:AcCOAmit_12 + PYRmit_23 + PYRmit_23
                target_emus = target_fragments[target_fragment]['atommap']
                # AcCOAmit_12 + PYRmit_23 + PYRmit_23 => [AcCOAmit_12, PYRmit_23, PYRmit_23]
                emus = target_emus.replace(" ", "").split("+")
                emus_in_fragment.update(emus)

            #
            # Dic to find precursor of EMUs of interest
            # r16: IsoCit_3 => aKG_2 as {'aKG_2': 'IsoCit_3'}
            #
            emus_in_matrix = {}
            for product_emu in matrix:
                if len(matrix[product_emu][product_emu]) != 1: continue #Only EMUs generated from only one precursor EMU
                if product_emu in emus_in_fragment:continue # ignore direct precursor
                if product_emu.split('_')[0] in self.symmetry: continue #Ignore EMUs of symmetry metaboiltes
                for substrate_emu in matrix[product_emu]:
                    if substrate_emu != product_emu:
                        emus_in_matrix[product_emu] = str(substrate_emu)
                        # r16: IsoCit_3 => aKG_2 as {'aKG_2': 'IsoCit_3'}
            #
            #for each candidate EMUs for deletion
            #
            emus_to_remove_list = emus_in_matrix.keys()
            for emus_to_remove in emus_to_remove_list:

                matrix_temp = deepcopy(matrix) #generate temporary matrix
                if self.configuration['callbacklevel'] == 7:
                    print('searchng', emus_to_remove)
                counter = 0
                for product_emu in matrix:
                    if product_emu==emus_to_remove: continue
                    for substrate_emu in matrix[product_emu]:
                        if substrate_emu==emus_to_remove:
                            #
                            # if new partial reaction still exists.
                            #
                            if emus_in_matrix[substrate_emu] in matrix_temp[product_emu]:
                                matrix_temp[product_emu][emus_in_matrix[substrate_emu]].extend(matrix[product_emu][substrate_emu][:])
                                if self.configuration['callbacklevel'] == 7:
                                    print("removed + added", product_emu, substrate_emu, emus_to_remove, emus_in_matrix[substrate_emu], matrix[product_emu][substrate_emu][:], matrix_temp[product_emu][emus_in_matrix[substrate_emu]])
                            else:
                                matrix_temp[product_emu][emus_in_matrix[substrate_emu]] = matrix[product_emu][substrate_emu][:]
                                if self.configuration['callbacklevel'] == 7:
                                    print("removed ", product_emu, substrate_emu, emus_to_remove, emus_in_matrix[substrate_emu], matrix[product_emu][substrate_emu][:], matrix_temp[product_emu][emus_in_matrix[substrate_emu]] )
                            matrix_temp[product_emu].pop(substrate_emu)
                            counter = counter + 1
                        if '+' in substrate_emu:
                            substrate_emu_temp = []
                            for emu in substrate_emu.split('+'):
                                if emu == emus_to_remove:
                                    substrate_emu_temp.append(emus_in_matrix[emu])
                                else:
                                    substrate_emu_temp.append(emu)
                            substrate_emu_temp = "+".join(substrate_emu_temp)
                            if substrate_emu_temp != substrate_emu:
                                #
                                # if new partial reaction still exists.
                                #
                                if substrate_emu_temp in matrix_temp[product_emu]:
                                    matrix_temp[product_emu][substrate_emu_temp].extend(matrix[product_emu][substrate_emu][:])
                                    if self.configuration['callbacklevel'] == 7:
                                        print("removed + added +", product_emu, substrate_emu, substrate_emu_temp, matrix[product_emu][substrate_emu][:], matrix_temp[product_emu][substrate_emu_temp])
                                else:
                                    matrix_temp[product_emu][substrate_emu_temp] = matrix[product_emu][substrate_emu][:]
                                    if self.configuration['callbacklevel'] == 7:
                                        print("removed +", product_emu, substrate_emu, substrate_emu_temp, matrix[product_emu][substrate_emu][:], matrix_temp[product_emu][substrate_emu_temp])
                                matrix_temp[product_emu].pop(substrate_emu)
                                hit_emu_reactionset.remove(substrate_emu)
                                hit_emu_reactionset.update([substrate_emu_temp])
                                counter = counter + 1
                if counter > 0:
                    matrix_temp.pop(emus_to_remove)
                    #
                    # delete from emus_in_matrix
                    #
                    for product_emu in emus_in_matrix:
                        # replace precursor data
                        if emus_in_matrix[product_emu] == emus_to_remove:
                            emus_in_matrix[product_emu] = str(emus_in_matrix[emus_to_remove])
                    # discord product data
                    emus_in_matrix.pop(emus_to_remove)
                    hit_emu_reactionset.remove(emus_to_remove)
                    if self.configuration['callbacklevel'] == 7:
                        print(emus_to_remove, "was removed")
                matrix = deepcopy(matrix_temp)



        #
        #start making function
        #
        string += "import numpy as numpy\n"
        string += "import scipy as scipy\n"
        string += "\n"
        string += "def "
        string += "calmdv"
        string += "(r, target_emu_list, mdv_carbon_sources):\n"
        for i in range(len(self.reaction_ids)):
            string += "\t" + self.reaction_ids[i] + " = r[" + str(i) + "]\n"
            #string += "\t" + self.reaction_ids[i] + " = r[" + str(i) + "]\n"
        string += "\temu_list = {}\n"
        string += "\temu_list_isotope_corrected = {}\n"
        string += "\tX_list = []\n"
        #
        #print matrix
        # AX = BY

        for emu_size in range(1, int(emu_size_max)+1):
            #
            # List up all EMUs of intermediated in the X of emu_size
            # from hit_emu_reactionset, side of emu is equal to emu_size, only one metabolite and not carbon source.
            emu_intermediate_of_this_layer = [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and (not ('+' in x))
                                                                                and not (compound_of_EMU(x) in self.carbon_source)]
            emu_intermediate_of_this_layer.sort()

            if len(emu_intermediate_of_this_layer) == 0:
                continue
            #
            emu_sourse_of_this_layer = []
            # from hit_emu_reactionset, side of emu is equal to emu_size, more than one metabolite
            emu_sourse_of_this_layer += [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and ('+' in x)]
            # from hit_emu_reactionset, side of emu is equal to emu_size, one metabolite and carbon source
            emu_sourse_of_this_layer += [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and (not ('+' in x))
                                                                                and (compound_of_EMU(x) in self.carbon_source)]

            emu_sourse_of_this_layer.sort()
            #
            # corresponding row of each emu in matrix X
            # Store EMU name <=> row number relationship
            #
            for i in range(len(emu_intermediate_of_this_layer)):
                if emu_size in emu_intermediate:
                    emu_intermediate[emu_size].setdefault(emu_intermediate_of_this_layer[i],i)
                else:
                    emu_intermediate[emu_size] = {emu_intermediate_of_this_layer[i]:i}
            # Store position in linerized matrix X <=> emu
            for i, emu in enumerate(emu_intermediate_of_this_layer):
                for num in range(int(emu_size)+1):
                    self.emu_order_in_X.append((emu, num))

            #
            for i in range(len(emu_intermediate_of_this_layer)):
                if emu_size in emu_intermediate:
                    emu_intermediate[emu_size].setdefault(emu_intermediate_of_this_layer[i],i)
                else:
                    emu_intermediate[emu_size] = {emu_intermediate_of_this_layer[i]:i}

            #Output A
            # number of line
            size_i = len(emu_intermediate_of_this_layer)
            # number of lows
            size_s = len(emu_sourse_of_this_layer)
            string += "\tA_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_i)+",))\n"
            for product in emu_intermediate_of_this_layer:
                string += "\t#"+product+"\n"
                for substrate in emu_intermediate_of_this_layer:
                    cell = []
                    if substrate in matrix[product]:
                        if substrate == product:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append(reaction_id.replace("*", "") + "/ 2.0")
                                else:
                                    cell.append(reaction_id)
                            cell_string = "-1.0 * (" + ' + '.join(cell) + ")"
                            string += "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            if self.configuration['callbacklevel'] == 8:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)

                        else:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append(reaction_id)
                            cell_string = ' + '.join(cell)
                            string += "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            if self.configuration['callbacklevel'] == 8:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)

            string += "\tA_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"

        # Output B
            string += "\tB_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_s)+",))\n"
            for product in emu_intermediate_of_this_layer:
                for substrate in emu_sourse_of_this_layer:
                    cell = []
                    if self.configuration['callbacklevel'] == 8:
                        print(emu_size, product, substrate, matrix[product])
                    if substrate in matrix[product]:
                        if substrate == product:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("-0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append("-1.0 * " + reaction_id)
                        else:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("-0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append("-1.0 * " + reaction_id)
                    if len(cell) != 0:
                        cell_string = ' + '.join(cell)
                        if self.configuration['callbacklevel'] == 8:
                            print(product,'\t',substrate,"\tB_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_s + emu_sourse_of_this_layer.index(substrate))+"] = " + cell_string)

                        string += "\tB_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_s + emu_sourse_of_this_layer.index(substrate))+"] = " + cell_string + "\n"

            string += "\tB_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_s)+"))\n"
        # Output Y
            string += "\tY_" + str(emu_size) + " =  numpy.zeros(("+str(len(emu_sourse_of_this_layer) * (emu_size + 1))+",))\n"
            for substrate in emu_sourse_of_this_layer:
                row = []
                #print substrate
                if '+' in substrate:# for rection with multiple substrates
                    subs_matrix=[]
                    list_emusize_plus1=[]
                    n = len(substrate.split('+'))
                    for subs in substrate.split('+'):   #multiple substrate reaction
                        list_emusize_plus1.append(size_of_EMU(subs)+1)

                        if compound_of_EMU(subs) in self.carbon_source:# Carbon source
                            #sub1_list = generate_carbonsource_MID(sub1)
                            compound = compound_of_EMU(subs)
                            size = size_of_EMU(subs);
                            subs_matrix.append(["mdv_carbon_sources[\"" + subs + "\"][" + str(x) + "]" for x in range(size+1)])

                        else:
                            compound = compound_of_EMU(subs)
                            size = size_of_EMU(subs);
                            subs_matrix.append(["X_" + str(size) + "[" + str((emu_intermediate[int(size)][subs]) * (size+1) + x) + "]" for x in range(size+1)])
                    sum_length=0
                    for subs_list in subs_matrix:
                        sum_length=sum_length +len(subs_list)
                    equation = [[] for x in range(sum_length)];
                    conbinaion_list=permute(list_emusize_plus1)
                    for conbination in conbinaion_list:
                       equation[sum(conbination)].extend([("*").join([str(subs_matrix[i][conbination[i]]) for i in range(len(conbination))])])

                    for numberofisotope in range (len(equation)-n+1):
                        if self.configuration['callbacklevel'] == 8:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]))
                        string += "\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"

                else:
                    compound = compound_of_EMU(substrate)
                    size = size_of_EMU(substrate);
                    for numberofisotope in range(size+1):
                        if self.configuration['callbacklevel'] == 8:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "]")
                        string += "\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "]" + "\n"

            string += "\tY_" + str(emu_size) + ".resize(("+str(size_s)+","+ str(emu_size + 1)+"))\n"
            # Calc X
            string += "\tnonzero_row = numpy.nonzero(numpy.diag(A_" + str(emu_size) + ") < -0.001)[0]\n"
            string += "\tif len(nonzero_row) == A_" + str(emu_size) + ".shape[0]:\n"
            string += "\t\tX_" + str(emu_size) + " = numpy.linalg.solve(A_" + str(emu_size) + ", numpy.dot(B_" + str(emu_size) + ", Y_" + str(emu_size) + "))\n"

            #string += "\t\tX_" + str(emu_size) + " = numpy.dot(numpy.linalg.inv(A_" + str(emu_size) + "), numpy.dot(B_" + str(emu_size) + ", Y_" + str(emu_size) + "))\n"

            string += "\telse:\n"
            string += "\t\teye = numpy.eye(A_" + str(emu_size) + ".shape[0])[nonzero_row,:]\n"
            string += "\t\teyeT = numpy.eye(A_" + str(emu_size) + ".shape[0])[:,nonzero_row]\n"
            #string += "\t\tprint 'yieldee'\n"
            string += "\t\tX_" + str(emu_size) + " = numpy.dot(eyeT, numpy.linalg.solve(numpy.dot(numpy.dot(eye, A_" + str(emu_size) + "), eyeT), numpy.dot(eye, numpy.dot(B_" + str(emu_size) + ", Y_" + str(emu_size) + "))))\n"
            #string += "\tprint X_" + str(emu_size) + "\n"
            string += "\tX_list.extend(list(X_" + str(emu_size)+".ravel()))\n"
            string += "\tX_" + str(emu_size) + " = X_" + str(emu_size) + ".reshape(("+str(size_i * (emu_size + 1) )+",))\n"


        # Calc MV
        for target_fragment in target_fragments:
            #Ignore
            if target_fragments[target_fragment]['use'] != 'use':
                continue
            target_emus = target_fragments[target_fragment]['atommap']
            #msms data
            if target_fragments[target_fragment]['type'] == "msms":
                precursor, neutralloss, product = target_emus.replace(' ','').split('+')
                # get sise information
                size_precursor = size_of_EMU(precursor);
                size_product = size_of_EMU(product);
                size_neutralloss = size_of_EMU(neutralloss);
                product_carbon_number = list(product.split('_')[1]);
                # Mask of product carbon
                product_mask = [];
                for i in (range(size_precursor)):
                    product_mask.append("0")
                for i in (product_carbon_number):
                    product_mask[size_precursor - int(i)] = "1"
                # Number of IDV
                numberofidv = 2 ** size_precursor
                # intitialize mdv_matrix
                mdv_matrix = [];
                for pre in (range(0, size_precursor + 1)):
                    temp = []
                    for pro in (range(0, size_product + 1)):
                        temp.append('')
                    mdv_matrix.append(temp)
                # Check IDV numbers for precursor and product ions
                # Number of MRM series
                mrm_count = 0
                for i in (sorted(range(numberofidv))):
                    #Generate IDV
                    idv="{0:0>{1}{2}}".format(i, size_precursor, 'b')
                    #Generate IDV of product ions
                    idv_masked_by_product = format((int(idv,2) & int("".join(product_mask),2)), 'b')
                    #Numbers of 13C in precursor and product ion
                    precursor_isotope_number = sum([1 for x in idv if x == "1"])
                    product_isotope_number = sum([1 for x in idv_masked_by_product if x == "1"])
                    if (mdv_matrix[precursor_isotope_number][product_isotope_number] == ''):
                        mdv_matrix[precursor_isotope_number][product_isotope_number] = mrm_count
                        mrm_count = mrm_count + 1
                #print mdv_matrix
                #Prepare matrix A and B
                string += "\tA_" + precursor + product + "= numpy.zeros(("+str((mrm_count)*(mrm_count))+"))\n"
                string += "\tB_" + precursor + product + "= numpy.zeros(("+str(mrm_count)+",))\n"
                equation_count = 0

                #precursor ion
                for i in range(size_precursor+1):
                    if (equation_count >= mrm_count):
                        break
                    for j in range(size_product+1):
                        if (mdv_matrix[i][j] != ''):
                            string += "\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][j]) + "] = 1\n"
                    string += "\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_precursor)+"[" + str(emu_intermediate[int(size_precursor)][precursor] * (size_precursor + 1) + i)+"]\n"
                    equation_count = equation_count + 1

                #product ion
                for j in range(size_product+1):
                    if (equation_count >= mrm_count):
                        break
                    for i in range(size_precursor+1):
                        #print i, j, mdv_matrix[i][j]
                        if (mdv_matrix[i][j] != ''):
                            string += "\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][j]) + "] = 1\n"
                    string += "\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_product)+"[" + str(emu_intermediate[int(size_product)][product] * (size_product + 1) + j)+"]\n"
                    equation_count = equation_count + 1

                #neutral loss
                for j in range(size_neutralloss+1):
                    if (equation_count >= mrm_count):
                        break
                    for i in range(size_product+1):
                        #print i, j, mdv_matrix[i][j]
                        if (mdv_matrix[i][i + j] != ''):
                            string += "\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][i + j]) + "] = 1\n"
                    string += "\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_neutralloss)+"[" + str(emu_intermediate[int(size_neutralloss)][neutralloss] * (size_neutralloss + 1) + j)+"]\n"
                    equation_count = equation_count + 1

                string +="\tA_" + precursor + product + ".resize(("+str(mrm_count)+","+ str(mrm_count)+"))\n"
                string += "\tMS_" + precursor + product + " = numpy.linalg.solve(A_" + precursor + product + ", B_" + precursor + product + ")\n"
                string += "\temu_list['" + target_fragment +"'] = [" + ','.join("MS_" + precursor + product + "["+str(x)+"]" for x in range(mrm_count)) + "]\n"
            #gcms
            else:
                emus = target_emus.replace(' ','').split('+')
                # one EMU
                if (len(emus) == 1):
                    emu = emus[0]
                    compound = compound_of_EMU(emu)
                    size = size_of_EMU(emu);
                    row = ["X_"+str(size)+"["+ str(emu_intermediate[int(size)][emu] * (size + 1) + x) + "]"for x in range(size + 1)]
                    string += "\temu_list['" + target_fragment +"'] = [" + ','.join(row) + "]\n"
                    if self.configuration['add_naturalisotope_in_calmdv'] == "yes":
                        if not self.target_fragments[target_fragment]["formula"] == "":
                            string += "\tmdvtemp = [" + ','.join(row) + "]\n"
                            string += "\tmdvcorrected = mdvtemp[:]\n"
                            for i in range(len(row)):
                                textdatatemp = []
                                for j in range(len(row)):
                                    textdatatemp.append(str(self.target_fragments[target_fragment]['natural_isotope_addition'][i,j]) + "* mdvtemp["+str(j)+"]" )
                                string += "\tmdvcorrected[" + str(i) +"] =" + '+'.join(textdatatemp) + "\n"

                            string += "\temu_list['" + target_fragment +"'] = mdvcorrected\n"

                # multiple EMUs
                else:
                    # number of emu
                    numberofemus = len(emus)
                    # list of carbon numbers
                    sizeofemus = [size_of_EMU(x) for x in emus]
                    # list of compounds in emu
                    compoundofemus = [compound_of_EMU(x) for x in emus]
                    # sum of carbon numbers in emus
                    sumofemus = sum(sizeofemus)
                    # All conbination of EMUs
                    emu_list = generate_iterative_list(sizeofemus)
                    row = []
                    for i in range(sumofemus + 1):
                        list_of_emuset = []
                        for emuset in emu_list:
                            if (sum(emuset) == i):
                                #print emuset
                                temp_one_emu = "*".join(["X_"+str(sizeofemus[j])+"["+ str(emu_intermediate[int(sizeofemus[j])][emus[j]] * (sizeofemus[j] + 1) + emuset[j]) + "]" for j in range(len(emuset))])
                                list_of_emuset.append(temp_one_emu)
                        row.append("+".join(list_of_emuset))
                    string += "\temu_list['" + target_fragment +"'] = [" + ','.join(row) + "]\n"
                    if self.configuration['add_naturalisotope_in_calmdv'] == "yes":
                        if not self.target_fragments[target_fragment]["formula"] == "":
                            string += "\tmdvtemp = [" + ','.join(row) + "]\n"
                            string += "\tmdvcorrected = mdvtemp[:]\n"
                            for i in range(len(row)):
                                textdatatemp = []
                                for j in range(len(row)):
                                    textdatatemp.append(str(self.target_fragments[target_fragment]['natural_isotope_addition'][i,j]) + "* mdvtemp["+str(j)+"]" )
                                string += "\tmdvcorrected[" + str(i) +"] =" + '+'.join(textdatatemp) + "\n"

                            string += "\temu_list['" + target_fragment +"'] = mdvcorrected\n"

                string += "\temu_list['X_list'] = X_list\n"
        #generate mdv
        string +=  "\tmdv = []\n"
        string +=  "\tfor emu in target_emu_list:\n"
        string +=  "\t\tif emu in emu_list:\n"
        string +=  "\t\t\tmdv.extend(emu_list[emu])\n"
        string +=  "\treturn(numpy.array(mdv), emu_list)\n"



        ##############################
        # start making diffmdv function
        # diffmdv(r, met, timepoints, target_emu_list, mdv_carbon_sources)
        # r: list of metabolic flux level
        # met: list of metabolite concentration
        # timepoints: list of time points of time course data (sec)
        # target_emu_list: list of emu of target fragments
        # mdv_carbon_source: dict of carbon source MDV (cs.generate_dict())
        #
        #
        string += "\n"
        string += "def diffmdv"
        string += "(r, met, timepoints, target_emu_list, mdv_carbon_sources, y0temp = []):\n"



        #######################
        #
        # start making diffmdv_dxdt
        # Generate python script that difine diffmdv_dxdt(t,y,p, mode = 'no')
        # if mode = 'no'
        #   calc dx = C-1 (AX-BY)
        #
        # if mode = "func"
        #  C-1A and C-1B were calced and then used for the construction of python
        #  script that define diffmdv_dxdt_f(t,y,p)、。
        #


        string += "\n"
        string += "\n"
        string += "\tdef diffmdv_dxdt(y,t,p, mode = 'no'):\n"
        string += "\t\tdx = []\n"

        #
        # for Cython
        #
        cython_string += "\n"
        cython_string += "\n"
        cython_string += "\tdef diffmdv_dxdt_cython(y,t,p):\n"
        cython_string += "\t\tdx = []\n"

        #
        # mode = func
        #
        string += "\t\tstring = 'import numpy as numpy\\n'\n"
        string += "\t\tstring += 'def diffmdv_dxdt_f(y,t,p):\\n'\n"
        string += "\t\tstring += '\tdx = numpy.array([]) \\n'\n"
        #
        # Pointer dictionally that inditate positions of MDVs of each EMU in array y
        #
        emu_intermediates_to_p = {}
        y_position_count = 0;
        #
        # Collection of EMUs in each emu_size
        #
        for emu_size in range(1, int(emu_size_max)+1):
            #
            # start position of this layer's EMU in array y
            #
            start = y_position_count
            # from hit_emu_reactionset, side of emu is equal to emu_size, only one metabolite and not carbon source.
            emu_intermediate_of_this_layer = [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and (not ('+' in x))
                                                                                and not (compound_of_EMU(x) in self.carbon_source)]
            emu_intermediate_of_this_layer.sort()
            #
            # Continue if there is no EMU
            #
            if len(emu_intermediate_of_this_layer) == 0:
                continue
            #
            # positions of MDVs of each EMU in array y
            #
            for emu in emu_intermediate_of_this_layer:
                emu_intermediates_to_p[emu] = y_position_count
                y_position_count = y_position_count + emu_size + 1
            for i, emu in enumerate(emu_intermediate_of_this_layer):
                for num in range(int(emu_size)+1):
                    self.emu_order_in_y.append((emu, num))

            #
            # end position of this layer's EMU in array y
            #
            end = y_position_count
            #
            # Input X from array y
            # X_1 = numpy_array(y[100:200])
            #
            string += "\t\tX_"
            string_temp = str(emu_size) + " = numpy.array(y[" + str(start) +":"+ str(end)+ "])"
            string += string_temp
            string += "\n"
            cython_string += "\t\tX_"
            cython_string_temp = str(emu_size) + " = numpy.array(y[" + str(start) +":"+ str(end)+ "])"
            cython_string += cython_string_temp
            cython_string += "\n"
            #
            string += "\t\tstring += '\\tX_" + string_temp + "\\n'\n"
        #
        # size of y
        #
        size_of_y = y_position_count
        #
        # relationship between flux and array p
        #
        p_position_count = 0;
        for i in range(len(self.reaction_ids)):
            #
            # r1 = p[0]
            #
            string += "\t\t" + self.reaction_ids[i] + " = p[" + str(p_position_count) + "]\n"
            cython_string += "\t\t" + self.reaction_ids[i] + " = p[" + str(p_position_count) + "]\n"
            p_position_count = p_position_count + 1

        #
        # relationship between metabolite and array p
        #
        for metabolite in self.metabolite_ids:
            #
            # Cit = p[0]
            #
            string += "\t\t" + metabolite + " = p[" + str(p_position_count) + "]\n"
            cython_string += "\t\t" + metabolite + " = p[" + str(p_position_count) + "]\n"
            p_position_count = p_position_count + 1
        #
        # relationship between reversible and array p
        #
        for intermediate in self.reversible_ids:
            #
            # Cit = p[0]
            #
            #string += "\t\t" + metabolite + " = p[" + str(p_position_count) + "]\n"
            #cython_string += "\t\t" + metabolite + " = p[" + str(p_position_count) + "]\n"
            p_position_count = p_position_count + 1
        #
        # relationship between carbon source and array p
        #
        string += "\t\tmdv_carbon_sources = {}\n"
        cython_string += "\t\tmdv_carbon_sources = {}\n"
        emu_carbon_source_to_p = {}
        #
        #
        cs_position_count = 0

        for emu_size in range(1, int(emu_size_max)+1):
            #
            # Carbon source in this laryer
            #
            emu_sourse_of_this_layer = [x for x in hit_emuset if (size_of_EMU(x) == emu_size)
                                                                                and (compound_of_EMU(x) in self.carbon_source)]
            emu_sourse_of_this_layer.sort()
            for emu_sourse in emu_sourse_of_this_layer:
                emu_carbon_source_to_p[emu_sourse] = cs_position_count
                #
                # tmdv_carbon_sources[AcCoA_12] = p[14:16]
                #
                string += "\t\tmdv_carbon_sources[\""+emu_sourse+"\"] = p["+str(p_position_count)+":"+str(p_position_count + emu_size + 1)+"]\n"
                cython_string += "\t\tmdv_carbon_sources[\""+emu_sourse+"\"] = p["+str(p_position_count)+":"+str(p_position_count + emu_size + 1)+"]\n"
                p_position_count = p_position_count + emu_size + 1
                cs_position_count = cs_position_count + emu_size + 1
        #
        #
        size_of_p = cs_position_count

        #
        #　AX = BY
        # Essentially identical with that for calmdv
        #
        for emu_size in range(1, int(emu_size_max)+1):
            emu_intermediate_of_this_layer = [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and (not ('+' in x))
                                                                                and not (compound_of_EMU(x) in self.carbon_source)]
            emu_intermediate_of_this_layer.sort()
            if len(emu_intermediate_of_this_layer) == 0:
                continue
            emu_sourse_of_this_layer = []
            emu_sourse_of_this_layer += [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and ('+' in x)]
            emu_sourse_of_this_layer += [x for x in hit_emu_reactionset if (size_of_EMU(x) == emu_size)
                                                                                and (not ('+' in x))
                                                                                and (compound_of_EMU(x) in self.carbon_source)]

            emu_sourse_of_this_layer.sort()

            #　corresponding row of each emu in matrix X
            #
            for i in range(len(emu_intermediate_of_this_layer)):
                if emu_size in emu_intermediate:
                    emu_intermediate[emu_size].setdefault(emu_intermediate_of_this_layer[i],i)
                else:
                    emu_intermediate[emu_size] = {emu_intermediate_of_this_layer[i]:i}
        #Out put A
            size_i = len(emu_intermediate_of_this_layer)
            size_s = len(emu_sourse_of_this_layer)
            #
            # A_1 = numpy.zero((50,))
            #
            string += "\t\tA_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_i)+",))\n"
            cython_string += "\t\tA_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_i)+",))\n"
            for product in emu_intermediate_of_this_layer:
                for substrate in emu_intermediate_of_this_layer:
                    cell = []
                    if substrate in matrix[product]:
                        if substrate == product:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append(reaction_id.replace("*", "") + "/ 2.0")
                                else:
                                    cell.append(reaction_id)
                            cell_string = "-1.0 * (" + ' + '.join(cell) + ")"
                            string += "\t\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            cython_string += "\t\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"

                            if self.configuration['callbacklevel'] == 8:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)

                        else:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append(reaction_id)
                            cell_string = ' + '.join(cell)
                            string += "\t\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            cython_string += "\t\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"

                            if self.configuration['callbacklevel'] == 8:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)
            #
            # Risize A A_1.resize((10,5))
            #
            string += "\t\tA_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"
            cython_string += "\t\tA_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"

            #
            # Resice X X_1.resize((10,5))
            #
            string += "\t\tX_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(emu_size + 1)+"))\n"
            cython_string += "\t\tX_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(emu_size + 1)+"))\n"
            #
            # mode = func
            #
            string += "\t\tstring += '\\tX_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(emu_size + 1)+"))\\n'\n"
        #Output C
            size_i = len(emu_intermediate_of_this_layer)
            size_s = len(emu_sourse_of_this_layer)
            string += "\t\tC_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_i)+",))\n"
            cython_string += "\t\tC_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_i)+",))\n"
            for product in emu_intermediate_of_this_layer:
                for substrate in emu_intermediate_of_this_layer:
                    cell = []
                    if substrate in matrix[product]:
                        if substrate == product:
                            metabolite = compound_of_EMU(product)
                            cell_string = "1.0 / " +  metabolite
                            string += "\t\tC_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            cython_string += "\t\tC_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"

                            if self.configuration['callbacklevel'] == 8:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string )

                        else:
                            pass

            string += "\t\tC_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"
            cython_string += "\t\tC_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"

        #Output B
            string += "\t\tB_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_s)+",))\n"
            cython_string += "\t\tB_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_s)+",))\n"
            for product in emu_intermediate_of_this_layer:
                for substrate in emu_sourse_of_this_layer:
                    cell = []
                    if substrate in matrix[product]:
                        if substrate == product:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("-0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append("-1.0 * " + reaction_id)
                        else:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("-0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append("-1.0 * " + reaction_id)
                    if len(cell) != 0:
                        cell_string = ' + '.join(cell)
                        if self.configuration['callbacklevel'] == 8:
                            print(product,'\t',substrate,"\tB_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_s + emu_sourse_of_this_layer.index(substrate))+"] = " + cell_string)

                        string += "\t\tB_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_s + emu_sourse_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                        cython_string += "\t\tB_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_s + emu_sourse_of_this_layer.index(substrate))+"] = " + cell_string + "\n"

            cython_string += "\t\tB_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_s)+"))\n"

            string += "\t\tB_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_s)+"))\n"
        #Output Y
            string += "\t\tY_" + str(emu_size) + " =  numpy.zeros(("+str(len(emu_sourse_of_this_layer) * (emu_size + 1))+",))\n"
            cython_string += "\t\tY_" + str(emu_size) + " =  numpy.zeros(("+str(len(emu_sourse_of_this_layer) * (emu_size + 1))+",))\n"
            string += "\t\tstring += '\\tY_" + str(emu_size) + " =  numpy.zeros(("+str(len(emu_sourse_of_this_layer) * (emu_size + 1))+",))\\n'\n"
            for substrate in emu_sourse_of_this_layer:#炭素源あるいは一つ前のそうから
                row = []
                #print substrate

                if '+' in substrate:
                    subs_matrix=[]
                    subs_matrix_forfunc=[]
                    list_emusize_plus1=[]
                    n = len(substrate.split('+'))
                    for subs in substrate.split('+'):
                        list_emusize_plus1.append(size_of_EMU(subs)+1)

                        if compound_of_EMU(subs) in self.carbon_source:
                            #sub1_list = generate_carbonsource_MID(sub1)
                            compound = compound_of_EMU(subs)
                            size = size_of_EMU(subs);
                            subs_matrix.append(["mdv_carbon_sources[\"" + subs + "\"][" + str(x) + "]" for x in range(size+1)])
                            subs_matrix_forfunc.append(["'+str(mdv_carbon_sources[\"" + subs + "\"][" + str(x) + "])+'" for x in range(size+1)])

                        else:
                            compound = compound_of_EMU(subs)
                            size = size_of_EMU(subs);
                            subs_matrix.append(["y[" + str(emu_intermediates_to_p[subs]+x) +  "]" for x in range(size+1)])
                            subs_matrix_forfunc.append(["y[" + str(emu_intermediates_to_p[subs]+x) + "]" for x in range(size+1)])



                    sum_length=0
                    for subs_list in subs_matrix:
                        sum_length=sum_length +len(subs_list)
                    #空のリストのリストを作製しておく。
                    equation = [[] for x in range(sum_length)];
                    equation_forfunc = [[] for x in range(sum_length)];
                    conbinaion_list=permute(list_emusize_plus1)
                    for conbination in conbinaion_list:
                       equation[sum(conbination)].extend([("*").join([str(subs_matrix[i][conbination[i]]) for i in range(len(conbination))])])
                       equation_forfunc[sum(conbination)].extend([("*").join([str(subs_matrix_forfunc[i][conbination[i]]) for i in range(len(conbination))])])

                    for numberofisotope in range (len(equation)-n+1):
                        if self.configuration['callbacklevel'] == 8:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]))
                        string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"
                        cython_string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"
                        string += "\t\tstring += '\\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation_forfunc[numberofisotope]) + "\\n'\n"

                else:
                    compound = compound_of_EMU(substrate)
                    size = size_of_EMU(substrate);
                    for numberofisotope in range(size+1):
                        if self.configuration['callbacklevel'] == 8:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "]" )
                        string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "]" + "\n"
                        cython_string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "]" + "\n"
                        string += "\t\tstring += '\\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + "'+str(mdv_carbon_sources[\"" + substrate + "\"][" + str(numberofisotope) + "])+'" + "\\n'\n"


            string += "\t\tY_" + str(emu_size) + ".resize(("+str(size_s)+","+ str(emu_size + 1)+"))\n"
            cython_string += "\t\tY_" + str(emu_size) + ".resize(("+str(size_s)+","+ str(emu_size + 1)+"))\n"

            string += "\t\tstring += '\\tY_" + str(emu_size) + ".resize(("+str(size_s)+","+ str(emu_size + 1)+"))\\n'\n"
            #
            # calc dx
            # Representative python script generated when mode == 'func'
            # if mode == 'func':
            #   CA_1 = numpy.dot(C_1, A_1)
            #   CB_1 = numpy.dot(C_1, B_1)
            #   temp = numpy.rabel(CA_1)
            #   string += '\\tCA_1 = numpy.zeros(50)\\n'
            #   for i in range(len(temp)):
            #       if temp[i] == 0.0:
            #           continue
            #       string += '\\tCA_1['+str(i)+'] = '+str(temp[i])+'\\n'
            #   string += '\\tCA_1.resize((50,20))\\n'
            string += "\t\tif mode == 'func':\n"
            string += "\t\t\tCA_"+str(emu_size)+"= numpy.dot(C_"+str(emu_size)+",A_"+str(emu_size)+")\n"
            string += "\t\t\tCB_"+str(emu_size)+"= numpy.dot(C_"+str(emu_size)+",B_"+str(emu_size)+")\n"
            string += "\t\t\ttemp = numpy.ravel(CA_"+str(emu_size)+")\n"
            string += "\t\t\tstring += '\\tCA_"+str(emu_size)+" = numpy.zeros('+str(len(temp))+')\\n'\n"
            string += "\t\t\tfor i in range(len(temp)):\n"
            string += "\t\t\t\tif temp[i] == 0.0:\n"
            string += "\t\t\t\t\tcontinue\n"
            string += "\t\t\t\tstring += '\\tCA_"+str(emu_size)+"['+str(i)+'] = '+str(temp[i])+'\\n'\n"
            string += "\t\t\tstring += '\\tCA_"+str(emu_size)+".resize(("+str(size_i)+","+str(size_i)+"))\\n'\n"
            #
            # Repeat for CB
            #
            string += "\t\t\ttemp = numpy.ravel(CB_"+str(emu_size)+")\n"
            string += "\t\t\tstring += '\\tCB_"+str(emu_size)+" = numpy.zeros('+str(len(temp))+')\\n'\n"
            string += "\t\t\tfor i in range(len(temp)):\n"
            string += "\t\t\t\tif temp[i] == 0.0:\n"
            string += "\t\t\t\t\tcontinue\n"
            string += "\t\t\t\tstring += '\\tCB_"+str(emu_size)+"['+str(i)+'] = '+str(temp[i])+'\\n'\n"
            string += "\t\t\tstring += '\\tCB_"+str(emu_size)+".resize(("+str(size_i)+","+str(size_s)+"))\\n'\n"
            #
            # string += '\\tdx_1 = numpy.dot(CA_1,X_1) - numpy.dot(CB_1,Y_1)\\n'
            # string += '\\tdx = numpy.append(dx, numpy.ravel(dx_1))\\n'
            #
            string += "\t\t\tstring += '\\tdx_"+str(emu_size)+" = numpy.dot(CA_"+str(emu_size)+",X_"+str(emu_size)+") - numpy.dot(CB_"+str(emu_size)+",Y_"+str(emu_size)+")\\n'\n"
            string += "\t\t\tstring += '\\tdx = numpy.append(dx, numpy.ravel(dx_"+str(emu_size)+"))\\n'\n"
            #
            # dx_1 = numpy.dot(C_1, (numpy.dot(A_1,X_1) - numpy.dot(B_1,Y_1)))
            # dx_1.tolist()
            # dx.extend([x for sublist in dx_1 for x in sublist])
            #
            string += "\t\tdx_"+str(emu_size)+" = numpy.dot(C_"+str(emu_size) +", (numpy.dot(A_"+str(emu_size)+",X_"+str(emu_size)+") - numpy.dot(B_"+str(emu_size)+",Y_"+str(emu_size)+")))\n"
            string += "\t\tdx_"+str(emu_size)+".tolist()\n"
            string += "\t\tdx.extend([x for sublist in dx_"+str(emu_size)+" for x in sublist])\n"
            cython_string += "\t\tdx_"+str(emu_size)+" = numpy.dot(C_"+str(emu_size) +", (numpy.dot(A_"+str(emu_size)+",X_"+str(emu_size)+") - numpy.dot(B_"+str(emu_size)+",Y_"+str(emu_size)+")))\n"
            cython_string += "\t\tdx_"+str(emu_size)+".tolist()\n"
            cython_string += "\t\tdx.extend([x for sublist in dx_"+str(emu_size)+" for x in sublist])\n"

        #
        # return
        #
        string += "\t\tstring += '\\treturn(dx)\\n'\n"
        string += "\t\tif mode == 'func':\n"
        string += "\t\t\treturn(string)\n"
        string += "\t\treturn(dx)\n"
        cython_string += "\t\treturn(dx)\n"


        #
        #return making diffmdc
        #
        string += "\timport scipy as scipy\n"
        string += "\temu_list = {}\n"

        #
        # Definition of p0 (list of parameters)
        #
        # append list of metabolic flux level
        # p0.extend(r)
        string += "\tp0 = list(r)\n"
        #
        # append list of metabolite concentration
        # p0.extend(met)
        string += "\tp0.extend(met)\n"
        #
        # Set MDV
        # csp = [0.0] * 50
        string += "\tcsp = [0.0] * "+str(size_of_p)+"\n"
        #
        # Generate csp from emu_carbon_source_to_p
        #
        for (carbon_source, position) in emu_carbon_source_to_p.items():
            # Check EMU size
            size = size_of_EMU(carbon_source)
            # csp[10:15] = list(mdv_carbon_sources["AcCoA_12"])
            string += "\tcsp["+str(position)+":"+ str(position + size + 1)+"] = list(mdv_carbon_sources[\""+carbon_source+"\"])\n"
        # p0.extend(csp)
        string += "\tp0.extend(csp)\n"
        #
        # Initial EMU
        # y0 = [0.0] * 300
        string += "\ty0 = [0.0] * "+str(size_of_y)+"\n"
        #
        # Generate initial MDV from emu_intermediates_to_p
        # Considering natural isotope。
        #

        for (emu, position) in emu_intermediates_to_p.items():

            #
            # y0[13] = 1.0
            #
            (stem, pos) = emu.split("_")
            number_of_carbons = len(pos)
            #string += "\ty0["+str(position)+"] = 1.0\n"
            string += "\ty0["+str(position)+"] = "
            string += str(0.98893**number_of_carbons)
            string += "#"+str(emu)+" "+str(position)
            string += "\n"
            string += "\ty0["+str(position+1)+"] = "
            string += str((0.98893**(number_of_carbons-1.0))*0.01106*number_of_carbons)
            string += "\n"
            if number_of_carbons > 2:
                string += "\ty0["+str(position+2)+"] = "
                string += str(1.0- 0.98893**number_of_carbons - (0.98893**(number_of_carbons-1.0))*0.01106*number_of_carbons)
                string += "\n"
            #print emu,position,stem, pos
        # emu_list = {}\
        string += "\tif len(y0temp) == "+str(size_of_y)+":\n"
        string += "\t\ty0 = y0temp[:]\n"

        """
        for (emu, position) in emu_intermediates_to_p.items():
            (stem, pos) = emu.split("_")
            number_of_carbons = len(pos)
            for posx in range(number_of_carbons + 1):
                string += "\t\ty0["+str(position + posx)+"] = "
                string += "y0temp["+str(position + posx)+"]"
                string += "#"+str(emu)+" "+str(position)
                string += "\n"
        """
        string += "\temu_list = {}\n"
        #
        # Generate  diffmdv_dxdt(y0, 0, p0, mode='func')
        #scipy
        string += "\tdiffmdv_func = diffmdv_dxdt(y0, 0, p0, mode='func')\n"
        #assimulo
        #string += "\tdiffmdv_func = diffmdv_dxdt(0, y0, p0, mode='func')\n"
        #
        # Save diff_mdv
        #
        if self.configuration['callbacklevel'] >= 3:
            string += "\tf = open('diffmdv_func.py', 'w')\n"
            string += "\tf.write(diffmdv_func)\n"
            string += "\tf.close()\n"
        #
        #temp
        #
        #
        string += "\tlocals_dic = locals()\n"
        string += "\texec(diffmdv_func, globals(), locals_dic)\n"
        string += "\tdiffmdv_dxdt_f = locals_dic['diffmdv_dxdt_f']\n"
        string += "\timport scipy as scipy\n"
        #
        # Perform integration
        #
        # Odeint scipy
        #
        #string += "\ty = odeint(diffmdv_dxdt_f, y0, timepoints, args=(p0,), rtol=1e-3, atol=1e-3)\n"
        if mode == "cython":
            string += "\ty = scipy.integrate.odeint(diffmdv_dxdt_cython, y0, timepoints, args=(p0,), rtol=1e-2, atol=1e-2)\n"
        else:
            string += "\ty = scipy.integrate.odeint(diffmdv_dxdt_f, y0, timepoints, args=(p0,), rtol=1e-3, atol=1e-3)\n"
        #

        #
        # assimulo
        #
        #string += "\texp_mod = assimulo.problem.Explicit_Problem(diffmdv_dxdt_f, y0)\n"
        #string += "\texp_mod = assimulo.problem.Explicit_Problem()\n"
        #string += "\texp_mod.name = 'Simple CVode Example'\n"
        #string += "\texp_mod.rhs = diffmdv_dxdt_f\n"
        #string += "\texp_mod.y0 = y0\n"
        #string += "\texp_mod.p0 = p0\n"
        #string += "\texp_sim = assimulo.solvers.CVode(exp_mod)\n"

        #string += "\texp_sim.iter  = 'Newton'\n"
        #string += "\texp_sim.discr = 'BDF'\n"
        #string += "\texp_sim.verbosity = 50\n"
        #string += "\texp_sim.maxord = 3\n"
        #string += "\texp_sim.usesens = False\n"
        #string += "\texp_sim.sensmethod = 'SIMULTANEOUS'\n"
        #string += "\texp_sim.suppress_sens = True\n"
        #string += "\texp_sim.atol = 1e-3\n"
        #string += "\texp_sim.rtol = 1e-3\n"
        #
        # Pefrom integration
        #
        #string += "\tt1, y = exp_sim.simulate(max(timepoints),0,timepoints)\n"
        #　
        #　Initialize list to store time course data
        #  emu_list['Glu_12345'] =[[time point 1,,,,], [time point 2,,,,], [time point 3,,,,], [time point 3,,,,],]
        #
        for emu in target_fragments:
            string += "\temu_list[\'"+emu+"\'] = []\n"
        #for i, time in enumerate(t1):
        #assimulo
        #string += "\tfor i, time in enumerate(t1):\n"
        string += "\tfor i, time in enumerate(timepoints):\n"
        #
        # Calc MDVs of target fradments
        #
        for target_fragment in target_fragments:
            #
            # Ignore
            #
            if target_fragments[target_fragment]['use'] != 'use':
                continue
            target_emus = target_fragments[target_fragment]['atommap']
            #
            # "msms" type
            #
            if target_fragments[target_fragment]['type'] == "msms":
                precursor, neutralloss, product = target_emus.replace(' ','').split('+')
                # Get sizes
                size_precursor = size_of_EMU(precursor);
                size_product = size_of_EMU(product);
                size_neutralloss = size_of_EMU(neutralloss);
                product_carbon_number = list(product.split('_')[1]);
                product_mask = [];
                for i in (range(size_precursor)):
                    product_mask.append("0")
                for i in (product_carbon_number):
                    product_mask[size_precursor - int(i)] = "1"
                numberofidv = 2 ** size_precursor
                mdv_matrix = [];
                for pre in (range(0, size_precursor + 1)):
                    temp = []
                    for pro in (range(0, size_product + 1)):
                        temp.append('')
                    mdv_matrix.append(temp)
                mrm_count = 0
                for i in (sorted(range(numberofidv))):
                    idv="{0:0>{1}{2}}".format(i, size_precursor, 'b')
                    idv_masked_by_product = format((int(idv,2) & int("".join(product_mask),2)), 'b')
                    precursor_isotope_number = sum([1 for x in idv if x == "1"])
                    product_isotope_number = sum([1 for x in idv_masked_by_product if x == "1"])
                    if (mdv_matrix[precursor_isotope_number][product_isotope_number] == ''):
                        mdv_matrix[precursor_isotope_number][product_isotope_number] = mrm_count
                        mrm_count = mrm_count + 1
                #print mdv_matrix
                string += "\t\tA_" + precursor + product + "= numpy.zeros(("+str((mrm_count)*(mrm_count))+"))\n"
                string += "\t\tB_" + precursor + product + "= numpy.zeros(("+str(mrm_count)+",))\n"
                equation_count = 0

                #precursor ion
                for i in range(size_precursor+1):
                    if (equation_count >= mrm_count):
                        break
                    for j in range(size_product+1):
                        if (mdv_matrix[i][j] != ''):
                            string += "\t\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][j]) + "] = 1\n"
                    #string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_precursor)+"[" + str(emu_intermediate[int(size_precursor)][precursor] * (size_precursor + 1) + i)+"]\n"
                    string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = y[i]["+ str(emu_intermediates_to_p[precursor] + i) + "]\n"

                    equation_count = equation_count + 1

                #product ion
                for j in range(size_product+1):
                    if (equation_count >= mrm_count):
                        break
                    for i in range(size_precursor+1):
                        #print i, j, mdv_matrix[i][j]
                        if (mdv_matrix[i][j] != ''):
                            string += "\t\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][j]) + "] = 1\n"
                    #string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_product)+"[" + str(emu_intermediate[int(size_product)][product] * (size_product + 1) + j)+"]\n"
                    string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = y[i]["+ str(emu_intermediates_to_p[product] + j) + "]\n"

                    equation_count = equation_count + 1

                #neutral loss
                for j in range(size_neutralloss+1):
                    if (equation_count >= mrm_count):
                        break
                    for i in range(size_product+1):
                        #print i, j, mdv_matrix[i][j]
                        if (mdv_matrix[i][i + j] != ''):
                            string += "\t\tA_" + precursor + product+"[" + str(equation_count * (mrm_count) + mdv_matrix[i][i + j]) + "] = 1\n"
                    #string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = X_"+str(size_neutralloss)+"[" + str(emu_intermediate[int(size_neutralloss)][neutralloss] * (size_neutralloss + 1) + j)+"]\n"
                    string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = y[i]["+ str(emu_intermediates_to_p[neutralloss] + j) + "]\n"

                    equation_count = equation_count + 1

                string +="\t\tA_" + precursor + product + ".resize(("+str(mrm_count)+","+ str(mrm_count)+"))\n"
                string += "\t\tMS_" + precursor + product + " = numpy.linalg.solve(A_" + precursor + product + ", B_" + precursor + product + ")\n"
                string += "\t\temu_list['" + target_fragment +"'].append([" + ','.join("MS_" + precursor + product + "["+str(x)+"]" for x in range(mrm_count)) + "])\n"
            #
            # "gcms"
            #
            else:
                emus = target_emus.replace(' ','').split('+')
                # only one EMU
                if (len(emus) == 1):
                    emu = emus[0]
                    compound = compound_of_EMU(emu)
                    size = size_of_EMU(emu);
                    #row = ["X_"+str(size)+"["+ str(emu_intermediate[int(size)][emu] * (size + 1) + x) + "]"for x in range(size + 1)]
                    row = ["y[i]["+ str(emu_intermediates_to_p[emu] + x) + "]"for x in range(size + 1)]

                    string += "\t\temu_list['" + target_fragment +"'].append([" + ','.join(row) + "])\n"
                    if self.configuration['add_naturalisotope_in_calmdv'] == "yes":
                        if not self.target_fragments[target_fragment]["formula"] == "":
                            string += "\t\tmdvtemp = [" + ','.join(row) + "]\n"
                            string += "\t\tmdvcorrected = mdvtemp[:]\n"
                            for i in range(len(row)):
                                textdatatemp = []
                                for j in range(len(row)):
                                    textdatatemp.append(str(self.target_fragments[target_fragment]['natural_isotope_addition'][i,j]) + "* mdvtemp["+str(j)+"]" )
                                string += "\t\tmdvcorrected[" + str(i) +"] =" + '+'.join(textdatatemp) + "\n"

                            string += "\t\temu_list['" + target_fragment +"'][-1] = mdvcorrected[:]\n"



                # multiple EMU
                else:
                    numberofemus = len(emus)
                    sizeofemus = [size_of_EMU(x) for x in emus]
                    compoundofemus = [compound_of_EMU(x) for x in emus]
                    sumofemus = sum(sizeofemus)
                    emu_list = generate_iterative_list(sizeofemus)
                    row = []
                    for i in range(sumofemus + 1):
                        list_of_emuset = []
                        for emuset in emu_list:
                            if (sum(emuset) == i):
                                #print emuset
                                temp_one_emu = "*".join(["X_"+str(sizeofemus[j])+"["+ str(emu_intermediate[int(sizeofemus[j])][emus[j]] * (sizeofemus[j] + 1) + emuset[j]) + "]" for j in range(len(emuset))])

                                temp_one_emu = "*".join(["y[i]["+ str(emu_intermediates_to_p[emus[j]] + emuset[j]) + "]" for j in range(len(emuset))])
                                #string += "\t\tB_" + precursor + product+"[" + str(equation_count )+ "] = y[i]["+ str(emu_intermediates_to_p[precursor] + i) + "]\n"

                                list_of_emuset.append(temp_one_emu)
                        row.append("+".join(list_of_emuset))
                    string += "\t\temu_list['" + target_fragment +"'].append([" + ','.join(row) + "])\n"
                    if self.configuration['add_naturalisotope_in_calmdv'] == "yes":
                        if not self.target_fragments[target_fragment]["formula"] == "":
                            string += "\t\tmdvtemp = [" + ','.join(row) + "]\n"
                            string += "\t\tmdvcorrected = mdvtemp[:]\n"
                            for i in range(len(row)):
                                textdatatemp = []
                                for j in range(len(row)):
                                    textdatatemp.append(str(self.target_fragments[target_fragment]['natural_isotope_addition'][i,j]) + "* mdvtemp["+str(j)+"]" )
                                string += "\t\tmdvcorrected[" + str(i) +"] =" + '+'.join(textdatatemp) + "\n"

                            string += "\t\temu_list['" + target_fragment +"'][-1] = mdvcorrected[:]\n"


        string += "\tmdv = []\n"
        #string += "\tfor i, time in enumerate(t1):\n"
        string += "\tfor i, time in enumerate(timepoints):\n"
        string += "\t\tfor emu in target_emu_list:\n"
        string += "\t\t\tif emu in emu_list:\n"
        string += "\t\t\t\tif emu_list[emu]==[]:\n"
        string += "\t\t\t\t\t continue\n"
        string += "\t\t\t\telse:\n"
        string += "\t\t\t\t\tmdv.extend(emu_list[emu][i])\n"



        string +=  "\treturn(numpy.array(mdv), emu_list)\n"

        return (string, cython_string)

    def set_configuration(self, *args, **kwargs):
        """
        Set configurations of model parameters.

        Parameters
        ----------
        callbacklevel:
            'callbacklevel' determines details of callbacks from metabolic model. default: 0
            It is useful for error shooting.
            0: none, 1: important things, 2: detailed
        iteration_max = 100,
            Maximal number of interations in model fitting. default: 100


        Examples
        --------
        >>> model.set_configuration(callbacklevel = 1)
        >>> model.set_configuration(iteration_max = 1000)


        See Also
        --------

        """

        for word in kwargs:
            #if word in self.configuration:
            self.configuration[word] = kwargs[word]
            if self.configuration['callbacklevel'] >= 6:
                print(word,'is set at',kwargs[word])

    def set_constrain(self, group, id, type, value = 0.1, stdev = 1.0):
        """
        Setter of a metabolic constrains.
        Please perform update() after model modification.

        Parameters
        ----------
        group :group of parameters in "reaction", "reversible", "metabolite"
        id : reaction id to be fixed (str)
        type: type of constrain "fitting", "fixed", "free", "pseudo"
            "fitting" : For the calclation of RSS by using "value" and stdev
            "fixed"   : Fixed to "value". Not used for the calculation of RSS.
            "free"    : Free between lower and upper boundary. Not used for the calculation of RSS.
            "pseudo"  : Pseudo reaction is a special "free" reaction to consider G-index.
            The reaction is ignored in a flux stoichiometry. However, the reaction is considered
            in the MDV calculation.


        value : flixed level of metabolic flux (float)
        stdev : standard deviation of metabolic flux level (optional the information will be used when the reaction is 'fitting' type)

        Examples
        --------
        >>> model.set_fixed_reaction("reaction", 'pgi',"fixed", 100.0,1)
        Glucuse intake reaction rate is fixed to 100.0

        See Also
        --------
        get_constrain

        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 0:
                print('ERROR:', group,' is wrong group')
                return False

        if type not in ["fitting", "fixed", "free", "pseudo"]:
            if self.configuration['callbacklevel'] >= 0:
                print('ERROR:', type, ' is wrong type')
                return False

        if id in dic_temp:
            dic_temp[id]['type'] = type
            dic_temp[id]['value'] = float(value)
            if float(stdev) >= 0:
                dic_temp[id]['stdev'] = float(stdev)
            if self.configuration['callbacklevel'] >= 6:
                print(id,' is set as ',type, group, ' with value at ',value, ' and stdev at ', stdev, '.')
            return True
        if self.configuration['callbacklevel'] >= 0:
            print('ERROR:', id,' not existed in the model')
        return False

    def get_constrain(self, group, id):
        """
        Getter of constrains of metabolic reaction.

        Examples
        --------
        >>> type, value, stdev = model.set_fixed_reaction('HEX1', 100.0,1)


        Parameters
        ----------
        group :group of parameters in "reaction", "reversible", "metabolite"
        id : reaction id to be fixed (str)


        Reterns
        ----------
        type: type of constrain "fitting", "fixed", "free", "pseudo"
        value : flixed level of metabolic flux (float)
        stdev : standard deviation of metabolic flux level (optional the information will be used when the reaction is 'fitting' type)



        See Also
        --------
        set_constrain

        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 0:
                print('ERROR:', group,' is wrong group')
                return False, False, False

        if id in dic_temp:
            types = dic_temp[id]['type']
            value = float(dic_temp[id]['value'])
            stdev = float(dic_temp[id]['stdev'])
            return types, value, stdev
        if self.configuration['callbacklevel'] >= 0:
            print('ERROR:', id,' not existed in the model')
        return False, False, False

    def set_boundary(self, group, id, lb, ub):
        """
        Set a lower and upper boundaries of metablic states
        Please perform update() after model modification.

        Parameters
        ----------
        group : metabolite, reaction, reversible (str)
        id :  id of state parameter (str)
        lb : lower boundary (float)
        ub : upper boundary (float)


        Examples
        --------
        >>> model.set_boundary("reaction", 'HEX1', 0.001, 100.0)

        See Also
        --------
        set_constrain
        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 0:
                print('ERROR:', group,' is wrong group')
                return False


        if id in dic_temp:
            dic_temp[id]['lb'] = float(lb)
            dic_temp[id]['ub'] = float(ub)
            if self.configuration['callbacklevel'] >= 6:
                print("lower and upper boundaries of", id,' are set to be ', lb, "and", ub)
            return True
        if self.configuration['callbacklevel'] >= 0:
            print('ERROR:', id,' not existed in the model')
        return False


    def set_active_reaction(self, id, type = "use"):
        """
        Set a metabolic reaction as "unuse" type. The reaction is completely ignored during the model construction.
        This is an experimental function.
        Please perform reconstruct() after model modification.

        Parameters
        ----------
        id : reaction id to be activated or inactivated (str)

        Examples
        --------
        >>> model.set_active_reaction('r2')

        See Also
        --------


        """
        if id in self.reactions:
            if type == "use":
                self.reactions[id]['use'] = 'use'
                if self.configuration['callbacklevel'] >= 6:
                    print(id,' was activated in the model.')
                return True
            else:
                self.reactions[id]['use'] = 'no'
                if self.configuration['callbacklevel'] >= 6:
                    print(id,' was inactivated in the model.')
                return True
        if self.configuration['callbacklevel'] >= 0:
            print('ERROR:', id,' not existed in the model')
        return False


    def set_constraints_from_state_dict(self, dict):
        """
        Reaction types, value, stdev, lb, and ub are set from the state dict data

        Parameters
        ----------
        flux_dict: Dictionary of flux. Information in 'value', 'stdev', 'lb', 'ub' and 'type' field was used.

        Reterns
        ----------


        Examples
        --------
        >>> model.set_constraints_from_state_dict(flux_dict)

        See Also
        --------

        """
        #
        # preparation of data
        #

        for i, (group, id) in enumerate(self.vector["ids"]):
            if id in dict[group]:
                type = dict[group][id]['type']
                value = dict[group][id]['value']
                stdev = dict[group][id]['stdev']
                lb = dict[group][id]['lb']
                ub = dict[group][id]['ub']
                self.set_constrain(group, id, type, value = value, stdev = stdev)
                self.set_boundary(group, id, lb, ub)

        self.update()
        return True


    def generate_state_dict(self, tmp_r):
        """
        Generage state dict from given vector.
        Parameters
        ----------
        tmp_r: tmp_r = numpy.dot(matrixinv, Rm_intial)

        Returns
        ----------
        state_dict : Dictionary of metabolic state data.

        Examples
        --------
        >>> state_dict = model.generate_state_dict(tmp_r)

        See Also
        --------
        """
        #
        # preparation of data
        #
        flux_dict = {}
        for i, id in enumerate(self.reaction_ids):
            flux_dict[id] = {
                'value': tmp_r[self.reactions[id]['position_in_tmp_r']],
                'stdev': self.reactions[id]['stdev'],
                'type': self.reactions[id]['type'],
                'reversible':self.reactions[id]['reversible'],
                'order':self.reactions[id]['order'],
                'lb':self.reactions[id]['lb'],
                'ub':self.reactions[id]['ub'],
            }
        conc_dict = {}
        for i, id in enumerate(self.metabolites):
            conc_dict[id] = {
                'value': self.metabolites[id]['value'],
                'stdev': self.metabolites[id]['stdev'],
                'type':self.metabolites[id]['type'],
                'order':self.metabolites[id]['order'],
                'lb':self.metabolites[id]['lb'],
                'ub':self.metabolites[id]['ub'],
            }
            if id in self.metabolite_ids:
                conc_dict[id]['value'] = tmp_r[self.metabolites[id]['position_in_tmp_r']]
        reversible_dict = {}
        for i, id in enumerate(self.reversible_ids):
            reversible_dict[id] = {
                'value':tmp_r[self.reversible[id]['position_in_tmp_r']],
                'stdev':self.reversible[id]['stdev'],
                'type':self.reversible[id]['type'],
                'order':self.reversible[id]['order'],
                'lb':self.reversible[id]['lb'],
                'ub':self.reversible[id]['ub'],
            }

        return {"reaction":flux_dict, "metabolite":conc_dict, "reversible": reversible_dict}


    def generate_carbon_source_templete(self):
        """
        Generate a templete instance for storing a labelling information of carbon source.
        Carbon source information is derived from the metabolic model (//Metabolites)

        Parameters
        ----------
        nothing


        Examples
        --------
        >>> carbon_source_idv = model.generate_carbon_source_templete()

        See Also
        --------
        CarbonSource.py


        """
        #
        #
        cs = {}
        #for each line in mass data
        for compound in self.carbon_source:
            cs[compound] = {
            'IDV': self.carbon_source[compound]['IDV'][:],
            'size': self.carbon_source[compound]['size']
            }
        # Initial data is full 12C without natural 13C
        for compound in cs:
            cs[compound]['IDV'][0] = 1.0
        return carbonsource.CarbonSource(cs)

    def generate_mdv(self, flux, carbon_sources, timepoint = [], startidv = []):
        """
        Generate a MdvData instance including a mass distribution vector (MDV) data
        under 'flux' and 'carbon_sources' condiction.

        Parameters
        ----------
        flux: Dictionary of metabolic flux distribution
        carbon_sources: Instance of carbon source object
        timepoint: Array of time points of measurement in IDV
        startidv:

        Returns
        ----------
        mdv : MdvData instance

        Examples
        --------
        >>> mdv = model.generate_mdv(flux, mdv_carbon_sources)

        See Also
        --------
        generate_carbonsource_MDV

        History
        --------
        140912 calc_MDV_from_flux instead of calmdv is used.

        """

        tmp_r = [flux[group][id]['value'] for (group, id) in self.vector["ids"]]
        mdv_carbon_sources = carbon_sources.generate_dict()
        #Call calmdv via calc_MDV_from_flux function in mfapy.optimize
        if len(timepoint) == 0:
            mdv_exp, mdv_hash = optimize.calc_MDV_from_flux(tmp_r, self.target_fragments.keys(), mdv_carbon_sources, self.func)
            mdv_data = mdv.MdvData(self.target_fragments)
            for fragment, item in mdv_hash.items():
                for number in range(len(item)):
                    if mdv_data.has_data(fragment, number):
                        mdv_data.set_data(fragment, number, mdv_hash[fragment][number], 1.0, "use")
            return mdv_data
        else:
            startidv_temp = []
            if len(startidv) == len(self.emu_order_in_y):
                startidv_temp = list(startidv)
            mdv_exp, mdv_hash = optimize.calc_MDV_from_flux(tmp_r, self.target_fragments.keys(), mdv_carbon_sources, self.func, timepoint = timepoint, y0temp = startidv_temp)
            #
            #mdv_timecourseインスタンスの生成
            #
            mdv_timecourse = mdv.MdvTimeCourseData()
            for point in timepoint:
                mdv_timecourse.add_time_point(point,mdv.MdvData(self.target_fragments))
            for fragment, item in mdv_hash.items():
                for i,ratio in enumerate(item):
                    for number in range(len(ratio)):
                        if mdv_timecourse.has_data(timepoint[i], fragment, number):
                            mdv_timecourse.set_data(timepoint[i], fragment, number, mdv_hash[fragment][i][number], 1.0, "use")
            return mdv_timecourse


    def set_experiment(self, name, mdv, carbon_sources, startidv = []):
        """
        Set an 'experiment' to metabolic model. Here an 'experiment' indicated
        a set of a carbon source labeling pattern and a measured MDV data.
        Parallel labeling experiment can be performed by setting multiple sets of 'experiment'.


        Parameters
        ----------
        name : name of the experiment (unique)
        mdv : MDVdata instance including measured MDV data of target fragments.
        carbon_sources: Instance of carbon source object
        startidv: array of idv as starting isotope distribution for INST


        Examples
        --------
        >>> model.set_experiment('ex1', mdv1, cs1)
        >>> model.set_experiment('ex2', mdv2, cs2, startidv)

        See Also
        --------


        """

        ids, mdv_exp_original, mdv_std_original, mdv_use, target_emu_list, rawdata = mdv.generate_observed_mdv()
        number_of_measurement = mdv.get_number_of_measurement()
        mdv_carbon_sources = carbon_sources.generate_dict()
        self.experiments[name]={
            'mode':"ST",
            'mdv_exp_original':list(mdv_exp_original),
            'mdv_std_original':list(mdv_std_original),
            'mdv_use':list(mdv_use),
            'mdv_ids':list(ids),
            'target_emu_list':list(target_emu_list),
            'mdv_carbon_sources':mdv_carbon_sources,
            'number_of_measurement':number_of_measurement
        }
        if mdv.mode == "timecourse":
            self.experiments[name]["mode"] = "INST"
            self.experiments[name]["timepoint"] = mdv.get_timepoints()
            #
            # Set initial y0 in diffmdv
            #
            self.experiments[name]["y0"] = [0.0] * len(self.emu_order_in_y)
            if len(startidv) == len(self.emu_order_in_y):
                self.experiments[name]["y0"] = list(startidv)
            else:
                for position, emu in enumerate(self.emu_order_in_y):
                    #
                    if emu[1] == 0:
                        (stem, pos) = emu[0].split("_")
                        number_of_carbons = len(pos)
                        self.experiments[name]["y0"][position] = 0.98893**number_of_carbons
                        self.experiments[name]["y0"][position + 1] = (0.98893**(number_of_carbons-1.0))*0.01106*number_of_carbons
                        if number_of_carbons > 2:
                            self.experiments[name]["y0"][position + 2] = 1.0- 0.98893**number_of_carbons - (0.98893**(number_of_carbons-1.0))*0.01106*number_of_carbons
        if self.configuration['callbacklevel'] >= 5:
            print("Set experiment: ", name,' was added to the metabolic model.')
        return True

    def calc_idv(self, flux, carbon_sources):
        """
        Calc IDV from given flux and carbon_sources.
        Examples
        --------
        >>> idv = model.calc_idv(flux, carbon_sources)

        Parameters
        ----------
        flux : dic of metabolic state of initial state
        carbon_sources: Instance of carbon source object

        Parameters
        ----------
        idv : array of idv


        See Also
        --------
        set_experiment()
        generate_mdv()



        """

        calmdv = self.func["calmdv"]
        matrixinv=self.matrixinv
        Rm_initial= self.vector["Rm_initial"]
        stoichiometric_num = self.numbers['independent_start']
        reaction_num= self.numbers['total_number']
        Rm_ind = [flux[group][id]["value"] for (group, id) in self.vector['independent_flux']]
        Rm = numpy.array(list(Rm_initial))
        Rm[stoichiometric_num: reaction_num] = list(Rm_ind)
        tmp_r = numpy.dot(matrixinv, Rm)
        target_emu_list = list(self.target_fragments)
        mdv_carbon_sources = carbon_sources.mdv_carbon_sources

        mdv_original_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)

        X = mdv_hash['X_list']
        for i, x in enumerate(X):
            if x <= 0:
                X[i] = 0.0
            if x >= 1.0:
                X[i] = 1.0

        X_dict = dict(zip(self.emu_order_in_X, X))
        y0 = [0.0] * len(self.emu_order_in_y)

        for position, emu in enumerate(self.emu_order_in_y):
            y0[position] = X_dict[emu]
        return y0


    def clear_experiment(self):
        """
        Clear all 'experiment set(s)' in the metabolic model.


        Parameters
        ----------


        Examples
        --------
        >>> model.clear_experiment()


        See Also
        --------
        model.set_experiment()

        History
        --------
        Newly developed at 4/9/2014


        """
        names = list(self.experiments.keys())
        self.experiments = {}
        if self.configuration['callbacklevel'] >= 5:
            print("Clear experiment: ", names,' are removed from the metabolic model.')

        return True

    def generate_flux_distribution(self, template = []):
        """
        Metabolic flux distribution is randomly generated.
        The metabolic flux distribution can be a initial flux data for model fitting.

        Parameters
        ----------
        template: Dictionary of metabolic state. When template is available, metabolic state most similar to the template is generated. The function is used in the grid search.

        Returns
        ----------
        flux : Dictionary of metabolic flux data.
        independent_flux : List of independent flux

        Examples
        --------
        >>> flux, independent_flux = model.generate_intial_flux_distribution()

        See Also
        --------

        """
        #
        # Set parameters
        #
        numbers = copy.deepcopy(self.numbers)
        vectors = copy.deepcopy(self.vector)
        configuration = copy.deepcopy(self.configuration)

        matrixinv = self.matrixinv
        initial_search_iteration_max = configuration["initial_search_iteration_max"]

        ub = [self.reactions[x]['ub'] for x in self.reaction_ids]
        lb = [self.reactions[x]['lb'] for x in self.reaction_ids]
        if len(template) > 0:
            template = [template[type][id]["value"] for (type, id) in self.vector["ids"]]

        tmp_r, Rm_temp, Rm_ind, state = optimize.initializing_Rm_fitting(numbers, vectors, matrixinv, template ,initial_search_iteration_max)

        return self.generate_state_dict(tmp_r), state


    def generate_initial_states(self, iterations = 100, initial_states = 1, method = 'normal', template = []):
        """
        Initial metabolic states are randomly generated for "iterations" times from which
        better states with lower RSS (number_of_initial_states) were selected.
        The metabolic flux distribution is a initial flux data for model fitting.

        Parameters
        ----------
        method:
            'normal' : Sequencial generation of flux distributions.
            'parallel' : Parallel generation of flux distributions using parallel python.
        iterations:
            Number of trials to generate initial metabolic state. A trial sometimes failed due to a bad optimization state.
        initial_states:
            Number of initial metabolic state(s) to be generated.
        template:
            Templete metabolic state(s) dict. When template is available, initial metabolic state(s) similiar to
            the templete metabolic state are generated. This function is used in the grid search fucntion.

        Returns
        ----------
        flux list : List of dists of metabolic state data (number_of_initial_states > 1) .
        flux : Dist of metabolic state data (number_of_initial_states == 1) .
        state : State of generation results

        Examples
        --------
        >>> state, dict = model.generate_initial_states(50,1)
        >>> state, dict = model.generate_initial_states(50, 2, template = flux)

        See Also
        --------

        Fitting the metabolic model to multiple experiments


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 0:
                print('No experiment was set to the modelnames.')
            return False
        #
        # Mas=number of MKL thread control
        #
        try:
            import mkl
            mkl.set_num_threads(1)
        except:
            if self.configuration['callbacklevel'] > 2:
                print("mkl-service is not installed this python!")
        #
        # Set parameters
        #
        numbers = copy.deepcopy(self.numbers)
        vectors = copy.deepcopy(self.vector)
        configuration = copy.deepcopy(self.configuration)

        matrixinv = self.matrixinv
        initial_search_iteration_max = configuration["initial_search_iteration_max"]

        ub = [self.reactions[x]['ub'] for x in self.reaction_ids]
        lb = [self.reactions[x]['lb'] for x in self.reaction_ids]
        fluxes = []
        rsses = []

        if len(template) > 0:
            template = [template[type][id]["value"] for (type, id) in self.vector["ids"]]
        if method == "parallel":
            #for i in range(iterations):
            #    tmp_r, Rm_temp, Rm_ind, state = optimize.initializing_Rm_fitting(numbers, vectors, matrixinv, template ,initial_search_iteration_max)
            #    if state == "Determined":
            #        flux_temp_r = self.generate_state_dict(tmp_r)
            #        fluxes.append(flux_temp_r)
            #        rsses.append(self.calc_rss(flux_temp_r))
            #Set ncpus
            if 'ncpus' in self.configuration:
                ncpus = self.configuration['ncpus']
            else:
                ncpus = 1

            #Set callbacklevel
            if 'ppservers' in self.configuration:
                ppservers = self.configuration['ppservers']
            else:
                ppservers = ("",)
            #
            # tuple of all parallel python servers to connect with
            #
            try:
                import pp
                job_server = pp.Server(ncpus = ncpus, ppservers=ppservers)
            except:
                print("This function requires Parallel Python!")
                return False

            if (self.configuration['callbacklevel'] >= 6):
                print("Number of active nodes:", job_server.get_active_nodes())

            jobs = []


            for i in range(iterations):
                parameters = (numbers, vectors, matrixinv, template ,initial_search_iteration_max)
                jobs.append([i, job_server.submit(optimize.initializing_Rm_fitting, parameters,
                 (optimize.calc_protrude_scipy,),
                 ("numpy","scipy.optimize","mkl"))])


            for j, job in jobs:
                results = job()
                if results == None:
                    continue
                tmp_r, Rm_temp, Rm_ind, state = results
                if state == "Determined":
                    flux_temp_r = self.generate_state_dict(tmp_r)
                    fluxes.append(flux_temp_r)
                    rsses.append(self.calc_rss(flux_temp_r))
            if (self.configuration['callbacklevel'] >= 6):
                job_server.print_stats()
            job_server.destroy()

        else:
            for i in range(iterations):
                #print(template)
                tmp_r, Rm_temp, Rm_ind, state = optimize.initializing_Rm_fitting(numbers, vectors, matrixinv, template ,initial_search_iteration_max)
                if state == "Determined":
                    flux_temp_r = self.generate_state_dict(tmp_r)
                    fluxes.append(flux_temp_r)
                    rsses.append(self.calc_rss(flux_temp_r))
                    #print(self.calc_rss(flux_temp_r))
        order = sorted(list(range(len(fluxes))), key=lambda x: rsses[x])
        if self.configuration['callbacklevel'] >= 6:
            print("initial states were generated.", rsses)
        fluxes = [fluxes[x] for x in order]

        state = int(len(fluxes))
        if self.configuration['callbacklevel'] >= 6:
            print(state, "initial states were generated.")
        if len(fluxes) == 0:
            return state, fluxes
        else:
            fluxes = fluxes[:initial_states]

        if initial_states == 1:
            return state, fluxes[0]
        else:
            return state, fluxes


    def fitting_flux(self, method = 'SLSQP', flux = [], output = 'result'):
        """
        Fitting the metabolic model to multiple experiments

        Parameters
        ----------
        method:
            'SLSQP': A sequential least squares programming algorithm implemented by scipy
            "LN_PRAXIS": Gradient-free local optimization via the "principal-axis method" implemented by nlopt
            "GN_CRS2_LM": Controlled random searchimplemented for global optimizatoin by nlopt
            'deep': Repeated execution of SLSQP & LN_PRAXIS
        flux: Initial flux distribution generated by self.generate_initial_states().
        When flux is a dict of one state, single fitting trial is exected.
        When flux is a array of multiple dists of states, muptiple fitting trial is exected by using the parallel python.

        output: Output method. (default, output = 'result')
            'result': Fitting problem is solved inside of the function.
            'for_parallel': Fitting problems is NOT solved inside of the function. In this mode, a tupple of
            parameters for fit_r_mdv_pyopt() functions are generated for parallel computing.

        Parameters in self.configuration:
            iteration_max: Maximal number of interation in the optimizers. Example: self.set_configuration(iteration_max = 1000)
            number_of_repeat: Number of repeated execution by 'deep' functions. model.set_configuration(number_of_repeat = 2)

        Examples
        --------
        #
        # Single fitting trial
        #
        >>> state, flux_initial = model.generate_initial_states(10, 1)
        >>> state, RSS_bestfit, flux_opt = model.fitting_flux(method = 'deep', flux = flux_initial)
        >>> results= [("deep", flux_opt)]
        >>> model.show_results(results, pool_size = "off")
        #
        # Multipe fitting trials by using parallel python
        #
        >>> state, flux_initial = model.generate_initial_states(50, 4)
        >>> state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = "SLSQP", flux = flux_initial)
        >>> results= [("deep", flux_opt[0])]
        >>> model.show_results(results, pool_size = "off")
        #
        # Generating parameters for parallel python
        #
        >>> parameters = model.fitting_flux(method = "SLSQP", flux = flux_initial, output = "for_parallel")

        See Also
        --------
        generate_initial_states()


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 0:
                print('No experiment was set to the modelnames.')
            return False
        if flux == []:
            if self.configuration['callbacklevel'] >= 0:
                print('Metabolic state data is required.')
            return False
        #
        # Set parameters
        #
        numbers = copy.deepcopy(self.numbers)
        vectors = copy.deepcopy(self.vector)
        configuration = copy.deepcopy(self.configuration)
        matrixinv = self.matrixinv
        #
        #
        #
        r_depended = list(self.r_depended)
        #
        #
        #
        if output == "for_parallel":
            if method == "SLSQP":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "SLSQP")
            elif method == "LN_PRAXIS":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "LN_PRAXIS")
            elif method == "GN_CRS2_LM":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "GN_CRS2_LM")
            elif method == "deep":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux)
            else:
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "SLSQP")
            return parameters



        if isinstance(flux, dict):
            if method == "SLSQP":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "SLSQP")
            elif method == "COBYLA":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "COBYLA")
            elif method == "LN_COBYLA":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_COBYLA")
            elif method == "LN_BOBYQA":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_BOBYQA")
            elif method == "LN_NEWUOA":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_NEWUOA")
            elif method == "LN_PRAXIS":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_PRAXIS")
            elif method == "LN_SBPLX":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_SBPLX")
            elif method == "LN_NELDERMEAD":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "LN_NELDERMEAD")
            elif method == "GN_DIRECT_L":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "GN_DIRECT_L")
            elif method == "GN_CRS2_LM":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "GN_CRS2_LM")
            elif method == "GN_ESCH":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "GN_ESCH")
            elif method == "GN_IRES":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "GN_IRES")


            elif method == "deep":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_deep(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux)
            else:
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "SLSQP")
            return state, kai, self.generate_state_dict(opt_flux)


        else:
            #Set callbacklevel
            if 'ncpus' in self.configuration:
                ncpus = self.configuration['ncpus']
            else:
                ncpus = 1

            #Set callbacklevel
            if 'ppservers' in self.configuration:
                ppservers = self.configuration['ppservers']
            else:
                ppservers = ("",)
            #
            # tuple of all parallel python servers to connect with
            #
            try:
                import pp
                job_server = pp.Server(ncpus = ncpus, ppservers=ppservers)
            except:
                print("This function requires Parallel Python!")
                return False

            if (self.configuration['callbacklevel'] >= 4):
                print("Number of active nodes by pp", job_server.get_active_nodes())

            jobs = []

            if method == "SLSQP":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, "SLSQP")
                    jobs.append([i, job_server.submit(optimize.fit_r_mdv_scipy, parameters,
                     (optimize.calc_MDV_residue_scipy,),
                     ("numpy","scipy","scipy.integrate"))])
            elif method == "LN_PRAXIS":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, "LN_PRAXIS")
                    jobs.append([i, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                     (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                     ("numpy","nlopt","scipy","scipy.integrate"))])
            elif method == "GN_CRS2_LM":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, "GN_CRS2_LM")
                    jobs.append([i, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                     (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                     ("numpy","nlopt","scipy","scipy.integrate"))])
            elif method == "deep":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp)
                    jobs.append([i, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                     (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                     ("numpy","nlopt","scipy","scipy.integrate"))])
            else:
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp)
                    jobs.append([i, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                     (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                     ("numpy","nlopt","scipy","scipy.integrate"))])

            state_list = []
            kai_list = []
            flux_list = []
            Rm_ind_sol_list = []
            for j, job in jobs:
                results = job()
                if results == None:
                    continue
                state, rss, flux, Rm_ind_sol = results
                if len(flux) == 0:
                    continue
                state_list.append(state)
                kai_list.append(rss)
                flux_list.append(self.generate_state_dict(flux))
                Rm_ind_sol_list.append(Rm_ind_sol)
                if (self.configuration['callbacklevel'] >= 3):
                    print('RSS', j,':', rss, state)
            if (self.configuration['callbacklevel'] >= 4):
                job_server.print_stats()
            job_server.destroy()
            order = list(range(len(state_list)))
            order.sort(key = lambda x: kai_list[x])
            states = [state_list[x] for x in order]
            kais = [kai_list[x] for x in order]
            fluxes = [flux_list[x] for x in order]
            return states, kais, fluxes
    def pp_pretreatment(self):
        """
        Pitch a fitting_flux job to parallel python

        Parameters
        ----------

        Examples
        --------


        See Also
        --------



        """

        #Set callbacklevel
        if 'ncpus' in self.configuration:
            ncpus = self.configuration['ncpus']
        else:
            ncpus = 1

        #Set callbacklevel
        if 'ppservers' in self.configuration:
            ppservers = self.configuration['ppservers']
        else:
            ppservers = ("",)
        #
        # tuple of all parallel python servers to connect with
        #
        try:
            import pp
            job_server = pp.Server(ncpus = ncpus, ppservers=ppservers)
        except:
            print("This function requires Parallel Python!")
            return False

        if (self.configuration['callbacklevel'] >= 4):
            print("Number of active nodes by pp", job_server.get_active_nodes())

        return job_server


    def pitch_fitting_flux_job(self, job_server, method = 'SLSQP', flux = [], jobs = [], label = ""):
        """
        Pitch a fitting_flux job to parallel python

        Parameters
        ----------
        method:
            'SLSQP': A sequential least squares programming algorithm implemented by scipy
            "LN_PRAXIS": Gradient-free local optimization via the "principal-axis method" implemented by nlopt
            "GN_CRS2_LM": Controlled random searchimplemented for global optimizatoin by nlopt
            'deep': Repeated execution of SLSQP & LN_PRAXIS
        flux: Initial flux distribution generated by self.generate_initial_states().
        When flux is a dict of one state, single fitting trial is exected.
        When flux is a array of multiple dists of states, muptiple fitting trial is exected by using the parallel python.

        output: Output method. (default, output = 'result')
            'result': Fitting problem is solved inside of the function.
            'for_parallel': Fitting problems is NOT solved inside of the function. In this mode, a tupple of
            parameters for fit_r_mdv_pyopt() functions are generated for parallel computing.

        Parameters in self.configuration:
            iteration_max: Maximal number of interation in the optimizers. Example: self.set_configuration(iteration_max = 1000)
            number_of_repeat: Number of repeated execution by 'deep' functions. model.set_configuration(number_of_repeat = 2)

        Examples
        --------
        #
        # Single fitting trial
        #
        >>> state, flux_initial = model.generate_initial_states(10, 1)
        >>> state, RSS_bestfit, flux_opt = model.fitting_flux(method = 'deep', flux = flux_initial)
        >>> results= [("deep", flux_opt)]
        >>> model.show_results(results, pool_size = "off")
        #
        # Multipe fitting trials by using parallel python
        #
        >>> state, flux_initial = model.generate_initial_states(50, 4)
        >>> state, RSS_bestfit, flux_opt_slsqp = model.fitting_flux(method = "SLSQP", flux = flux_initial)
        >>> results= [("deep", flux_opt[0])]
        >>> model.show_results(results, pool_size = "off")
        #
        # Generating parameters for parallel python
        #
        >>> parameters = model.fitting_flux(method = "SLSQP", flux = flux_initial, output = "for_parallel")

        See Also
        --------
        generate_initial_states()


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 0:
                print('No experiment was set to the modelnames.')
            return False
        if flux == []:
            if self.configuration['callbacklevel'] >= 0:
                print('Metabolic state data is required.')
            return False



        if isinstance(flux, dict):
            flux = [flux]

        if method == "SLSQP":
            for i, flux_temp in enumerate(flux):
                parameters = self.fitting_flux(method = 'SLSQP', flux = flux_temp, output = 'for_parallel')
                jobs.append([label, job_server.submit(optimize.fit_r_mdv_scipy, parameters,
                 (optimize.calc_MDV_residue_scipy,),
                 ("numpy","scipy","scipy.integrate"))])
        elif method == "LN_PRAXIS":
            for i, flux_temp in enumerate(flux):
                parameters = self.fitting_flux(method = 'LN_PRAXIS', flux = flux_temp, output = 'for_parallel')
                jobs.append([label, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                 (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                 ("numpy","nlopt","scipy","scipy.integrate"))])
        elif method == "GN_CRS2_LM":
            for i, flux_temp in enumerate(flux):
                parameters = self.fitting_flux(method = 'GN_CRS2_LM', flux = flux_temp, output = 'for_parallel')
                jobs.append([label, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                 (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                 ("numpy","nlopt","scipy","scipy.integrate"))])
        elif method == "deep":
            for i, flux_temp in enumerate(flux):
                parameters = self.fitting_flux(method = 'deep', flux =flux_temp, output = 'for_parallel')
                jobs.append([label, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                 (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                 ("numpy","nlopt","scipy","scipy.integrate"))])
        else:
            for i, flux_temp in enumerate(flux):
                parameters = self.fitting_flux(method = 'deep', flux = flux_temp, output = 'for_parallel')
                jobs.append([label, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                 (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                 ("numpy","nlopt","scipy","scipy.integrate"))])



    def pp_posttreatment(self, job_server, jobs):
        """
        Pitch a fitting_flux job to parallel python

        Parameters
        ----------

        Examples
        --------


        See Also
        --------



        """
        label_list = []
        state_list = []
        kai_list = []
        flux_list = []
        Rm_ind_sol_list = []
        for label, job in jobs:
            results = job()
            if results == None:
                continue
            state, rss, flux, Rm_ind_sol = results
            if len(flux) == 0:
                continue
            label_list.append(label)
            state_list.append(state)
            kai_list.append(rss)
            flux_list.append(self.generate_state_dict(flux))
            Rm_ind_sol_list.append(Rm_ind_sol)
            if (self.configuration['callbacklevel'] >= 3):
                print('RSS', j,':', rss, state)
        if (self.configuration['callbacklevel'] >= 4):
            job_server.print_stats()
        job_server.destroy()
        order = list(range(len(state_list)))
        order.sort(key = lambda x: kai_list[x])
        states = [state_list[x] for x in order]
        kais = [kai_list[x] for x in order]
        fluxes = [flux_list[x] for x in order]
        lables = [label_list[x] for x in order]
        return lables, states, kais, fluxes


    def show_results(self, input, flux = "on", rss = "on", mdv = "on", pool_size = "on", filename = "", format = "csv"):
        """
        List of reactions in the model with metabolic flux data.

        Examples
        --------
        >>> free_reversible_reations = get_free_reversible_reactions()

        Parameters
        ----------
        input: Tuple of ("name", dic of metabolic statet) [('name1', state1),('name2', state2),('name3', state3),]
        flux: show flux data "on"/"off"
        rss: show rss data "on"/"off"
        mdv: show mdv data "on"/"off"
        pool_size: show metabolite data "on"/"off"
        filename: Results are saved when file name is not ""/
        format: "csv" or "text"

        Returns
        --------



        """
        reaction_header = []
        rssd = []
        reaction = []
        mdv_header = []
        mdvd = []
        metabolites_header = []
        metabolites = []
        #
        # Reaction header
        #
        reaction_header.extend(['Id', 'Reaction'])
        reaction_header.extend([z[0] for z in input])
        reaction_header.extend(['Atom_mapping','lb','ub','Used','Type','Value', 'Stdev', 'Reversible'])
        #print reaction_header
        #
        # Metabolite header header
        #
        metabolites_header.extend(['Id', 'C_number'])
        metabolites_header.extend([z[0] for z in input])
        metabolites_header.extend(['Excreted','lb','ub','Symmetry','Type','Value', 'Stdev', 'carbonsource'])
        #
        # RSS
        #
        rssd.append([])
        rssd[0].extend(['RSS', ''])
        for fluxd in input:
            rssd[0].append(self.calc_rss(fluxd[1]))
        rssd.append([])
        rssd.append([])
        rssd[1].extend(['Thres', ''])
        rssd[2].extend(['p_value', ''])
        for fluxd in input:
            pvalue, rss_thres = self.goodness_of_fit(fluxd[1], alpha = 0.05)
            rssd[1].append(rss_thres)
            rssd[2].append(pvalue)
        #
        #metabolites
        #
        for i, metid in enumerate(self.metabolite_ids):
            metabolites.append([])
            metabolites[i].append(metid)# 1 id
            metabolites[i].append(self.metabolites[metid]['C_number'])# 2 reaction
            metabolites[i].extend([fluxdd["metabolite"][metid]['value'] for (name, fluxdd) in input])# 4 values
            metabolites[i].append(self.metabolites[metid]['excreted'])# 5 atom mapping
            metabolites[i].append(self.metabolites[metid]['lb'])# 6 lb
            metabolites[i].append(self.metabolites[metid]['ub'])# 7 ub
            metabolites[i].append(self.metabolites[metid]['symmetry'])# 8 use or not
            metabolites[i].append(self.metabolites[metid]['type'])# 9 type
            metabolites[i].append(self.metabolites[metid]['value'])# 10 met value (measured)
            metabolites[i].append(self.metabolites[metid]['stdev'])# 11 met std
            metabolites[i].append(self.metabolites[metid]['carbonsource'])# 12 revserible#id
        #
        # Reaction and fluxes
        #
        for i, rid in enumerate(self.reaction_ids):
            reaction.append([])
            reaction[i].append(rid)# 1 id
            reaction[i].append(self.reactions[rid]['stoichiometry'])# 2 reaction
            reaction[i].extend([fluxdd["reaction"][rid]['value'] for (name, fluxdd) in input])# 4 values
            reaction[i].append(self.reactions[rid]['atommap'])# 5 atom mapping
            reaction[i].append(self.reactions[rid]['lb'])# 6 lb
            reaction[i].append(self.reactions[rid]['ub'])# 7 ub
            reaction[i].append(self.reactions[rid]['use'])# 8 use or not
            reaction[i].append(self.reactions[rid]['type'])# 9 type
            reaction[i].append(self.reactions[rid]['value'])# 10 flux value
            reaction[i].append(self.reactions[rid]['stdev'])# 11 flux std
            reaction[i].append(self.reactions[rid]['reversible'])# 12 revserible#id
        reversible_list = self.reversible.keys()
        length = len(reaction)
        for i, rid in enumerate(self.reversible.keys()):
            forward = self.reversible[rid]['forward']
            reverse = self.reversible[rid]['reverse']
            reaction.append([])
            reaction[length + i].append(rid)# 1 id
            reaction[length + i].append(str(forward + "<=>" + reverse))# 2 reaction
            reaction[length + i].extend([fluxdd["reversible"][rid]['value'] for (name, fluxdd) in input])# 4 values
            reaction[length + i].append('')# 5 atom mapping
            reaction[length + i].append(self.reversible[rid]['lb'])# 6 lb
            reaction[length + i].append(self.reversible[rid]['ub'])# 7 ub
            reaction[length + i].append('')# 8 use or not
            reaction[length + i].append(self.reversible[rid]['type'])# 9 type
            reaction[length + i].append(self.reversible[rid]['value'])# 10 flux value
            reaction[length + i].append(self.reversible[rid]['stdev'])# 11 flux std
            reaction[length + i].append('')# 12 revserible#id
        #print reaction
        #
        #
        # MDVs header
        #
        mdv_header.extend(['Experiment', 'Fragment_Num',"Time"])
        mdv_header.extend([fluxd[0] for fluxd in input])
        mdv_header.extend(['Use','Ratio','Stdev'])
        #
        # MDVs
        #
        for ex_id in sorted(self.experiments.keys()):
            mdv_data = []
            target_fragments_temp = self.target_fragments.keys()
            mdv_carbon_sources_temp = self.experiments[ex_id]['mdv_carbon_sources']
            if self.experiments[ex_id]['mode'] == "ST":
                for (name, fluxd) in input:
                    tmp_r = [fluxd["reaction"][x]['value'] for x in self.reaction_ids]
                    mdv_exp, mdv_hash = optimize.calc_MDV_from_flux(tmp_r, target_fragments_temp, mdv_carbon_sources_temp, self.func)
                    mdv_data.append((name, mdv_hash))
                #if len(self.experiments[ex_id]["timepoint"]) == 0:
                for i in range(len(self.experiments[ex_id]['mdv_ids'])):
                    mdvd.append([])
                    mdvd[len(mdvd)-1].append(ex_id) # 1 ex
                    fragment_number = self.experiments[ex_id]['mdv_ids'][i]
                    mdvd[len(mdvd)-1].append(fragment_number) # 2 fragment
                    mdvd[len(mdvd)-1].append(int(0)) # 3 Timecourse
                    number = int(fragment_number.split("_")[-1])
                    fragment = "_".join(fragment_number.split("_")[:-1])
                    mdvd[len(mdvd)-1].extend([mdvdd[1][fragment][number] for mdvdd in mdv_data])
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_use'][i]) # 6 use or not
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_exp_original'][i]) # 7 ratio
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_std_original'][i]) # 7 stdev

            elif self.experiments[ex_id]['mode'] == "INST":
                timepoints = self.experiments[ex_id]['timepoint']
                timepoint_id={}
                for n,point in enumerate(timepoints):
                    timepoint_id[point]=n
                for (name, fluxd) in input:
                    tmp_r = [fluxd[group][id]["value"] for (group, id) in self.vector["ids"]]
                    mdv_exp, mdv_hash = optimize.calc_MDV_from_flux(tmp_r, target_fragments_temp, mdv_carbon_sources_temp, self.func, timepoint = timepoints, y0temp = self.experiments[ex_id]['y0'])
                    mdv_data.append((name, mdv_hash))
                timepoint_list= sorted(timepoints*int(len(self.experiments[ex_id]['mdv_ids'])/len(timepoints)))

                #timepoint_list= sorted(timepoints*(len(self.experiments[ex_id]['mdv_ids'])))
                for j in range(len(self.experiments[ex_id]['mdv_ids'])):
                    mdvd.append([])
                    mdvd[len(mdvd)-1].append(ex_id) # 1 ex
                    fragment_number = self.experiments[ex_id]['mdv_ids'][j]
                    number = int(fragment_number.split("_")[-1])
                    fragment = "_".join(fragment_number.split("_")[:-1])
                    mdvd[len(mdvd)-1].append(fragment_number) # 2 fragment
                    mdvd[len(mdvd)-1].append(timepoint_list[j])   # 3 timepoint
                    mdvd[len(mdvd)-1].extend([mdvdd[1][fragment][timepoint_id[timepoint_list[j]]][number] for mdvdd in mdv_data])
                    mdvd[len(mdvd)-1].append(ex_id) # 1 ex
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_use'][j]) # 6 use or not
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_exp_original'][j]) # 7 ratio
                    mdvd[len(mdvd)-1].append(self.experiments[ex_id]['mdv_std_original'][j]) # 7 stdev

        #
        #
        # output results
        #
        #
        if len(filename) > 0:
            import csv
            if format == "csv":
                with open(filename, 'w', newline='') as f:
                    writer = csv.writer(f, dialect='excel')
                    #writer = csv.writer(f)
                    if rss == "on":
                        writer.writerow(reaction_header)
                        writer.writerows(rssd)
                    if flux == "on":
                        writer.writerow(reaction_header)
                        writer.writerows(reaction)
                    if pool_size == "on":
                        writer.writerow(metabolites_header)
                        writer.writerows(metabolites)
                    if mdv == "on":
                        writer.writerow(mdv_header)
                        writer.writerows(mdvd)
                f.close()
            if format == "text":
                with open(filename, 'w', newline='') as f:
                    writer = csv.writer(f, delimiter='\t')
                    if rss == "on":
                        writer.writerow(reaction_header)
                        writer.writerows(rssd)
                    if flux == "on":
                        writer.writerow(reaction_header)
                        writer.writerows(reaction)
                    if pool_size == "on":
                        writer.writerow(metabolites_header)
                        writer.writerows(metabolites)
                    if mdv == "on":
                        writer.writerow(mdv_header)
                        writer.writerows(mdvd)
        else:
            text = ""
            if rss == "on":
                #print reaction header
                text = text + "{0:15.15s}".format(reaction_header[0])
                text = text + "{0:25.25s}".format(reaction_header[1])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(reaction_header[2+i]))
                text = text + "\n"
                for data in rssd:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.25s}".format(data[1])
                    for i in range(len(input)):
                        text = text + "{0:>8.2f}".format(data[2+i])
                    text = text + "\n"

            if flux == "on":
                #print flux data
                text = text + "{0:15.15s}".format(reaction_header[0])
                text = text + "{0:25.25s}".format(reaction_header[1])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(reaction_header[2+i]))
                text = text + " " + "{0:25.25s}".format(reaction_header[-8])
                text = text + "{0:>6.5s}".format(reaction_header[-7])
                text = text + "{0:>6.5s}".format(reaction_header[-6])
                text = text + " " + "{0:5.4s}".format(reaction_header[-5])
                text = text + "{0:8.7s}".format(reaction_header[-4])
                text = text + "{0:>6.5s}".format(reaction_header[-3])
                text = text + "{0:>6.5s}".format(reaction_header[-2])
                text = text + "\n"
                for data in reaction:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.25s}".format(data[1])
                    for i in range(len(input)):
                        text = text + "{0:>8.1f}".format(data[2+i])
                    text = text + " " + "{0:25.25s}".format(data[-8])
                    text = text + "{0:6.1f}".format(float(data[-7]))
                    text = text + "{0:6.1f}".format(float(data[-6]))
                    text = text + " " + "{0:5.4s}".format(data[-5])
                    text = text + "{0:8.7s}".format(data[-4])
                    text = text + "{0:6.1f}".format(data[-3])
                    text = text + "{0:6.1f}".format(data[-2])
                    text = text + "\n"

            if pool_size == "on":
                #print reaction header
                text = text + "{0:15.15s}".format(metabolites_header[0])
                text = text + "{0:5.5s}".format(metabolites_header[1])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(metabolites_header[2+i]))
                text = text + " "+"{0:6.6s}".format(metabolites_header[-7])
                text = text + " "+"{0:6.6s}".format(metabolites_header[-6])
                text = text + " "+"{0:9.9s}".format(metabolites_header[-5])
                text = text + "{0:7.7s}".format(metabolites_header[-4])
                text = text + "{0:>6.6s}".format(metabolites_header[-3])
                text = text + "{0:>6.6s}".format(metabolites_header[-2])
                text = text + "\n"
                for data in metabolites:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:5.0f}".format(data[1])
                    for i in range(len(input)):
                        text = text + "{0:>8.3f}".format(data[2+i])
                    text = text + " "+"{0:6.3f}".format(float(data[-7]))
                    text = text + " "+"{0:6.1f}".format(float(data[-6]))
                    text = text + " "+"{0:9.9s}".format(data[-5])
                    text = text + "{0:7.7s}".format(data[-4])
                    text = text + "{0:6.1f}".format(data[-3])
                    text = text + "{0:6.1f}".format(data[-2])
                    text = text + "\n"

            if mdv == "on":
                text = text + "{0:15.15s}".format(mdv_header[0])
                text = text + "{0:15.10s}".format(mdv_header[1])
                text = text + "{0:5.5s}".format(mdv_header[2])
                for i in range(len(input)):
                    text = text + "{0:8.7s}".format(str(mdv_header[3+i]))
                text = text + "{0:>6.5s}".format(str(mdv_header[-3]))
                text = text + " " + "{0:>6.6s}".format(mdv_header[-2])
                text = text + "{0:>6.6s}".format(mdv_header[-1])
                text = text + "\n"
                for data in mdvd:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:15.15s}".format(data[1])
                    text = text + "{0:5.2f}".format(data[2])
                    for i in range(len(input)):
                        text = text + "{0:8.4f}".format(data[3+i])
                    text = text + "{0:>6.5s}".format(str(data[-3]))
                    text = text + " " + "{0:6.4f}".format(data[-2])
                    text = text + "{0:6.2f}".format(data[-1])
                    text = text + "\n"
            print(text)

    def calc_rss(self, flux, mode = "flux"):
        """
        Determine a residual sum of square (RSS) between estimated MDVs of 'flux' and measured MDVs in 'experiment(s)'

        Parameters
        ----------
        flux: Dictionary of metabolic flux distribution or independent flux vector

        mode: "independent" Calculation of RSS from independent flux vector


        Returns
        ----------
        rss : Residual sum of square (RSS)

        Examples
        --------
        >>> rss = model.calc_rss(flux_opt)
        >>> print rss
        7564.123

        See Also
        --------


        History
        --------
        140912 Call of calmdv funcion is removed.


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 0:
                print('No experiment was set to the modelnames.')
            return False
        if flux == []:
            if self.configuration['callbacklevel'] >= 0:
                print('Metabolic state data is required.')
            return False

        # zero independent flux
        if mode == "independent":
            Rm_ind = flux
        else:
            Rm_ind = [flux[group][id]["value"] for (group, id) in self.vector['independent_flux']]
        #
        # MDV vector of all experiments
        #

        mdv_exp_original = list(self.vector["value"])
        mdv_std_original = list(self.vector["stdev"])
        mdv_use = list(self.vector["use"])
        for experiment in sorted(self.experiments.keys()):
            mdv_exp_original.extend(self.experiments[experiment]['mdv_exp_original'])
            mdv_std_original.extend(self.experiments[experiment]['mdv_std_original'])
            mdv_use.extend(self.experiments[experiment]['mdv_use'])
        mdv_exp = numpy.array([y for x, y in enumerate(mdv_exp_original) if mdv_use[x] != 0])
        spectrum_std = numpy.array([y for x, y in enumerate(mdv_std_original) if mdv_use[x] != 0])
        #
        # Covariance matrix
        #
        covinv = numpy.zeros((len(spectrum_std),len(spectrum_std)))
        for i, std in enumerate(spectrum_std):
            if std <= 0.0:
                print("Error in ", i, std)
            covinv[i,i] = 1.0/(std**2)

        rss = optimize.calc_MDV_residue(Rm_ind,
            stoichiometric_num = self.numbers['independent_start'],
            reaction_num= self.numbers['total_number'],
            matrixinv=self.matrixinv,
            experiments=self.experiments,
            mdv_exp=mdv_exp,
            mdv_use=mdv_use,
            covinv=covinv,
            Rm_initial= self.vector["Rm_initial"],
            lb = copy.copy(self.vector["lb"]),
            ub = copy.copy(self.vector["ub"]),
            reac_met_number = self.numbers['reac_met_number'],
            calmdv = self.func["calmdv"],
            diffmdv = self.func["diffmdv"]
            )
        return rss

    def goodness_of_fit(self, flux, alpha = 0.05):
        """
        Calculate goodness-of-fit of a given flux distribution

        Examples
        --------
        >>> p-value = model.goodness_of_fit()

        Parameters
        ----------
        flux: Dictionary of metabolic flux distribution
        alpha: Confidence level. Default is 0.05.

        Returns
        --------
        pvalue: p-value determined by chi-squere test.
        rss_thres: Threshold RSS at the alpha levels

        Examples
        --------
        >>> pvalue, rss_thres= model.goodness_of_fit(flux_opt_fitted, alpha = 0.05)


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 0:
                print('No experiment was set to the modelnames.')
            return False
        if flux == []:
            if self.configuration['callbacklevel'] >= 0:
                print('Metabolic state data is required.')
            return False

        from scipy.stats import chi2, f

        degree_of_freedom = self.get_degree_of_freedom()*1.0
        number_of_measurement = self.get_number_of_independent_measurements()*1.0
        RSS = self.calc_rss(flux)*1.0
        thres = chi2.ppf(1.0-alpha, number_of_measurement - degree_of_freedom)
        pvalue = chi2.sf(RSS, number_of_measurement - degree_of_freedom)

        if self.configuration["callbacklevel"] >= 8:
            print("number_of_measurement", number_of_measurement)
            print("degree_of_freedom", degree_of_freedom)
            print("RSS", RSS)
            print("thres", thres)
            print("pvalue", pvalue)
        if self.configuration["callbacklevel"] >= 1:
            if number_of_measurement <= degree_of_freedom:
                print("number_of_measurement", number_of_measurement, "is smaller than degree_of_freedom", degree_of_freedom)
        return pvalue, thres


    def get_thres_confidence_interval(self, flux, alpha = 0.05, dist = 'F_dist'):
        """
        Determine a threshold level of RSS for searching confidence interval

        Parameters
        ----------
        flux: Dictionary of best fitted flux distribution.
        alpha: confidence level. For 95% confience interval, please set alpha = 0.05
        dist:
            'F-dist' (default): Using F-distribution
            'chai-dist':Using Chai-square distribution

        Returns
        --------
        thres: a threshold level of RSS for searching confidence interval
        number_of_measurements: number of measurements
        degree_of_freedom: degree of freecom"

        Examples
        --------
        >>> thres, number_of_reactions, degree_of_freedom = model.get_thres_confidence_interval(RSS, alpha = 0.05, dist = 'F_dist')

        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 1:
                print('No experiment was set to the modelnames.')
            return False
        RSS = self.calc_rss(flux)
        from scipy.stats import chi2, f
        n = self.get_number_of_independent_measurements()
        p = self.get_degree_of_freedom()

        if dist == 'chai-dist':
            return RSS + chi2.ppf(1-alpha, 1), n, p
        return RSS * (1 + f.ppf(1-alpha, 1, n-p)/(n-p)), n, p


    def get_degree_of_freedom(self):
        """
        Getter of degree of freedom of the model

        Examples
        --------
        >>> degree_of_freedom = model.get_degree_of_freedom()

        Returns
        --------
        degree_of_freedom: degree of freedom.

        """
        return self.numbers["independent_number"]

    def get_number_of_independent_measurements(self):
        """
        Getter of number of independent measurements of the model

        Examples
        --------
        >>> number_of_independent_measurements = model.get_number_of_independent_measurements()

        Returns
        --------
        number_of_independent_measurements: number of independent measurements.

        History
        --------
        29.9.2014 Modified to consider "fitting" flux information.

        """
        n = 0;
        for experiment in self.experiments:
            n = n + self.experiments[experiment]['number_of_measurement']
        #
        # Count "fitting" reactions
        #
        fitting_reactions = []
        for id in self.reactions:
            if self.reactions[id]['type'] == 'fitting':
                fitting_reactions.append(id)
        for id in self.metabolites:
            if self.metabolites[id]['type'] == 'fitting':
                fitting_reactions.append(id)
        for id in self.reversible:
            if self.reversible[id]['type'] == 'fitting':
                fitting_reactions.append(id)
        return n + len(fitting_reactions)

    def generate_ci_templete(self, targets = "normal"):
        """
        Generate a blank dictionary to keep confidence interval search results.
        A dictionaty generated by the


        Parameters
        ----------
        targets:
        'normal': Free reversible and inreversible reactions
        'independent': Independent reactions
        'all': all reactions
        'without_reversible_reactions': inreversible reactions


        Examples
        --------
        >>> ci = model.generate_ci_templete(targets = [("reaction","r29_g6pdh"),("reaction","r27_pc"),("reaction","r28_mae")])
        >>> ci = model.generate_ci_templete(targets = 'normal')

        See Also
        --------
        search_ci.


        """
        ci = {
        'record':{},
        'data':{}
        }
        #
        # Upper and lower boundaries
        #
        ub = self.vector["ub"]
        lb = self.vector["lb"]
        for (state, id) in self.vector['ids']:
            if state == "reaction":
                data = self.reactions[id]
            elif state == "metabolite":
                data = self.metabolites[id]
            elif state == "reversible":
                data = self.reversible[id]
            ci['data'][(state,id)] = {
            'forward': id,
            'reverse': '',
            'upper_boundary': ub[data['position_in_tmp_r']],
            'lower_boundary': lb[data['position_in_tmp_r']],
            'flux_data':[],
            'rss_data':[],
            'raw_flux_data':[],
            'upper_boundary_state':"not deternined",
            'lower_boundary_state':"not deternined",
            'score':'',
            'use':'no',
            'type': data['type'],
            }
            if state == "reversible":
                ci['data'][(state,id)]['forward'] = data['forward']
                ci['data'][(state,id)]['reverse'] = data['reverse']
            elif state == "reaction":
                ci['data'][(state,id)]['reversible'] = data['reversible']


        if targets == 'normal':
            for (state, id) in self.vector['ids']:
                if state == "reaction":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        if ci['data'][(state,id)]['reversible'] == "no":
                            ci['data'][(state,id)]['use'] = 'on'
                if state == "metabolite":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
                if state == "reversible":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
        elif targets == 'without_reversible_reactions':
                if state == "reaction":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        if ci['data'][(state,id)]['reversible'] == "no":
                            ci['data'][(state,id)]['use'] = 'on'
                if state == "metabolite":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
        elif targets == 'all':
            for (state, id) in self.vector['ids']:
                if state == "reaction":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
                if state == "metabolite":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
                if state == "reversible":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
        elif targets == 'independent':
            for (state, id) in self.vector["independent_flux"]:
                if state == "reaction":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        if ci['data'][(state,id)]['reversible'] == "no":
                            ci['data'][(state,id)]['use'] = 'on'
                if state == "metabolite":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
                if state == "reversible":
                    if ci['data'][(state,id)]['type'] != "fixed":
                        ci['data'][(state,id)]['use'] = 'on'
        else:
            for state, id in targets:
                ci['data'][(state,id)]['use'] = 'on'

        return ci



    def search_ci(self, ci, flux, method = 'grid', alpha = 0.05, dist = 'F_dist'):
        """
        Confidence intervals were determined based on information in 'ci'

        Parameters
        ----------
        ci: dictionary for storeing confidence interval results. ci is generated by self.generate_ci_templete()
        flux: Distribution of best fitted flux
        method:
            'grid': Grid search method
        alpha: (default 0.05) Confidence interval level. 0.05 means a 95% confidence interval level
        dist:
            'F-dist' (default): Using F-distribution
            'chai-dist':Using Chai-square distribution

        Examples
        --------
        >>> ci = model.generate_ci_templete(targets = 'independent')
        >>> ci = model.search_confidence_interval_parallel(ci, flux, method = 'grid')


        Returns
        --------
        ci:dictionary of confidence interval results.

        """

        #Set callbacklevel
        if 'callbacklevel' in self.configuration:
            callbacklevel = self.configuration['callbacklevel']
        else:
            callbacklevel = 0
        #
        # Set parameters
        #
        import copy
        numbers = copy.deepcopy(self.numbers)
        vectors = copy.deepcopy(self.vector)
        configuration = copy.deepcopy(self.configuration)
        matrixinv = self.matrixinv
        r_depended = list(self.r_depended)
        #
        #Set ncpus
        if 'ncpus' in self.configuration:
            ncpus = self.configuration['ncpus']
        else:
            ncpus = 1

        #Set ppservers
        if 'ppservers' in self.configuration:
            ppservers = self.configuration['ppservers']
        else:
            ppservers = ("",)
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            print('No experiment was set to the modelnames.')
            return False
        #
        # Check best fitted flux
        #
        if len(flux) == 0:
            print('No best fit flux was setted.')
            return False
        #
        # Calc threadfold
        #
        thres, number_of_measurements, degree_of_freedom  = self.get_thres_confidence_interval(flux, alpha = alpha, dist = dist)
        #
        # Check Parallel python
        #
        try:
            import pp
            job_server = pp.Server(ncpus = ncpus, ppservers=ppservers)
        except:
            print("This function requires Parallel Python!")
            return False

        if (callbacklevel >= 2):
            print("Number of active nodes is " + str(job_server.get_active_nodes()))
        ci['record']['flux'] = flux
        ci['record']['thres'] = thres
        ci['record']['rss'] = self.calc_rss(flux)

        jobs = []

        #######################################################
        #
        # Grid search
        #
        ########################################################
        if method == "grid":
            starttime=time.perf_counter()
            step = 0.0
            data_tmp={}
            job_number = self.configuration["grid_search_iterations"]
            initial_search_repeats_in_grid_search = self.configuration["initial_search_repeats_in_grid_search"]
            #
            # confidence interval
            #
            rss_bestfit = self.calc_rss(flux)
            jobs=[]
            for group, rid in sorted(ci['data'].keys()):
                if ci['data'][(group, rid)]['use'] != 'on': continue
                flux_opt_rid = flux[group][rid]['value']


                data={(group, rid):{"flux_data":[flux_opt_rid],"rss_data":[rss_bestfit],"raw_flux_data":[[flux["reaction"][i]["value"] for i in self.reaction_ids]]}}
                data_tmp.update(data)

                flux_lower =ci['data'][(group, rid)]["lower_boundary"]
                flux_upper =ci['data'][(group, rid)]["upper_boundary"]



                # Grid number is 20
                if (flux_upper-flux_lower)>=100:
                    n=15
                else:
                    n=10.0
                n=15

                if self.configuration['callbacklevel'] >= 1:
                    print("Setting:", rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower," ub ",flux_upper, "n=", n)
                #
                # Store original metabolic constrains
                #
                temp_type, temp_value, temp_stdev = self.get_constrain(group, rid)
                #
                # Set initial search to upward
                #
                temp_array_initial_fluxes = []
                counter_of_missed_initial_state = 0
                #
                for i in range(int(n+1)):
                    if  counter_of_missed_initial_state > 5:
                        break
                    step = (flux_upper-flux_opt_rid) * (0.5**(n-i))
                    fixed_flux = flux_opt_rid + step
                    if fixed_flux > flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux < flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    #


                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if self.configuration['callbacklevel'] >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)

                    for i in range(job_number):
                        #state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux)
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            #
                            # Set large data
                            #
                            data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                            data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                            data_tmp[(group, rid)]["raw_flux_data"].append([])
                            if self.configuration['callbacklevel'] >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                            counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        #parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                        #functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                        #jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                        #if self.configuration['callbacklevel'] >= 2:
                        #    print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                #
                # Set initial search to upward
                #
                counter_of_missed_initial_state = 0
                #
                for i in range(int(n+1)):
                    if  counter_of_missed_initial_state > 5:
                        break
                    step = (flux_opt_rid - flux_lower) * (0.5**(n-i))
                    fixed_flux = flux_opt_rid - step
                    if fixed_flux > flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux < flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if self.configuration['callbacklevel'] >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)                    #

                    for i in range(job_number):
                        #state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux)
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            #
                            # Set large data
                            #
                            data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                            data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                            data_tmp[(group, rid)]["raw_flux_data"].append([])
                            if self.configuration['callbacklevel'] >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                            counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        #parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                        #functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                        #jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                        #if self.configuration['callbacklevel'] >= 2:
                        #    print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in temp_array_initial_fluxes:
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if self.configuration['callbacklevel'] >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                    jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                    if self.configuration['callbacklevel'] >= 2:
                        print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)

                #
                # Return to original position
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()
                if self.configuration['callbacklevel'] >= 2:
                    print("Waiting for pp response", rid, "flux_opt ")
            #
            #
            # Retrive results
            #
            for (group, rid), flux_value, job in jobs:
                results = job()
                if results == None:
                    #continue
                    #
                    # Large value is used when falied
                    #
                    rss = thres * 10.0
                    opt_flux = []
                    state = "failed"
                else:
                    state, rss, opt_flux, Rm_ind_sol = results
                #
                #  Large value is used when falied
                #
                if len(opt_flux) == 0:
                    #continue
                    rss = thres * 10.0
                    opt_flux = []

                if self.configuration['callbacklevel'] >= 2:
                    print("Finished:", rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, ",n ", n)

                data_tmp[(group, rid)]["flux_data"].append(flux_value)
                data_tmp[(group, rid)]["rss_data"].append(rss)
                data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)

            #
            # Cals again around thres
            #
            jobs_right=[]
            jobs_left=[]
            for group, rid in ci['data'].keys():
                if ci['data'][(group, rid)]['use'] != 'on': continue
                forward_id = ci['data'][(group, rid)]['forward']
                reverse_id = ci['data'][(group, rid)]['reverse']
                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']
                # Detect just before the lower boundary
                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[1])
                #
                # Determinie lower bondary
                #

                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if i == 0:
                        rss_previous = rss * 1.0
                        flux_previous_left  = flux_temp * 1.0
                        continue
                    if flux_temp >= flux_previous_left -0.01:
                        continue
                    if rss > thres:
                        #flux_next_left = flux * 1.0
                        continue
                    rss_previous = rss * 1.0
                    flux_previous_left = flux_temp * 1.0

                flux_next_left = ci['data'][(group, rid)]["lower_boundary"]
                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if flux_temp >= flux_previous_left:
                        continue
                    if rss > thres:
                        flux_next_left = flux_temp * 1.0
                        break



                # Detect just before the upper boundary
                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[1])

                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if i == 0:
                        rss_previous = rss * 1.0
                        flux_previous_right = flux_temp * 1.0
                        continue
                    if flux_temp <= flux_previous_right+0.01:
                        continue
                    if rss > thres:
                        #flux_next_right = flux * 1.0
                        break
                    rss_previous = rss * 1.0
                    flux_previous_right = flux_temp * 1.0

                flux_next_right = ci['data'][(group, rid)]["upper_boundary"]
                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if flux_temp <= flux_previous_right:
                        continue
                    if rss > thres:
                        flux_next_right = flux_temp * 1.0
                        break

                if self.configuration['callbacklevel'] >= 2:
                    print("flux_previous_left ",flux_previous_left,", flux_previous_right ",flux_previous_right)
                    print("flux_next_left ",flux_next_left,", flux_next_right ",flux_next_right)

                if (flux_next_right-flux_next_left)>=40:
                    n=20
                elif (flux_next_right-flux_next_left)>=5:
                    n=10
                else:
                    n=5
                #
                # Store original metabolic constrains
                #
                temp_type, temp_value, temp_stdev = self.get_constrain(group, rid)
                #
                temp_array_initial_fluxes = []
                for e in range (1,n):
                    #
                    #
                    # right direcition
                    step = abs(flux_previous_right - flux_next_right)/n
                    fixed_flux = flux_previous_right +step*e
                    if fixed_flux > flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux < flux_lower:
                        fixed_flux = flux_lower

                    #
                    # fix reacton
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    for i in range(job_number):
                        #state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux)
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            #
                            # Set large data
                            #
                            data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                            data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                            data_tmp[(group, rid)]["raw_flux_data"].append([])

                            if self.configuration['callbacklevel'] >= 2:
                                print("Passed:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        #
                        #
                        #
                        #parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                        #functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                        #jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                        #if self.configuration['callbacklevel'] >= 2:
                        #    print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)


                    #
                    # left reaction
                    #
                    step = abs(flux_previous_left - flux_next_left)/n
                    fixed_flux = flux_previous_left - step*e
                    if fixed_flux > flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux < flux_lower:
                        fixed_flux = flux_lower
                    #
                    # fixed reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    for i in range(job_number):
                        #state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux)
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            #
                            # Set large data
                            #
                            data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                            data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                            data_tmp[(group, rid)]["raw_flux_data"].append([])
                            if self.configuration['callbacklevel'] >= 2:
                                print("Passed:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        #
                        #
                        #
                        #parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                        #functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                        #jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                        #if self.configuration['callbacklevel'] >= 2:
                        #    print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in temp_array_initial_fluxes:
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if self.configuration['callbacklevel'] >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                    jobs.append([(group, rid),fixed_flux , job_server.submit(optimize.fit_r_mdv_scipy, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                    if self.configuration['callbacklevel'] >= 2:
                        print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i)
                #
                # Return to original position [need modificaton]
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()

            #
            #
            # Retrive results
            #
            for (group, rid), flux_value, job in jobs:
                results = job()
                if results == None:
                    #continue
                    #
                    # Large value is used when falied
                    #
                    rss = thres * 10.0
                    opt_flux = []
                    state = "failed"
                else:
                    state, rss, opt_flux, Rm_ind_sol = results
                #
                #  Large value is used when falied
                #
                if len(opt_flux) == 0:
                    #continue
                    rss = thres * 10.0
                    opt_flux = []

                if self.configuration['callbacklevel'] >= 2:
                    print("Finished:", rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, ",n ", n)

                data_tmp[(group, rid)]["flux_data"].append(flux_value)
                data_tmp[(group, rid)]["rss_data"].append(rss)
                data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)
            if self.configuration['callbacklevel'] >= 1:
                print("Data collection was finished")

            for group, rid in ci['data'].keys():
                if ci['data'][(group, rid)]['use'] != 'on': continue
                forward_id = ci['data'][(group, rid)]['forward']
                reverse_id = ci['data'][(group, rid)]['reverse']
                #
                # Retirieve results
                #
                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']

                #
                # Detect lower bounary of confidence interval
                #
                state_lower = "not determined"
                flux_lower = ci['data'][(group, rid)]['lower_boundary']
                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[1])
                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if i == 0:
                        rss_previous = rss * 1.0
                        flux_previous = flux_temp * 1.0
                        continue
                    if flux_temp >= flux_previous-0.01:
                        continue
                    if rss > thres:
                        flux_lower = flux_previous - (flux_previous - flux_temp) * (thres - rss_previous)/(rss - rss_previous)
                        state_lower = "determined"
                        break
                    rss_previous = rss * 1.0
                    flux_previous = flux_temp * 1.0
                #
                # Detect upper bounary of confidence interval
                #
                state_upper = "not determined"
                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[1])
                for i, tuples in enumerate(data_sorted):
                    rss = tuples[1]
                    flux_temp = tuples[0]
                    if i == 0:
                        rss_previous = rss * 1.0
                        flux_previous = flux_temp * 1.0
                        continue
                    if flux_temp <= flux_previous+0.01:
                        continue
                    if rss > thres:
                        flux_upper = flux_previous + (flux_temp - flux_previous) * (thres - rss_previous)/(rss - rss_previous)
                        state_upper = "determined"
                        break
                    rss_previous = rss * 1.0
                    flux_previous = flux_temp * 1.0
                if self.configuration['callbacklevel'] >= 1:
                    print(rid, 'Grid search maeda method. Lower boundary:', flux_lower, 'Upper boundery' ,flux_upper)
                ci['data'][(group, rid)]['upper_boundary'] = flux_upper
                ci['data'][(group, rid)]['lower_boundary'] = flux_lower
                ci['data'][(group, rid)]['upper_boundary_state'] = state_upper
                ci['data'][(group, rid)]['lower_boundary_state'] = state_lower
                #
                # Grids with too large rss are removed
                #
                flux_values_array = [flux_values_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < (thres * 2) ]
                raw_flux_array = [raw_flux_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < thres * 2 ]
                rss_values_array = [rss_values_array_tmp[i] for (i, rsst)  in enumerate(rss_values_array_tmp) if rsst < thres * 2 ]

                ci['data'][(group, rid)]['flux_data'] = flux_values_array
                ci['data'][(group, rid)]['rss_data'] = rss_values_array
                ci['data'][(group, rid)]['raw_flux_data'] = raw_flux_array



            #"
            # Show calc time
            #
            if self.configuration['callbacklevel'] >= 1:
                endtime=time.perf_counter()
                print("time elapsed ", endtime-starttime,"s")

            #return flux_lb
        if (callbacklevel >= 0):
            job_server.print_stats()

        job_server.destroy()
        return ci
    def load_states(self, filename, format = 'text'):
        """
        Load a text/csv file with 'reacton type' information to generate new states dict.

        Parameters
        ----------
        model : a instance of mfapy.model for templete
        filename : filename of flux data with following format.

        type	Id	type	flux_calue	flux_std	lb	ub
        reaction	v1	fixed	100	1	0.001	100
        reaction	v2	free	100	1	0.001	100
        reaction	v3	free	50	1	0.001	100
        metabolites	AcCoA	fixed	1	1	0.001	100
        metabolites	OAC	fixed	1	1	0.001	100
        metabolites	OACs	fixed	1	1	0.001	100
        metabolites	OACx	fixed	1	1	0.001	100
        reversible	FUM	free	1	1	0.001	100

        format:
            'csv' : CSV.
            'text' : tab-deliminated text.

        Reterns
        ----------
        Dict: Dictionary of states dict.

        Examples
        --------
        >>> states = model.load_states("test.csv", format = 'csv')

        See Also
        --------
        """
        #
        # preparation of data
        #
        flux_dict = {}
        for i, id in enumerate(self.reaction_ids):
            flux_dict[id] = {
                'value':0.0,
                'stdev': 1.0,
                'type':"free",
                'reversible':self.reactions[id]['reversible'],
                'order':self.reactions[id]['order'],
                'lb':self.reactions[id]['lb'],
                'ub':self.reactions[id]['ub'],
            }
        conc_dict = {}
        for i, id in enumerate(self.metabolites):
            conc_dict[id] = {
                'value':0.0,
                'stdev': 1.0,
                'type':"free",
                'lb':self.metabolites[id]['lb'],
                'ub':self.metabolites[id]['ub'],
            }
        reversible_dict = {}
        for i, id in enumerate(self.reversible):
            reversible_dict[id] = {
                'value':0.0,
                'stdev': 1.0,
                'type':"free",
                'order':self.reversible[id]['order'],
                'lb':self.reversible[id]['lb'],
                'ub':self.reversible[id]['ub'],
            }

        with open(filename, 'r') as f:
            import csv
            if format == "text":
                reader = csv.reader(f, delimiter='\t')
            elif format == "csv":
                reader = csv.reader(f, dialect='excel')
            else:
                print("Unknown format!")
                f.close()
                return False

            for i, row in enumerate(reader):
                if i == 0: continue
                if len(row) < 7:
                    continue

                state, rid, type, value, stdev, lb, ub, *over = row
                if state == "reaction":
                    dict = flux_dict
                elif state == "metabolite":
                    dict = conc_dict
                elif state == "reversible":
                    dict = reversible_dict
                else:
                    continue
                if rid in dict:
                    dict[rid]['value'] = float(value)
                    dict[rid]['stdev'] = float(stdev)
                    dict[rid]['type'] = type
                    dict[rid]['lb'] = float(lb)
                    dict[rid]['ub'] = float(ub)

        f.close()
        return {"reaction":flux_dict, "metabolite":conc_dict, "reversible": reversible_dict}


    def save_states(self, dict, filename, format = 'text'):
        """
        Save state dict to a text/csv file with 'type' information.

        Parameters
        ----------
        model : a instance of mfapy.model for templete
        filename : filename of MDV data with the format.

        format:
            'csv' : CSV.
            'text' : tab-deliminated text.

        Reterns
        ----------
        Boolean: True/False

        Examples
        --------
        >>> model.save_states(model,  "test.csv", flux_dict, format = 'csv')

        See Also
        --------
        """

        #
        # preparation of data
        #
        Data = []
        Data.append(['State','Id','type','value','stdev','lb','ub'])
        for rid in self.reaction_ids:
            dict_temp = dict["reaction"]
            Data.append(['reaction', rid, str(dict_temp[rid]['type']),
                        str(dict_temp[rid]['value']),
                        str(dict_temp[rid]['stdev']),
                        str(dict_temp[rid]['lb']),
                        str(dict_temp[rid]['ub']),
                        ])
        for rid in dict["metabolite"]:
            dict_temp = dict["metabolite"]
            Data.append(['metabolite', rid, str(dict_temp[rid]['type']),
                        str(dict_temp[rid]['value']),
                        str(dict_temp[rid]['stdev']),
                        str(dict_temp[rid]['lb']),
                        str(dict_temp[rid]['ub']),
                        ])
        for rid in self.reversible.keys():
            dict_temp = dict["reversible"]
            Data.append(['reversible', rid, str(dict_temp[rid]['type']),
                        str(dict_temp[rid]['value']),
                        str(dict_temp[rid]['stdev']),
                        str(dict_temp[rid]['lb']),
                        str(dict_temp[rid]['ub']),
                        ])


        try:
            with open(filename, 'w', newline='') as f:
                import csv
                if format == "text":
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerows(Data)
                if format == "csv":
                    writer = csv.writer(f, dialect='excel')
                    writer.writerows(Data)
            return True
        except:
            return False


    def load_mdv_data(self, filename, format = 'text',output = "normal"):
        """
        Load MDV data from the text file. This function generate new instance of mfapy.mdv class from a instance of mfapy.model.

        Parameters
        ----------
        model : a instance of mfapy.model for templete
        filename : filename of MDV data with following format.
        Name	Spectrum	Select	MDV	Std
        Ala57	m0	1	0.566990778	0.000774686
        Ala57	m1	1	0.148623963	0.000774686
        Ala57	m2	1	0.039467636	0.000774686
        Ala57	m3	1	0.244917622	0.000774686

        format : "text" (defalut) or "csv"
        output : "normal" (defalut) or "debug"
        Reterns
        ----------
        Boolean: True/False

        Examples
        --------
        >>> mdv = load_mdv_data('filename')

        See Also
        --------

        """
        #
        # mdv
        #
        mdvloaded = mdv.MdvData(self.target_fragments)
        #
        mdvloaded.load(filename, format, output)
        #
        #
        return(mdvloaded)

