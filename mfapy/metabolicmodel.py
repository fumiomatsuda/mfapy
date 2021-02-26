#!/usr/bin/python
# -*- coding: utf-8 -*-
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

"""metabolicmodel.py:MetabolicModel class in mfapy

This is a core module of mfapy.

The module includes::

    MetabolicModel class

Todo:
    * Rewriting model construction part.
    * Check callbacklevel.

"""

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
    """Class of MetabolicModel

    An instances of this class has fuctions to analyze mass spectrometry based 13C MFA and INST MFA data.

    An metabolicModel instance generates MDVData and CarbonSouece classes

    Returns:
        instance: instance of metabolic model


    Examples:
        >>> reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_1_toymodel_model.txt", format = "text")
        >>> model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)

    Attributes (incomplete)::

        self.reactions = reactions
        self.metabolites = metabolites
        self.reversible = reversible
        self.target_fragments = target_fragments
        self.symmetry = {}
        self.carbon_source = {}
        self.experiments = {}

    """

    def __init__(self, reactions, reversible, metabolites, target_fragments, mode = "normal"):
        '''Generator of new instance
        Args:
            reactions (dict): Dictionary describing metabolite reactions

            reversible (dict): Dictionary for defining reversible reactions

            metabolites (dict): Dictionary including metabolite information

            target_fragments (dict): Dictionary of target_fragments

            mode (str): "normal" or "debug"

        '''
        #
        # Keep original information
        #
        self.reactions = copy.deepcopy(reactions)
        self.metabolites = copy.deepcopy(metabolites)
        self.reversible = copy.deepcopy(reversible)
        self.target_fragments = copy.deepcopy(target_fragments)

        number_of_errors = self.modelcheck()
        if number_of_errors > 0:
            print("There are critical problem(s) in the metabolic model")
            return

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
    def modelcheck(self):
        """Method to check medel definition data.

        Args:
            nothing

        Returns:
            Nothing.

        Examples:
            >>> model.datacheck()

        See Also:
            reconstruct

        History:
            Newly developed for version 059
        """
        import re
        number_of_errors = 0
        #
        # Check metabolite names
        # 6PG => "m6PG"
        # CO2_in => CO2in
        #
        metabolitename_changed = {}
        metabolitename_in_model = {}
        metabolitename_in_metabolites = {}

        for name in sorted(self.metabolites.keys(), key=lambda x: self.metabolites[x]['order']):
            oldname = str(name)
            name = name.replace("_", "")
            if not re.match(r"^[A-Za-z]",name):
                name = "m"+name
            if not name == oldname:
                print("Caution:", oldname, " was changed to", name)
                self.metabolites[name] = self.metabolites.pop(oldname)

            metabolitename_changed[oldname] = name
            metabolitename_changed[name] = name
            metabolitename_in_metabolites[oldname] = name



        #
        # Compare metabolite names in reactions and stoichiometry with Metabolite list
        #
        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #
            reaction = self.reactions[id]['stoichiometry']
            if re.match(r".+-->.+",self.reactions[id]['reaction']):
                reaction = reaction + "+" + self.reactions[id]['reaction']
            reaction = re.sub(r'\{.+\}', "", reaction)
            temps = re.split(r"-->|\+", reaction)
            for metabolite in temps:
                metabolitename_in_model[metabolite] = 1
                if not metabolite in metabolitename_changed:
                    number_of_errors = number_of_errors + 1
                    print ("Error!:",id, metabolite, "is not listed in //Metabolites in the metabolic model")
        #
        # Compare metabolite names in target_fragment with Metabolite list
        #
        for id in sorted(self.target_fragments.keys(), key=lambda x: self.target_fragments[x]['order']):
            #
            # get reaction formula and substrate, product names
            #

            atommap = self.target_fragments[id]['atommap']
            if not re.match(r".+_.+",atommap):
                continue
            atommap.replace(" ","")
            fragments = atommap.split("+")
            fragments_array = []
            #
            # For each fragments
            #
            for fragment in fragments:
                temps = fragment.split("_")
                metabolite = "_".join(temps[:-1])
                if not metabolite in metabolitename_changed:
                    number_of_errors = number_of_errors + 1
                    print ("Error!:",id, metabolite, "is not listed in //Metabolites in the metabolic model")
                if not metabolite in metabolitename_in_model:
                    number_of_errors = number_of_errors + 1
                    print ("Error!:",id, metabolite, "does not exist in metabolid network")
        #
        # Chenk metabolites not used in the metablic network
        #
        for metabolitename in metabolitename_in_metabolites:
            if not metabolitename in metabolitename_in_model:
                 print ("Caution:",metabolitename, "was not used in the metabolid network")

            if number_of_errors > 0:
                return(number_of_errors)
        #
        # Check metabolite name in stoichiometry
        #
        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #
            reaction = self.reactions[id]['stoichiometry']
            reaction.replace(" ","")
            substrate, product = reaction.split("-->")
            #
            # Substrate
            #
            new_text = ""
            newarray = []
            for metabolite in substrate.split("+"):
                #
                # Check coefficient: {0.5}Pyr
                #
                if len(metabolite.split("}")) == 1:
                    # no coefficient
                    metabolite_name = metabolite
                    newarray.append(metabolitename_changed[metabolite_name])
                elif len(metabolite.split("}")) == 2:
                    # coefficient data
                    stnum, metabolite_name = metabolite.split("}")
                    newarray.append(stnum+"}"+metabolitename_changed[metabolite_name])
            new_text = "+".join(newarray)
            new_text = new_text+"-->"

            #
            # Product
            #
            newarray = []
            for metabolite in product.split("+"):
                if len(metabolite.split("}")) == 1:
                    metabolite_name = metabolite
                    newarray.append(metabolitename_changed[metabolite_name])
                elif len(metabolite.split("}")) == 2:
                    stnum, metabolite_name = metabolite.split("}")
                    newarray.append(stnum+"}"+metabolitename_changed[metabolite_name])
            new_text = new_text+"+".join(newarray)
            if not reaction == new_text:
                print("Caution!:",id, 'stoichiometry', reaction, "was changed to" , new_text)
                self.reactions[id]['stoichiometry'] = new_text

        #
        # Check metabolite name in reaction
        #
        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #

            reaction = self.reactions[id]['reaction']
            if not re.match(r".+-->.+",reaction):
                continue
            reaction.replace(" ","")
            substrate, product = reaction.split("-->")
            #
            # Substrate
            #
            new_text = ""
            newarray = []
            for metabolite in substrate.split("+"):
                newarray.append(metabolitename_changed[metabolite])
            new_text = "+".join(newarray)
            new_text = new_text+"-->"
            #
            # Product
            #
            newarray = []
            for metabolite in product.split("+"):
                newarray.append(metabolitename_changed[metabolite])
            new_text = new_text+"+".join(newarray)
            if not reaction == new_text:
                print("Caution!:",id, 'reaction', reaction, "was changed to" , new_text)
                self.reactions[id]['reaction'] = new_text

        #
        # Check metabolite name in target_fragment data
        #
        for id in sorted(self.target_fragments.keys(), key=lambda x: self.target_fragments[x]['order']):
            #
            # get reaction formula and substrate, product names
            #

            atommap = self.target_fragments[id]['atommap']
            if not re.match(r".+_.+",atommap):
                continue
            atommap.replace(" ","")
            fragments = atommap.split("+")
            fragments_array = []
            #
            # For each fragments
            #
            for fragment in fragments:
                temps = fragment.split("_")
                metabolite = "_".join(temps[:-1])
                atomdata = temps[-1]
                #print('atommap', id, atommap, fragment, temps, metabolite, atomdata)
                fragments_array.append(metabolitename_changed[metabolite]+"_"+atomdata)
            new_text = "+".join(fragments_array)
            if not atommap == new_text:
                print("Caution!:",id, 'target_fragments', atommap, "was changed to" , new_text)
                self.target_fragments[id]['atommap'] = new_text
        #
        # Check number of carbons in atom mappinng data
        #
        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #
            reaction = self.reactions[id]['reaction']
            if not re.match(r".+-->.+",reaction):
                continue
            meatbolites = re.split(r'-->|\+', reaction)
            atoms = re.split(r'-->|\+', self.reactions[id]['atommap'])
            if not len(meatbolites) == len(atoms):
                print("There is an problem in ", id,  reaction, self.reactions[id]['atommap'])
            for metabolite, atom in zip(meatbolites, atoms):
                if not len(atom) == self.metabolites[metabolite]['C_number']:
                    print("Error!", id, "Carbon number of ", metabolite, "does not match the metabolite's carbon number:", self.metabolites[metabolite]['C_number'])
                    number_of_errors = number_of_errors + 1
        #
        # Check Carbon number in target_fragment data
        #
        for id in sorted(self.target_fragments.keys(), key=lambda x: self.target_fragments[x]['order']):
            #
            # get reaction formula and substrate, product names
            #

            atommap = self.target_fragments[id]['atommap']
            if not re.match(r".+_.+",atommap):
                continue
            atommap.replace(" ","")
            fragments = atommap.split("+")
            fragments_array = []
            #
            # For each fragments
            #
            for fragment in fragments:
                temps = fragment.split("_")
                metabolite = "_".join(temps[:-1])
                atomdata = temps[-1]
                carbonnumber = self.metabolites[metabolite]['C_number']
                for number in atomdata.split(":"):
                    if int(number) > carbonnumber:
                        print("Error!", id, "Carbon number in ", fragment, "does not match the metabolite's carbon number:", carbonnumber)
                        number_of_errors = number_of_errors + 1
        #
        # Check Metabolic Network
        #
        metabolitename_in_model = {}
        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #
            reaction = self.reactions[id]['stoichiometry']
            reaction.replace(" ","")
            substrate, product = reaction.split("-->")
            #
            # Substrate
            #
            new_text = ""
            newarray = []
            for metabolite in substrate.split("+"):
                #
                # Check coefficient: {0.5}Pyr
                #
                if len(metabolite.split("}")) == 1:
                    # no coefficient
                    metabolite_name = metabolite
                    if not metabolite_name in metabolitename_in_model:
                        metabolitename_in_model[metabolite_name] = {'p':[], 's':[]}
                    metabolitename_in_model[metabolite_name]["s"].append(id)
                elif len(metabolite.split("}")) == 2:
                    # coefficient data
                    stnum, metabolite_name = metabolite.split("}")
                    if not metabolite_name in metabolitename_in_model:
                        metabolitename_in_model[metabolite_name] = {'p':[], 's':[]}
                    metabolitename_in_model[metabolite_name]["s"].append(id)
            #
            # Product
            #
            newarray = []
            for metabolite in product.split("+"):
                if len(metabolite.split("}")) == 1:
                    metabolite_name = metabolite
                    if not metabolite_name in metabolitename_in_model:
                        metabolitename_in_model[metabolite_name] = {"p":[], "s":[]}
                    metabolitename_in_model[metabolite_name]["p"].append(id)
                elif len(metabolite.split("}")) == 2:
                    stnum, metabolite_name = metabolite.split("}")
                    if not metabolite_name in metabolitename_in_model:
                        metabolitename_in_model[metabolite_name] = {"p":[], "s":[]}
                    metabolitename_in_model[metabolite_name]["p"].append(id)
        for metabolite in metabolitename_in_model:
            substrate_ids = metabolitename_in_model[metabolite]["s"]
            product_ids = metabolitename_in_model[metabolite]["p"]
            if self.metabolites[metabolite]['excreted'] == 'excreted':
                if len(substrate_ids) > 0:
                    print("Error!: Excreted metabolite", metabolite, "is a substrate in ", substrate_ids)
                    number_of_errors = number_of_errors + 1
                continue
            if self.metabolites[metabolite]['carbonsource'] == 'carbonsource':
                if len(product_ids) > 0:
                    print("Error!: Carbon source metabolite", metabolite, "is a product in ", product_ids)
                    number_of_errors = number_of_errors + 1
                continue
            if len(product_ids) == 0:
                print("Error!: Metabolite", metabolite, "is not produced in the metabolic network, but used in reaction(s) ", substrate_ids)
                number_of_errors = number_of_errors + 1
            if len(substrate_ids) == 0:
                print("Error!: Metabolite", metabolite, "is not used in the metabolic network, but produced  in reaction(s) ", product_ids)
                number_of_errors = number_of_errors + 1

        #
        # Check inconsistency in atom mapping
        #

        for id in sorted(self.reactions.keys(), key=lambda x: self.reactions[x]['order']):
            #
            # get reaction formula and substrate, product names
            #
            atommap = self.reactions[id]['atommap']
            if not re.match(r".+-->.+",atommap):
                continue
            substrate, product = atommap.split("-->")
            substrate_dict = {}
            for metabolite in substrate.split("+"):
                for letter in metabolite:
                    if not letter in substrate_dict:
                        substrate_dict[letter] = 0
                    substrate_dict[letter] = substrate_dict[letter] + 1
            product_dict = {}
            for metabolite in product.split("+"):
                for letter in metabolite:
                    if not letter in product_dict:
                        product_dict[letter] = 0
                    product_dict[letter] = product_dict[letter] + 1
            for letter in substrate_dict:
                if substrate_dict[letter] > 1:
                    print("Error!:", id,  "Symbol", letter, "is not unique in atom mapping data of substrate ", substrate)
                    number_of_errors = number_of_errors + 1
            for letter in product_dict:
                if product_dict[letter] > 1:
                    print("Error!:", id,  "Symbol", letter, "is not unique in atom mapping data of product ", product)
                    number_of_errors = number_of_errors + 1
                if not letter in substrate_dict:
                    print("Error!:", id,  "Symbol", letter, "in product atom mapping", product, "does not exist in the subtrate atom mapping", substrate)
                    number_of_errors = number_of_errors + 1

        return(number_of_errors)


    def update(self):
        """Method to generate stoichiometry matrix.

        A metabolic model have to be updated when any types (fixed, free, fitting, pseudo) of reactions, metabolites, and reversible reactions are modified.

        Args:
            Not required.

        Examples:
            >>> model.update()

        See Also:
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
            #if self.reactions[id]['use'] != 'use':
            #    continue
            #
            # Order of the reaction in self.reaction_ids is stored
            #
            self.reactions[id]["position_in_tmp_r"] = counter
            #
            # stoichiometry substrate, key => coefficient
            #
            self.reactions[id]["stoichiometry_metabolite_list"] = {}
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
                self.reactions[id]["stoichiometry_metabolite_list"][metabolite_name] = float(stnum) * -1.0
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
                self.reactions[id]["stoichiometry_metabolite_list"][metabolite_name] = float(stnum) * 1.0
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
        """Method to reconstruct the metabolic model to generate new stoichiometry matrix and a function to
        calculate MDV data.

        Args:
            mode (str): experimental

        Returns:
            calmdv_text (dict) The dictionary has two keys. "calmdv" and "diffmdv" are instance of functions for 13C-MFA and INST-13C-MFA to generate MDV from a given metabolic flux, respectively.

        Examples:
            >>> model.reconstruct()


        See Also:
            update()

        """

        self.update()
        calmdv_text, ccalmdv_text = self.generate_calmdv(mode)
        self.calmdv_text = calmdv_text
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
        """Constructer of the python function to calculate MDV of target fragments from a metabolic flux distribution.

        This function is performed in reconstruct()

        Args:
            mode (str): optimize EMU network reduction (experimental)

        Returns:
            string, cython_string (str) Texts of python code describing calmdv and diffmdv functions for python 3 and Cython (Experimental)

        Examples:
            >>> calmdv_text, ccalmdv_text = self.generate_calmdv(mode)

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
                            if self.configuration['callbacklevel'] == 7:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)

                        else:
                            for reaction_id in matrix[product].get(substrate):
                                if '*' in reaction_id:
                                    cell.append("0.5 * " + reaction_id.replace("*", ""))
                                else:
                                    cell.append(reaction_id)
                            cell_string = ' + '.join(cell)
                            string += "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string + "\n"
                            if self.configuration['callbacklevel'] == 7:
                                print(product,'\t',substrate, "\tA_" + str(emu_size) + "["+str(emu_intermediate_of_this_layer.index(product) * size_i + emu_intermediate_of_this_layer.index(substrate))+"] = " + cell_string)

            string += "\tA_" + str(emu_size) + ".resize(("+str(size_i)+","+ str(size_i)+"))\n"

        # Output B
            string += "\tB_" + str(emu_size) + " =  numpy.zeros(("+str(size_i * size_s)+",))\n"
            for product in emu_intermediate_of_this_layer:
                for substrate in emu_sourse_of_this_layer:
                    cell = []
                    if self.configuration['callbacklevel'] == 7:
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
                        if self.configuration['callbacklevel'] == 7:
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
                        if self.configuration['callbacklevel'] == 7:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]))
                        string += "\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"

                else:
                    compound = compound_of_EMU(substrate)
                    size = size_of_EMU(substrate);
                    for numberofisotope in range(size+1):
                        if self.configuration['callbacklevel'] == 7:
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
                pass
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

                            if self.configuration['callbacklevel'] == 7:
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

                            if self.configuration['callbacklevel'] == 7:
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

                            if self.configuration['callbacklevel'] == 7:
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
                        if self.configuration['callbacklevel'] == 7:
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
                        if self.configuration['callbacklevel'] == 7:
                            print(substrate ,'\t','',"\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]))
                        string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"
                        cython_string += "\t\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation[numberofisotope]) + "\n"
                        string += "\t\tstring += '\\tY_" + str(emu_size) + "["+str(emu_sourse_of_this_layer.index(substrate)*(emu_size + 1) + numberofisotope)+"] = " + '+'.join(equation_forfunc[numberofisotope]) + "\\n'\n"

                else:
                    compound = compound_of_EMU(substrate)
                    size = size_of_EMU(substrate);
                    for numberofisotope in range(size+1):
                        if self.configuration['callbacklevel'] == 7:
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

            # 200517 ver 056 modified
            H0ratio = 0.9893
            H1ratio = 0.0107

            (stem, pos) = emu.split("_")
            pos = pos.replace(':','') #200517 ver 056 modified 炭素移動の表記が変わったので
            number_of_carbons = len(pos)
            #string += "\ty0["+str(position)+"] = 1.0\n"
            string += "\ty0["+str(position)+"] = "
            string += str(H0ratio**number_of_carbons)
            string += "#"+str(emu)+" "+str(position)
            string += "\n"
            string += "\ty0["+str(position+1)+"] = "
            string += str((H0ratio**(number_of_carbons-1.0))*H1ratio*number_of_carbons)
            string += "\n"
            if number_of_carbons > 2:
                string += "\ty0["+str(position+2)+"] = "
                string += str(1.0- H0ratio**number_of_carbons - (H0ratio**(number_of_carbons-1.0))*H1ratio*number_of_carbons)
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
        if self.configuration['callbacklevel'] == 7:
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
                pass
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
        """Setter of a configuration of model parameters.

        Args:
            *args (list): "parameter_name = number"

        Returns:
            Nothing.

        Parameters::

            'callbacklevel': default = 0,
                'callbacklevel' determines a frequency level of callbacks from metabolic model.
                # 0: No call back
                # 1: Error report
                # 2: Simple pregress report for grid search
                # 3: Detailed pregress report for grid search
                # 4: Simple pregress report for flux fitting
                # 5: Detailed pregress report for flux fitting
                # 7: Detailed report for model constrution
                # 8: Debug
            'iteration_max': default = 1000,
                Maximal number of interations (steps) in each fitting task.
            'default_reaction_lb': default = 0.001,
            'default_reaction_ub': default = 5000.0,
            'default_metabolite_lb': default = 0.001,
            'default_metabolite_ub': default = 500.0,
            'default_reversible_lb': default = -300,
            'default_reversible_ub': default = 300.0,
                Default lower and upper boundaries used only when there is no data in the model definition file.
            'grid_search_iterations': default = 1,
                Number of trials for model fitting at each point in grid search.
            'initial_search_repeats_in_grid_search': default = 50,
                Number of initial metabolic flux disributions generated for each trial in grid search.
                Best initial metabolic flux disribution with smallest RSS is employed.
            'initial_search_iteration_max': default = 1000,
                Maximal number of interations (steps) allowed in each task to find feasible initial metabolic flux distribution.
            'add_naturalisotope_in_calmdv':default = "no"
                For INST-13C MFA. Natural isotope is added to MDVs of all intermediate at time = 0.
            'number_of_repeat':default = 3,
                For "deep" function of self.fitting flux. "Global" => "Local" optimization methods are repeated for n times.
            'ncpus':default = 3,
                Number of CPUs used in parallel processing
            'odesolver': default = "scipy" # or "sundials"
                Only "scipy" is avaliable due to memory leak problem of sundials

        Examples:

            >>> model.set_configuration(callbacklevel = 1)
            >>> model.set_configuration(iteration_max = 1000)


        See Also
        --------

        """

        for word in kwargs:
            #if word in self.configuration:
            self.configuration[word] = kwargs[word]
            if self.configuration['callbacklevel'] >= 4:
                print(word,'is set at',kwargs[word])

    def set_constrain(self, group, id, type, value = 0.1, stdev = 1.0):
        """Setter of a metabolic constrains.

        Please perform update() after model modification.

        Args:
            group (str):group of parameters ("reaction", "reversible", "metabolite")

            id (str): id in model definition

            type (str): type of constrain "fitting", "fixed", "free", "pseudo":
                * "fitting" : For the calclation of RSS by using "value" and stdev
                * "fixed"   : Fixed to "value". Not used for the calculation of RSS.
                * "free"    : Free between lower and upper boundary. Not used for the calculation of RSS.
                * "pseudo"  : Pseudo reaction is a special "free" reaction to consider G-index. The reaction is ignored in a stoichiometry matrix. However, the reaction is considered in the MDV calculation.

            value (float): Level of metabolic flux

            stdev (float): Standard deviation of metabolic flux level used for 'fitting' type reactions.

        Returns:
            Booleans: True/False.


        Examples:
            >>> model.set_fixed_reaction("reaction", 'pgi', "fixed", 100.0,1)


        See Also:
            get_constrain

        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 1:
                print('ERROR:', group,' is wrong group')
                return False

        if type not in ["fitting", "fixed", "free", "pseudo"]:
            if self.configuration['callbacklevel'] >= 1:
                print('ERROR:', type, ' is wrong type')
                return False

        if id in dic_temp:
            dic_temp[id]['type'] = type
            dic_temp[id]['value'] = float(value)
            if float(stdev) >= 0:
                dic_temp[id]['stdev'] = float(stdev)
            if self.configuration['callbacklevel'] >= 5:
                print(id,' is set as ',type, group, ' with value at ',value, ' and stdev at ', stdev, '.')
            return True
        if self.configuration['callbacklevel'] >= 1:
            print('ERROR:', id,' not existed in the model')
        return False

    def get_constrain(self, group, id):
        """Getter of constrains of metabolic reaction.


        Args:
            group (str) :group of parameters ("reaction", "reversible", "metabolite")

            id (str) : id in model definition


        Returns:
            * type (str) type of constrain "fitting", "fixed", "free", "pseudo"
            * value (float) Level of metabolic flux
            * stdev (float) Standard deviation of metabolic flux level used for 'fitting' type reactions.

        Examples:
            >>> type, value, stdev = model.get_constrain("reaction","pgi")

        See Also:
            set_constrain

        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 1:
                print('ERROR:', group,' is wrong group')
                return False, False, False

        if id in dic_temp:
            types = dic_temp[id]['type']
            value = float(dic_temp[id]['value'])
            stdev = float(dic_temp[id]['stdev'])
            return types, value, stdev
        if self.configuration['callbacklevel'] >= 1:
            print('ERROR:', id,' not existed in the model')
        return False, False, False

    def set_boundary(self, group, id, lb, ub):
        """ Setter of a lower and upper boundaries of metablic states

        Args:
            group (str): metabolite, reaction, reversible

            id (str): id of state parameter

            lb (float): lower boundary

            ub (float): upper boundary

        Returns:
            Booleans: True/False.


        Examples:
            >>> model.set_boundary("reaction", 'HEX1', 0.001, 100.0)

        See Also:
            set_constrain
        """
        if group == "reaction":
            dic_temp = self.reactions
        elif group == "metabolite":
            dic_temp = self.metabolites
        elif group == "reversible":
            dic_temp = self.reversible
        else:
            if self.configuration['callbacklevel'] >= 1:
                print('ERROR:', group,' is wrong group')
                return False


        if id in dic_temp:
            dic_temp[id]['lb'] = float(lb)
            dic_temp[id]['ub'] = float(ub)
            if self.configuration['callbacklevel'] >= 5:
                print("lower and upper boundaries of", id,' are set to be ', lb, "and", ub)
            return True
        if self.configuration['callbacklevel'] >= 1:
            print('ERROR:', id,' not existed in the model')
        return False


    def set_constraints_from_state_dict(self, dict):
        """Reaction types, value, stdev, lb, and ub are set from the state dict data at once

        Args:
            dict (dict): Dictionary of flux. Data in 'value', 'stdev', 'lb', 'ub' and 'type' fields are used.

        Examples:
            >>> model.set_constraints_from_state_dict(flux_dict)



        """
        #
        # preparation of data
        #
        counter = 0
        for i, (group, id) in enumerate(self.vector["ids"]):
            if id in dict[group]:
                type = dict[group][id]['type']
                value = dict[group][id]['value']
                stdev = dict[group][id]['stdev']
                lb = dict[group][id]['lb']
                ub = dict[group][id]['ub']
                self.set_constrain(group, id, type, value = value, stdev = stdev)
                self.set_boundary(group, id, lb, ub)
                counter = counter   + 1
                if self.configuration['callbacklevel'] >= 5:
                    print(id,' is set as ',type, group, ' with value at ',value, ' and stdev at ', stdev, " and boundaries", lb, ub,'.')
        self.update()
        if self.configuration['callbacklevel'] >= 3:
            print("Batch setting of", counter, "constrains.")
        return True


    def generate_state_dict(self, tmp_r):
        """Generator of a state dict from a given vector.

        Args:
            tmp_r (array): tmp_r = numpy.dot(matrixinv, Rm_intial)

        Returns:
            dict : Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

        Examples:
            >>> state_dict = model.generate_state_dict(tmp_r)


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
        return self.generate_carbon_source_template()

    def generate_carbon_source_template(self):
        """Generator of a templete of CarbonSourse instance. Labeling information of carbon sources will be added to the instance.

        Carbon source compounds are derived from the metabolic model (//Metabolites)

        Args:
            Not required.

        Examples:
            >>> carbon_source_idv = model.generate_carbon_source_templete()

        See Also:
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
        """Generator of a MdvData instance including MDV data generated from given flux and carbon sources.

        Args:
            flux (dict): Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

            carbon_sources (instance): Instance of CarbonSource

            timepoint (array): For INST-13C MFA. An array of time points of measurement in MDV

            startidv (array): array of idv as starting isotope distribution for INST

        Returns:
            instance* MdvData instance

        Examples:
            >>> mdv = model.generate_mdv(flux, mdv_carbon_sources)

        See Also:
            generate_carbonsource_MDV

        History:
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
        """Setter of an 'experiment' to metabolic model.

        Here an 'experiment' indicated a set of a carbon source labeling pattern and a measured MDV data.

        Parallel labeling experiment can be performed by setting multiple sets of 'experiment'.


        Args:
            name (str): name of the experiment (unique)

            mdv (instance): MDVdata instance including measured MDV data of target fragments.

            carbon_sources (instance): Instance of carbon source object

            startidv (array): array of idv as starting isotope distribution for INST


        Examples:
            >>> model.set_experiment('ex1', mdv1, cs1)
            >>> model.set_experiment('ex2', mdv2, cs2, startidv)


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
                        # 200517 ver 056 modified
                        H0ratio = 0.9893
                        H1ratio = 0.0107

                        (stem, pos) = emu[0].split("_")
                        pos = pos.replace(':','') #200517 ver 056 modified 炭素移動の表記が変わったので
                        number_of_carbons = len(pos)
                        self.experiments[name]["y0"][position] = H0ratio**number_of_carbons
                        self.experiments[name]["y0"][position + 1] = (H0ratio**(number_of_carbons-1.0))*H1ratio*number_of_carbons
                        if number_of_carbons > 2:
                            self.experiments[name]["y0"][position + 2] = 1.0- H0ratio**number_of_carbons - (H0ratio**(number_of_carbons-1.0))*H1ratio*number_of_carbons
        if self.configuration['callbacklevel'] >= 3:
            print("Set experiment: ", name,' was added to the metabolic model.')
        return True

    def calc_idv(self, flux, carbon_sources):
        """Calc IDV from given flux and carbon_sources.

        Isotopomer distribution vector (IDV) was used to calc MDV at time = 0 in the INST mode.

        Examples:
            >>> idv = model.calc_idv(flux, carbon_sources)

        Args:

            flux (dict): Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

            carbon_sources (instance): Instance of carbon source object

        Returns:
            array : array of idv


        See Also:
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
        """Clear all 'experiment set(s)' in the metabolic model.


        Args:
            Not required

        Examples:
            >>> model.clear_experiment()


        See Also:

            model.set_experiment()

        History:
            Newly developed at 4/9/2014


        """
        names = list(self.experiments.keys())
        self.experiments = {}
        if self.configuration['callbacklevel'] >= 3:
            print("Clear experiment: ", names,' are removed from the metabolic model.')

        return True

    def generate_state(self, template = []):
        """Generator of a random metabolic flux distribution.

        Used parameters in self.configuration().

            * 'initial_search_iteration_max': Maximal number of interation in the optimizers. Example: model.set_configuration(initial_search_iteration_max = 1000)

        Args:
            template (dict): Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels. When template is available, metabolic state most similar to the template is generated. The function is used in the grid search.

        Returns:
            * flux (dict) Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

            * independent_flux (array) List of independent flux

        Examples:
            >>> flux, independent_flux = model.generate_intial_flux_distribution()

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
        """Generator of random initial metabolic states for model fitting

        Initial metabolic states are randomly generated for "iterations" times from which better states with lower RSS (number_of_initial_states) were selected.


        Used parameters in self.configuration().

            * 'initial_search_iteration_max': Maximal number of interation in the optimizers. Example: model.set_configuration(initial_search_iteration_max = 1000)

            * 'ncpus': Number of CPUs used in parallel processing. Example: model.set_configuration(ncpus = 16)

        Args:
            method (str): computational method
                * 'normal' :  Generation of flux distributions by single thread.
                * 'parallel' : Parallel generation of flux distributions using parallel processing.

            iterations (int): Number of trials to generate initial metabolic state. A trial sometimes failed due to a bad optimization state.

            initial_states (int): Number of initial metabolic state(s) with lower RSS to be generated.

            template (array): Array of templete metabolic state(s) dictionaries. When template is available, initial metabolic state(s) similiar to the templete metabolic state are generated. This function is used in the grid search fucntion.

        Returns:
            * flux (array or dict) : Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

                * [number_of_initial_states > 1 ] Array of dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

                * [number_of_initial_states == 1 ] Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

            * state  (str or array) : stop condition
                * [number_of_initial_states > 1 ] Str of stop condition.

                * [number_of_initial_states == 1 ] Array of Str of stop condition.

        Examples:
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
            if self.configuration['callbacklevel'] >= 1:
                print('No experiment was set to the modelnames.')
            return False
        #
        # Mas=number of MKL thread control
        #
        try:
            import mkl
            mkl.set_num_threads(1)
        except:
            if self.configuration['callbacklevel'] >= 1:
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

        if self.configuration['callbacklevel'] >= 3:
            print("Generation of initial state(s) was started by",method,"method.")

        if method == "parallelpp":
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

            if (self.configuration['callbacklevel'] >= 2):
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
        #
        # joblib
        #
        elif method == "parallel":
            #Set ncpus
            if 'ncpus' in self.configuration:
                ncpus = self.configuration['ncpus']
            else:
                ncpus = 1
            #
            # tuple of all parallel python servers to connect with
            #
            try:
                from joblib import Parallel, delayed
            except:
                if self.configuration['callbacklevel'] >= 1:
                    print("This function requires joblib!")
                return False
            jobs = []
            for i in range(iterations):
                parameters = (numbers, vectors, matrixinv, template ,initial_search_iteration_max)
                jobs.append(parameters)
            result = Parallel(n_jobs=ncpus)([delayed(optimize.initializing_Rm_fitting)(numbers, vectors, matrixinv, template ,initial_search_iteration_max) for (numbers, vectors, matrixinv, template ,initial_search_iteration_max) in jobs])


            for results in result:
                if results == None:
                    continue
                tmp_r, Rm_temp, Rm_ind, state = results
                if state == "Determined":
                    flux_temp_r = self.generate_state_dict(tmp_r)
                    fluxes.append(flux_temp_r)
                    rsses.append(self.calc_rss(flux_temp_r))

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


        fluxes = [fluxes[x] for x in order]
        rsses = [rsses[x] for x in order]
        state = int(len(fluxes))

        if self.configuration['callbacklevel'] == 4:
            print(state, "initial state(s) was generated.")
        if self.configuration['callbacklevel'] >= 5:
            print(state, "initial state(s) was generated. RSS of best one is: ", rsses[0])


        if len(fluxes) == 0:
            return state, fluxes
        else:
            fluxes = fluxes[:initial_states]
        if initial_states == 1:
            return state, fluxes[0]
        else:
            return state, fluxes


    def fitting_flux(self, method = 'SLSQP', flux = [], output = 'result'):
        """Method for fitting a metabolic model to multiple experiments

        Used parameters in self.configuration().

            * 'iteration_max': Maximal number of interation in the optimizers. Example: model.set_configuration(iteration_max = 1000)

            * 'number_of_repeat': Number of repeated execution by 'deep' functions. Example: model.set_configuration(number_of_repeat = 2)

            * 'ncpus': Number of CPUs used in parallel processing. Example: model.set_configuration(ncpus = 16)

        Args:
            method (str): Algorism used for optimization.

                * 'deep': Repeated execution of SLSQP & LN_PRAXIS

                * 'SLSQP': Sequential Least SQuares Programming algorithm implemented in scipy

                * 'COBYLA': Constrained Optimization By Linear Approximation implemented by scipy

                * "LN_BOBYQA":The BOBYQA algorithm for bound constrained optimization without derivatives by nlopt

                * "LN_SBPLX": Sbplx (based on Subplex, (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces) by nlopt

                * "LN_NEWUOA": The NEWUOA software for unconstrained optimization without derivatives by nlopt

                * "LN_COBYLA": Constrained Optimization By Linear Approximation implemented by nlopt

                * "LN_PRAXIS": Gradient-free local optimization via the "principal-axis method" implemented by nlopt

                * "LN_NELDERMEAD": the original Nelder-Mead simplex algorithm by nlopt

                * "GN_ESCH": Evolutionary Algorithm for global optimization by nlopt

                * "GN_CRS2_LM": Controlled random searchimplemented for global optimizatoin by nlopt

                * "GN_DIRECT_L":  DIviding RECTangles algorithm for global optimization by nlopt

                * "GN_IRES": Improved Stochastic Ranking Evolution Strategy by nlopt

                * "GN_AGS": AGS NLP solver by nlopt


            flux (dict): Initial flux distribution generated by self.generate_initial_states():

                * When flux is a dict of one state, single fitting trial is executed.

                * When flux is a array of multiple dists of states, muptiple fitting trial is exected by using the parallel processing.

            output (str): Output method.

                * 'result' (default): Fitting problem is solved inside of the function.

                * 'for_parallel': Fitting problems is NOT solved inside of the function. In this mode, a tupple of parameters for fit_r_mdv_pyopt() functions are generated for parallel processinng.


        Returns:
            * state  (str or array) : stop condition

                * [number_of_initial_states == 1 ] Str of stop condition.

                * [number_of_initial_states > 1 ] Array of str of stop condition.


            * flux_opt (array or dict): Initial flux distribution generated by self.generate_initial_states()

                * [number_of_initial_states == 1 ] dict of optimzed metabolic state.

                * [number_of_initial_states > 1 ] Array of dict of optimzed metabolic state. Sorted by ascending order of their RSS levels


            * RSS_bestfit (array or float) : Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

                * [number_of_initial_states == 1 ] Dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.

                * [number_of_initial_states > 1 ] Array of dictionary of metabolic state inclucing metabolic flux and metabolite pool size levels.


        Examples::

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
            # Generating parameters for parallel processing
            #
            >>> parameters = model.fitting_flux(method = "SLSQP", flux = flux_initial, output = "for_parallel")

        See Also:
            generate_initial_states()


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 1:
                print('No experiment was set to the modelnames.')
            return False
        if flux == []:
            if self.configuration['callbacklevel'] >= 1:
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
            elif method == "COBYLA":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "COBYLA")
            elif method in ["LN_COBYLA","LN_BOBYQA","LN_NEWUOA","LN_PRAXIS","LN_SBPLX","LN_NELDERMEAD","GN_DIRECT_L","GN_AGS","GN_CRS2_LM","GN_ESCH","GN_IRES"]:
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, method)
            elif method == "deep":
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux)
            else:
                parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux, "SLSQP")
            return parameters



        if isinstance(flux, dict):
            if self.configuration['callbacklevel'] >= 3:
                print("Trying to start fitting task using",method,"method in single mode")
            if method == "SLSQP":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "SLSQP")
            elif method == "COBYLA":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "COBYLA")
            elif method in ["LN_COBYLA","LN_BOBYQA","LN_NEWUOA","LN_PRAXIS","LN_SBPLX","LN_NELDERMEAD","GN_DIRECT_L","GN_AGS","GN_CRS2_LM","GN_ESCH","GN_IRES"]:
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_nlopt(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = method)
            elif method == "deep":
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_deep(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux)
            else:
                state, kai, opt_flux, Rm_ind_sol = optimize.fit_r_mdv_scipy(configuration, self.experiments, numbers, vectors, self.matrixinv, self.func, flux, method = "SLSQP")


            if self.configuration['callbacklevel'] == 3:
                print("Fitting task in single mode was finished.")
            if self.configuration['callbacklevel'] == 4:
                print("Fitting task in single mode was finished with finishing state",state,".")
            if self.configuration['callbacklevel'] >= 5:
                print("Fitting task in single mode was finished with finishing state",state," RSS is: ", kai)


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
                from joblib import Parallel, delayed
            except:
                print("This function requires joblib")
                return False

            jobs = []

            if self.configuration['callbacklevel'] >= 3:
                print("Trying to start fitting task using", method, "method using parallel processing")

            if method in ["SLSQP","COBYLA"]:
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, method)
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_scipy)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) in jobs])

            elif method in ["GN_CRS2_LM","LN_COBYLA","LN_BOBYQA","LN_NEWUOA","LN_PRAXIS","LN_SBPLX","LN_NELDERMEAD","GN_DIRECT_L","GN_AGS","GN_CRS2_LM","GN_ESCH","GN_IRES"]:
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, method)
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_nlopt)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) in jobs])

                    #jobs.append([i, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                    # (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                    # ("numpy","nlopt","scipy","scipy.integrate"))])
            elif method == "deep":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp)
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_deep)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) in jobs])

                    #jobs.append([i, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                    # (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                    # ("numpy","nlopt","scipy","scipy.integrate"))])
            else:
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp)
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_deep)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) in jobs])

                    #jobs.append([i, job_server.submit(optimize.fit_r_mdv_deep, parameters,
                    # (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                    # ("numpy","nlopt","scipy","scipy.integrate"))])
            """
            elif method == "LN_PRAXIS":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, "LN_PRAXIS")
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_nlopt)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) in jobs])

                    #jobs.append([i, job_server.submit(optimize.fit_r_mdv_nlopt, parameters,
                    # (optimize.calc_MDV_residue_scipy, optimize.calc_MDV_residue_nlopt,optimize.fit_r_mdv_scipy,optimize.fit_r_mdv_nlopt),
                    # ("numpy","nlopt","scipy","scipy.integrate"))])
            elif method == "GN_CRS2_LM":
                for i, flux_temp in enumerate(flux):
                    parameters = (configuration, self.experiments, numbers, vectors, self.matrixinv, self.calmdv_text, flux_temp, "GN_CRS2_LM")
                    jobs.append(parameters)
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_nlopt)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp, method) in jobs])
            """

            state_list = []
            kai_list = []
            flux_list = []
            Rm_ind_sol_list = []
            for job in result:
                results = job
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
                    print('RSS',':', rss, state)

            order = list(range(len(state_list)))
            order.sort(key = lambda x: kai_list[x])
            states = [state_list[x] for x in order]
            kais = [kai_list[x] for x in order]
            fluxes = [flux_list[x] for x in order]


            if self.configuration['callbacklevel'] == 3:
                print(len(fluxes), "fitting tasks in paralell processing mode were finished.")
            if self.configuration['callbacklevel'] == 4:
                print(len(fluxes), "fitting tasks in paralell processing mode were finished with finishing state",states[0],".")
            if self.configuration['callbacklevel'] >= 5:
                print(len(fluxes), "fitting tasks in paralell processing mode were finished with finishing state",states[0]," RSS of best one is: ", kais[0])

            return states, kais, fluxes

        """
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
        """

    def show_results_in_map(self, flux, mapfilename, outputfilename, formattext = ".1f"):
        """Method to project metabolic state data into metabolic file (.GML).


        Args:
            flux (dict): dict of metabolic state

            mapfilename (str): name of blank GML file

            outputfilename (str): name of uutput GML file

            format (str): format string of "format" method of string

        Examples:
            >>> model.show_results_in_map(flux_opt1, "Example_2_cancer_map_new.gml", "Example_2_cancer_map_mapped.gml")
            >>> model.show_results_in_map(flux_opt1, "Example_2_cancer_map_new.gml", "Example_2_cancer_map_mapped.gml", formattext = ".2f")




        """
        kegg_id=[]
        Id2value={}
        def idsep(ids):
            idsori = ids
            ids = re.sub("\)\s*\(",")(",ids)
            ids = re.sub("^\(","((",ids)
            ids = re.sub("\)\(","))((",ids)
            ids = re.sub("\)$","))",ids)
            #idlist = re.findall("\(\([^:]*?:[^:]*?\)\)", ids)
            idlist = []
            for idtemp in re.findall("\(\(.*?\)\)", ids):
                idtemp = re.sub("\(\(", "(", idtemp)
                idtemp = re.sub("\)\)", ")", idtemp)
                idlist.append(idtemp)
            return(idlist)


        for rid in self.reaction_ids:
            ids = self.reactions[rid]['externalids']
            for id in idsep(ids):
                kegg_id.append(id)
                Id2value[id]=flux['reaction'][rid]['value']
        for rid in self.reversible_ids:
            ids = self.reversible[rid]['externalids']
            for id in idsep(ids):
                kegg_id.append(id)
                Id2value[id]=flux['reversible'][rid]['value']
        for rid in self.metabolite_ids:
            ids = self.metabolites[rid]['externalids']
            for id in idsep(ids):
                kegg_id.append(id)
                Id2value[id]=flux['metabolite'][rid]['value']


        node="off"
        text=""
        replace="off"
        memo_off=[]
        memo=""
        kegg=[]
        label=[]

        with open(mapfilename,'r') as f:
            for line in f:
                if "  node [\n" in line:
                    node="on"

                if node=="off":
                    text=text+line
                if node=="on":
                    memo=memo+line
                    if line=="  ]\n":
                        for id in kegg_id:
                            if id in memo:
                                temp = "{:"+formattext+"}"
                                templine = temp.format(Id2value[id])
                                memo_mod=re.sub('label "\d+(\.\d+)?"',"label "+'''"'''+templine+'''"''',memo)
                                text=text+memo_mod
                                replace="on"
                        if replace=="off":
                            text=text+memo
                        node='off'
                        replace="off"
                        memo=''



        with open(outputfilename, mode='w') as f:
            f.write(text)
        f.close()

    def show_results(self, input, flux = "on", rss = "on", mdv = "on", pool_size = "off", checkrss= "off", filename = "", format = "csv"):
        """Method to output metabolic state data.

        Args:
            input (tuple): Tuple of ("name", dic of metabolic statet) [('name1', state1),('name2', state2),('name3', state3),]

            flux (str): show flux data "on"/"off"

            rss (str): show rss data "on"/"off"

            mdv (str): show mdv data "on"/"off"

            pool_size (str): show metabolite data "on"/"off"

            check rss (str): show each RSS value "on"/"off"

            filename (str): Results are saved when file name is not ""/

            format (str): "csv" or "text"

        Examples:
            >>> model.show_results(results, flux = "on", rss = "on", mdv = "on", pool_size = "off",  checkrss = "on")# Output to console
            >>> model.show_results(results, filename = "show_results.csv", format = "csv")# Output to .csv file




        """
        reaction_header = []
        rssd = []
        reaction = []
        rssdata = []
        mdv_header = []
        mdvd = []
        metabolites_header = []
        metabolites = []

        # Reaction header
        #
        reaction_header.extend(['Id', 'Reaction', "External ids"])
        reaction_header.extend([z[0] for z in input])
        reaction_header.extend(['Type','Value', 'Stdev',' lb','ub', 'Atom_mapping','Reversible'])
        #print reaction_header
        #
        # Metabolite header header
        #
        metabolites_header.extend(['Id', 'C_number', "External ids"])
        metabolites_header.extend([z[0] for z in input])
        metabolites_header.extend(['Type','Value', 'Stdev','lb','ub','Excreted','Symmetry', 'carbonsource'])
        #
        # RSS
        #
        rssd.append([])
        rssd[0].extend(['RSS', '', ""])
        for fluxd in input:
            rssd[0].append(self.calc_rss(fluxd[1]))
        rssd.append([])
        rssd.append([])
        rssd[1].extend(['Thres', '', ""])
        rssd[2].extend(['p_value', '', ""])
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
            metabolites[i].append(self.metabolites[metid]['externalids'])# 2 reaction
            metabolites[i].extend([fluxdd["metabolite"][metid]['value'] for (name, fluxdd) in input])# 4 values
            metabolites[i].append(self.metabolites[metid]['type'])# 8 type
            metabolites[i].append(self.metabolites[metid]['value'])# 9 met value (measured)
            metabolites[i].append(self.metabolites[metid]['stdev'])# 10 met std
            metabolites[i].append(self.metabolites[metid]['lb'])# 6 lb
            metabolites[i].append(self.metabolites[metid]['ub'])# 7 ub
            metabolites[i].append(self.metabolites[metid]['excreted'])# 5 atom mapping
            metabolites[i].append(self.metabolites[metid]['symmetry'])# 11 symmetry
            metabolites[i].append(self.metabolites[metid]['carbonsource'])# 12 revserible
            #
            #metabolites rss
            #
            if not self.metabolites[metid]['type'] == "fitting":
                continue
            value = self.metabolites[metid]['value']
            stdev = self.metabolites[metid]['stdev']
            fluxlist = [fluxdd["metabolite"][metid]['value'] for (name, fluxdd) in input]
            temprssd = []
            temprssd.append(metid)# 1 id
            temprssd.append(self.metabolites[metid]['C_number'])# 2 reaction
            temprssd.append(self.metabolites[metid]['externalids'])# 2 reaction
            temprssd.extend([((x-value)/stdev)**2 for x in fluxlist])# 4 values
            temprssd.append(self.metabolites[metid]['type'])# 8 type
            temprssd.append(self.metabolites[metid]['value'])# 9 met value (measured)
            temprssd.append(self.metabolites[metid]['stdev'])# 10 met std
            rssdata.append(temprssd)


        #
        # Reaction and fluxes
        #
        for i, rid in enumerate(self.reaction_ids):
            reaction.append([])
            reaction[i].append(rid)# 1 id
            reaction[i].append(self.reactions[rid]['stoichiometry'])# 2 reaction
            reaction[i].append(self.reactions[rid]['externalids'])# 3 external ids
            reaction[i].extend([fluxdd["reaction"][rid]['value'] for (name, fluxdd) in input])# 4 values
            reaction[i].append(self.reactions[rid]['type'])# 8 type
            reaction[i].append(self.reactions[rid]['value'])# 9 flux value
            reaction[i].append(self.reactions[rid]['stdev'])# 10 flux std
            reaction[i].append(self.reactions[rid]['lb'])# 6 lb
            reaction[i].append(self.reactions[rid]['ub'])# 7 ub
            reaction[i].append(self.reactions[rid]['atommap'])# 5 atom mapping
            reaction[i].append(self.reactions[rid]['reversible'])# 11 revserible#id
            #
            # RSS
            #
            if not self.reactions[rid]['type'] == "fitting":
                continue

            value = self.reactions[rid]['value']
            stdev = self.reactions[rid]['stdev']
            fluxlist =[fluxdd["reaction"][rid]['value'] for (name, fluxdd) in input]
            temprssd = []
            temprssd.append(rid)# 1 id
            temprssd.append(self.reactions[rid]['stoichiometry'])# 2 reaction
            temprssd.append(self.reactions[rid]['externalids'])# 3 external ids
            temprssd.extend([((x-value)/stdev)**2 for x in fluxlist])# 4 values
            temprssd.append(self.reactions[rid]['type'])# 8 type
            temprssd.append(self.reactions[rid]['value'])# 9 flux value
            temprssd.append(self.reactions[rid]['stdev'])# 10 flux std
            rssdata.append(temprssd)

        reversible_list = self.reversible.keys()
        length = len(reaction)
        for i, rid in enumerate(self.reversible.keys()):
            forward = self.reversible[rid]['forward']
            reverse = self.reversible[rid]['reverse']
            reaction.append([])
            reaction[length + i].append(rid)# 1 id
            reaction[length + i].append(str(forward + "<=>" + reverse))# 2 reaction
            reaction[length + i].append(self.reversible[rid]['externalids'])# 3 reaction
            reaction[length + i].extend([fluxdd["reversible"][rid]['value'] for (name, fluxdd) in input])# 4 values
            reaction[length + i].append(self.reversible[rid]['type'])# 8 type
            reaction[length + i].append(self.reversible[rid]['value'])# 9 flux value
            reaction[length + i].append(self.reversible[rid]['stdev'])# 10 flux std
            reaction[length + i].append(self.reversible[rid]['lb'])# 6 lb
            reaction[length + i].append(self.reversible[rid]['ub'])# 7 ub
            reaction[length + i].append('')# 5 atom mapping
            reaction[length + i].append('')# 11 revserible#id
            #
            # RSS
            #
            if not self.reversible[rid]['type'] == "fitting":
                continue

            value = self.reversible[rid]['value']
            stdev = self.reversible[rid]['stdev']
            fluxlist =[fluxdd["reversible"][rid]['value'] for (name, fluxdd) in input]

            forward = self.reversible[rid]['forward']
            reverse = self.reversible[rid]['reverse']
            temprssd = []
            temprssd.append(rid)# 1 id
            temprssd.append(str(forward + "<=>" + reverse))# 2 reaction
            temprssd.append(self.reversible[rid]['externalids'])# 3 reaction
            temprssd.extend([((x-value)/stdev)**2 for x in fluxlist])# 4 values
            temprssd.append(self.reversible[rid]['type'])# 8 type
            temprssd.append(self.reversible[rid]['value'])# 9 flux value
            temprssd.append(self.reversible[rid]['stdev'])# 10 flux std
            rssdata.append(temprssd)

        #print reaction
        #
        #
        # MDVs header
        #
        mdv_header.extend(['Experiment', 'Fragment_Num', "Time"])
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
                    #
                    # RSS
                    #
                    if not self.experiments[ex_id]['mdv_use'][i] == 1:
                        continue

                    value = self.experiments[ex_id]['mdv_exp_original'][i]
                    stdev = self.experiments[ex_id]['mdv_std_original'][i]
                    fluxlist =[mdvdd[1][fragment][number] for mdvdd in mdv_data]

                    temprssd = []
                    temprssd.append(ex_id) # 1 ex
                    temprssd.append(fragment_number) # 2 fragment
                    temprssd.append(int(0)) # 3 Timecourse
                    temprssd.extend([((x-value)/stdev)**2 for x in fluxlist])
                    temprssd.append(self.experiments[ex_id]['mdv_use'][i]) # 6 use or not
                    temprssd.append(self.experiments[ex_id]['mdv_exp_original'][i]) # 7 ratio
                    temprssd.append(self.experiments[ex_id]['mdv_std_original'][i]) # 7 stdev
                    rssdata.append(temprssd)

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
                    # RSS
                    #
                    if not self.experiments[ex_id]['mdv_use'][j] == 1:
                        continue
                    value = self.experiments[ex_id]['mdv_exp_original'][j]
                    stdev = self.experiments[ex_id]['mdv_std_original'][j]
                    fluxlist =[mdvdd[1][fragment][timepoint_id[timepoint_list[j]]][number] for mdvdd in mdv_data]

                    temprssd = []
                    temprssd.append(ex_id) # 1 ex

                    temprssd.append(fragment_number) # 2 fragment
                    temprssd.append(timepoint_list[j])   # 3 timepoint
                    temprssd.extend([((x-value)/stdev)**2 for x in fluxlist])
                    temprssd.append(ex_id) # 1 ex
                    temprssd.append(self.experiments[ex_id]['mdv_use'][j]) # 6 use or not
                    temprssd.append(self.experiments[ex_id]['mdv_exp_original'][j]) # 7 ratio
                    temprssd.append(self.experiments[ex_id]['mdv_std_original'][j]) # 7 stdev
                    rssdata.append(temprssd)


        #
        #
        # output results
        #
        #
        if len(filename) > 0:
            import csv

            with open(filename, 'w', newline='') as f:
                if format == "csv":
                    writer = csv.writer(f, dialect='excel')
                else:
                    writer = csv.writer(f, delimiter='\t')
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
                if checkrss == "on":
                    writer.writerow(mdv_header)
                    writer.writerows(rssdata)


        else:
            text = ""
            if rss == "on":
                #print reaction header
                text = text + "{0:15.15s}".format(reaction_header[0])
                text = text + "{0:25.25s}".format(reaction_header[1])
                text = text + "{0:10.10s}".format(reaction_header[2])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(reaction_header[3+i]))
                text = text + "\n"
                for data in rssd:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.25s}".format(data[1])
                    text = text + "{0:10.10s}".format(data[2])
                    for i in range(len(input)):
                        text = text + "{0:>8.2f}".format(data[3+i])
                    text = text + "\n"

            if flux == "on":
                #print flux data
                text = text + "{0:15.15s}".format(reaction_header[0])
                text = text + "{0:25.25s}".format(reaction_header[1])
                text = text + "{0:10.10s}".format(reaction_header[2])
                step = 2
                for i in range(len(input)):
                    step = step + 1
                    text = text + "{0:>8.7s}".format(str(reaction_header[step]))
                text = text + " " + "{0:8.7s}".format(reaction_header[step + 1])
                text = text + "{0:>6.5s}".format(reaction_header[step + 2])
                text = text + "{0:>6.5s}".format(reaction_header[step + 3])
                text = text + "{0:>6.5s}".format(reaction_header[step + 4])
                text = text + "{0:>6.5s}".format(reaction_header[step + 5])
                text = text + " " + "{0:25.25s}".format(reaction_header[step + 6])


                text = text + "\n"
                for data in reaction:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.25s}".format(data[1])
                    text = text + "{0:10.10s}".format(data[2])
                    step = 2
                    for i in range(len(input)):
                        step = step + 1
                        text = text + "{0:>8.1f}".format(data[step])
                    text = text + " " + "{0:8.7s}".format(data[step + 1])
                    text = text + "{0:6.1f}".format(data[step + 2])
                    text = text + "{0:6.1f}".format(data[step + 3])
                    text = text + "{0:6.1f}".format(float(data[step + 4]))
                    text = text + "{0:6.1f}".format(float(data[step + 5]))
                    text = text + " " + "{0:25.25s}".format(data[step + 6])
                    text = text + "\n"

            if pool_size == "on":
                #print reaction header
                text = text + "{0:15.15s}".format(metabolites_header[0])
                text = text + "{0:25.5s}".format(metabolites_header[1])
                text = text + "{0:10.10s}".format(metabolites_header[2])
                step = 2
                for i in range(len(input)):
                    step = step + 1
                    text = text + "{0:>8.7s}".format(str(metabolites_header[step]))
                text = text + " "+"{0:8.7s}".format(metabolites_header[step + 1])
                text = text + "{0:>6.6s}".format(metabolites_header[step + 2])
                text = text + "{0:>6.6s}".format(metabolites_header[step + 3])

                text = text + " "+"{0:6.6s}".format(metabolites_header[step + 4])
                text = text + " "+"{0:6.6s}".format(metabolites_header[step + 5])
                text = text + " "+"{0:9.9s}".format(metabolites_header[step + 6])
                text = text + " "+"{0:9.9s}".format(metabolites_header[step + 7])
                text = text + " "+"{0:9.9s}".format(metabolites_header[step + 8])
                text = text + "\n"
                for data in metabolites:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:<25.0f}".format(data[1])
                    text = text + "{0:10.10s}".format(data[2])
                    step = 2
                    for i in range(len(input)):
                        step = step + 1
                        text = text + "{0:>8.3f}".format(data[step])
                    text = text + " "+"{0:8.7s}".format(data[step + 1])
                    text = text + "{0:6.1f}".format(data[step + 2])
                    text = text + "{0:6.1f}".format(data[step + 3])

                    text = text + " "+"{0:6.3f}".format(float(data[step + 4]))
                    text = text + " "+"{0:6.1f}".format(float(data[step + 5]))
                    text = text + " "+"{0:9.9s}".format(data[step + 6])
                    text = text + " "+"{0:9.9s}".format(data[step + 7])
                    text = text + " "+"{0:9.9s}".format(data[step + 8])

                    text = text + "\n"

            if mdv == "on":
                text = text + "{0:15.15s}".format(mdv_header[0])
                text = text + "{0:25.10s}".format(mdv_header[1])
                text = text + "{0:10.5s}".format(mdv_header[2])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(mdv_header[3+i]))
                text = text + "{0:>8.5s}".format(str(mdv_header[-3]))
                text = text + " " + "{0:>6.6s}".format(mdv_header[-2])
                text = text + "{0:>6.6s}".format(mdv_header[-1])
                text = text + "\n"
                for data in mdvd:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.15s}".format(data[1])
                    text = text + "{0:<10.2f}".format(data[2])
                    for i in range(len(input)):
                        text = text + "{0:8.4f}".format(data[3+i])
                    text = text + "{0:>8.5s}".format(str(data[-3]))
                    text = text + " " + "{0:6.4f}".format(data[-2])
                    text = text + "{0:6.2f}".format(data[-1])
                    text = text + "\n"
            if checkrss == "on":
                text = text + "{0:15.15s}".format(mdv_header[0])
                text = text + "{0:25.10s}".format(mdv_header[1])
                text = text + "{0:10.5s}".format(mdv_header[2])
                for i in range(len(input)):
                    text = text + "{0:>8.7s}".format(str(mdv_header[3+i]))
                text = text + "{0:>8.5s}".format(str(mdv_header[-3]))
                text = text + " " + "{0:>6.6s}".format(mdv_header[-2])
                text = text + "{0:>6.6s}".format(mdv_header[-1])
                text = text + "\n"
                for data in rssdata:
                    text = text + "{0:15.15s}".format(data[0])
                    text = text + "{0:25.15s}".format(data[1])
                    text = text + "{0:<10.2s}".format(str(data[2]))
                    for i in range(len(input)):
                        text = text + "{0:8.2f}".format(data[3+i])
                    text = text + "{0:>8.5s}".format(str(data[-3]))
                    text = text + " " + "{0:6.4f}".format(data[-2])
                    text = text + "{0:6.2f}".format(data[-1])
                    text = text + "\n"
            print(text)



    def show_flux_balance(self, results, metabolite, filename = "", format = "csv"):
        """List of metabolic flux levels of reactions related to given metabolite.

        Args:
            input: Tuple of ("name", dic of metabolic statet) [('name1', state1),('name2', state2),('name3', state3),]

            metabolite: ID of metabolite of interest

            filename: Results are saved when file name is not "".

            format: "csv" or "text"

        Examples:
            >>> model.show_metabolite_balance(results, "Fum")
            >>> model.show_metabolite_balance(results, "Fum", filename = "metabolite_balance.csv", format = "csv")



        """

        used_reaction = {}
        result_hash = {}
        result_hash_for_screen = {}
        for id in self.reversible:
            for text, flux in results:
                forward = self.reversible[id]['forward']
                reverse = self.reversible[id]['reverse']
                string = forward+"+"+reverse
                sumtemp = 0
                check = 0
                for reaction_id in string.split("+"):
                    if not reaction_id in self.reaction_ids:
                        continue
                    if metabolite in self.reactions[reaction_id]['stoichiometry_metabolite_list']:
                        check = check + 1
                        used_reaction[reaction_id] = 1
                        coefficient= self.reactions[reaction_id]['stoichiometry_metabolite_list'][metabolite]
                        #print(metabolite, id, reaction_id, coefficient, text, flux['reaction'][reaction_id]['value'])
                        sumtemp = sumtemp + coefficient * flux['reaction'][reaction_id]['value']
                if check > 0:
                    if not id in result_hash:
                        result_hash[id] = []
                        result_hash_for_screen[id]=[]
                    if sumtemp > 0:
                        result_hash[id].extend([sumtemp, 0])
                    else:
                        result_hash[id].extend([0, sumtemp*-1])
                    result_hash_for_screen[id].append(sumtemp)

        for id in self.reactions:
            if id in used_reaction:
                continue
            if not id in self.reaction_ids:
                continue
            if self.reactions[id]["type"] == "pseudo":
                continue
            #used_reaction[reaction_id] = 1
            if metabolite in self.reactions[id]['stoichiometry_metabolite_list']:
                if not id in result_hash:
                    result_hash[id] = []
                    result_hash_for_screen[id] = []
                coefficient= self.reactions[id]['stoichiometry_metabolite_list'][metabolite]
                for text, flux in results:
                    tempvalue = coefficient * flux['reaction'][id]['value']
                    #print(metabolite, id, coefficient, text, tempvalue)
                    if tempvalue > 0:
                        result_hash[id].extend([tempvalue, 0])
                    else:
                        result_hash[id].extend([0, tempvalue*-1])
                    result_hash_for_screen[id].append(tempvalue)
        if filename == "":

            text = ""
            text = text + "{0:10.10s}".format("Reaction")
            for texttemp, flux in results:
                text = text + "{0:10.10s}".format(texttemp)
            text = text + "\n"

            for id, values in result_hash_for_screen.items():
                text = text+ "{0:10.10s}".format(id)
                for value in values:
                    text = text+ "{0:>10.2f}".format(value)
                text = text + "\n"
            print(text)
        else:

            output = []
            header = ["Reaction"]
            for text, flux in results:
                header.append(text+"_production")
                header.append(text+"_consumption")
            output.append(header)
            for id, value in result_hash.items():
                output.append([id] +value)

            try:
                with open(filename, 'w', newline='') as f:
                    import csv
                    if format == "text":
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerows(output)
                    else:
                        writer = csv.writer(f, dialect='excel')
                        writer.writerows(output)
                return True
            except:
                return False





    def calc_rss(self, fluxes, mode = "flux"):
        """Getter to calc residual sum of square (RSS) level between estimated MDVs of 'flux' and measured MDVs in 'experiment(s)'

        Args:
            flux (dict): Dictionary of metabolic flux distribution or independent flux vector

            mode (str):
                * "independent":Calculation of RSS from independent flux vector
                * "flux":Calculation of RSS from state dictionary

        Returns:
            float : Residual sum of square (RSS)

        Examples:
            >>> rss = model.calc_rss(flux_opt)
            >>> print rss
            7564.123

        History:
            Call of calmdv funcion is removed.12/9/2014


        """
        #
        # Check experiment
        #
        if len(self.experiments.keys()) == 0:
            if self.configuration['callbacklevel'] >= 1:
                print('No experiment was set to the modelnames.')
            return False
        if fluxes == []:
            if self.configuration['callbacklevel'] >= 1:
                print('Metabolic state data is required.')
            return False
        if mode == "independent":
            fluxlist = [fluxes]
        else:
            if type(fluxes) == list:
                fluxlist = fluxes
            else:
                fluxlist = [fluxes]

        rss_list = []

        for flux in fluxlist:
            # zero independent flux
            if mode == "independent":
                Rm_ind = flux
            else:
                if type(flux) == list:
                    flux = flux[0]

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
                diffmdv = self.func["diffmdv"],
                callbacklevel = self.configuration['callbacklevel']
                )
            rss_list.append(rss)
            if self.configuration['callbacklevel'] >= 5:
                print("RSS of given state(s): ",rss_list)
        if mode == "independent":
            return rss_list[0]
        if type(fluxes) == list:
            return rss_list
        else:
            return rss_list[0]

    def goodness_of_fit(self, flux, alpha = 0.05):
        """Getter to calculate goodness-of-fit of a given flux distribution

        Args:
            flux (dict): Dictionary of metabolic flux distribution

            alpha (float): Confidence level. Default is 0.05.

        Returns:
            * pvalue (float) p-value determined by chi-squere test.
            * rss_thres (float) Threshold RSS at the alpha levels.

        Examples:
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

        if self.configuration["callbacklevel"] >= 4:
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
        """Getter to obtain a threshold level of RSS for searching confidence interval

        Args:
            flux (dict): Dictionary of best fitted flux distribution.

            alpha (float): confidence level. For 95% confience interval, please set alpha = 0.05

            dist (str): Probability distribution for test
                * 'F-dist' (default): Using F-distribution
                * 'chai-dist':Using Chai-square distribution

        Returns:
            * thres (float) a threshold level of RSS for searching confidence interval

            * number_of_measurements (int) number of measurements

            * degree_of_freedom (int) degree of freedom

        Examples:
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
        """Getter of degree of freedom of the model

        Examples:
            >>> degree_of_freedom = model.get_degree_of_freedom()

        Returns:
            int: degree of freedom.

        """
        #modified 056 200517
        if self.experiments[sorted(self.experiments.keys())[0]]["mode"] == "ST":
            free_metabolites = 0
            for i in range(len(self.metabolite_ids)):
                if self.metabolites[self.metabolite_ids[i]]["type"] == "free":
                    free_metabolites = free_metabolites + 1
            p = self.numbers["independent_number"] - free_metabolites
            return p
        if self.experiments[sorted(self.experiments.keys())[0]]["mode"] == "INST":
            return self.numbers["independent_number"]
    def get_number_of_independent_measurements(self):
        """Getter of number of independent measurements of the model

        Returns:
            int: number of independent measurements.

        History:
            29/9/2014 Modified to consider "fitting" flux information.

        Example:
            >>> number_of_independent_measurements = model.get_number_of_independent_measurements()

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

    def add_perturbation(self, flux):
        """New state dict is generarate from a given state dict by an addition of small random perturbation to an independent flux


        Args:
            state_dict (dict): Dictionary of metabolic state data.

        Returns:
            * check (boolean)
            * perturbed_state (dict)Dictionary of metabolic state data.
            * flux_vector (array) List of metabolic flux, metabolite concentration data.
            * reaction_name_list (array) List of ids

        Examples:
            >>> check, perturbed_state, flux_vector, reaction_name_list = model.add_perturbation(state_dict)

        History
            Newly created at 21/1/2021

        """
        independent_vector, lbi, ubi, independent_list = self.get_independents(flux)
        #print(independent_vector, lbi, ubi, independent_list)
        for j in range(10000):
            reac_num = numpy.random.randint(0, len(independent_vector))
            Rm_ind_next_temp = independent_vector[:]
            perturbation =  (numpy.random.rand() - 0.5) * 2 * (ubi[reac_num]-lbi[reac_num])/100
            Rm_ind_next_temp[reac_num] += perturbation
            if lbi[reac_num] < Rm_ind_next_temp[reac_num] < ubi[reac_num]:
                check, flux_vector, reaction_name = self.check_independents(Rm_ind_next_temp)
                if check == True:
                    break
        else:
            #
            # Failed
            #
            check, flux_vector, reaction_name = self.check_independents(independent_vector)
            return(False, flux, flux_vector, reaction_name)
        perturbed_state = self.generate_state_dict(flux_vector)
        return(True, perturbed_state, flux_vector, reaction_name)




    def get_independents(self, flux):
        """Getter of independent status (flux) vector, its lower and upper boundaries and a list of independent flux ids of given state (flux) dict.

        Examples:
            >>> independent_vector, lb, ub, independent_list = model.get_independents(state_dict)

        Args:
            flux (dict): Dictionary of metabolic state data.

        Returns:
            * independent_vector (array) list of independent flux values

            * lb (array) list of independent flux values

            * ub (array) list of independent flux values

            * independent_list (array) list of independent fluxes, metabolites, and reversible reactions.


        History:
            Newly created at 21/1/2021

        """
        Rm_ind = [flux[group][id]["value"] for (group, id) in self.vector['independent_flux']]

        return Rm_ind, self.vector['independent_lb'], self.vector['independent_ub'], self.vector['independent_flux']

    def check_independents(self, independents):
        """Checker method to test whether a given vector of independent fluxes is able to produce feasible status (flux) vector within lower and upper boundaries.

        Args:
            independents (array): list of independent fluxes generated by self.get_independents

        Returns:
            * check (Boolean) True means that a given vector of independent fluxes is able to produce feasible status (flux) vector within lower and upper boundaries.

            * flux_vector (array) List of metabolic flux, metabolite concentration data.

            * reaction_name (array) List of ids


        Examples:
            >>> check, flux_vector, reaction_name = self.check_independents(independents)

        History:

            Newly created at 21/1/2021

        """
        callbacklevel = self.configuration['callbacklevel']
        stoichiometric_num = self.numbers['independent_start']
        reaction_num= self.numbers['total_number']
        matrixinv=self.matrixinv
        Rm_initial= self.vector["Rm_initial"]
        lb = copy.copy(self.vector["lb"])
        ub = copy.copy(self.vector["ub"])

        Rm = numpy.array(list(Rm_initial))
        Rm[stoichiometric_num:reaction_num] = list(independents)
        tmp_r = numpy.dot(matrixinv, Rm)
        temp = sum([1 for x in numpy.array(lb) - tmp_r if x > 0]) + sum([1 for x in tmp_r - numpy.array(ub) if x > 0])
        reaction_name = [id for (group, id) in self.vector["ids"]]
        if temp > 0:
            return False, tmp_r, reaction_name
        return True, tmp_r, reaction_name


    def posterior_distribution(self, state_dict, number_of_steps = 1000000):
        """Generator of posterior distribution by the Metropolis-Hastings algorism

        Args:
            state_dict (dist): Dictionary of initial metabolic state

            number_of_steps (int): Length of Markov-chain
        Returns
            array: 2d array. Markov-chain of metabolic state vector

        Examples:
            >>> record = model.posterior_distribution(perturbed_state, number_of_steps = numberofsteps)

        History:

            Newly created at 21/1/2021
        """
        callbacklevel = self.configuration['callbacklevel']

        if 'ncpus' in self.configuration:
            numberofchains = self.configuration['ncpus']
        else:
            numberofchains = 1
        #
        # tuple of all parallel python servers to connect with
        #
        try:
            from joblib import Parallel, delayed
        except:
            print("This function requires joblib")
            return False

        def calcrsss(Rm_ind_original, rss_original, parameters, func):
            """Low level function for Metropolic Hasting method for joblib parallel execution.

            Args:
                Rm_ind_original (array): vector of independent flux of given initial state.

                rss_original (float): RSS of given initial state.

                parameters (dict): dict of parameters.

                func (dict): functions used in this method.

            Returns:
                float: RSS + Penalty score (When out side of the lower and upper boundaries)


            See Also:]

                fit_r_mdv_scipy

            """
            import mkl

            mkl.set_num_threads(1)
            Rm_initial = parameters['Rm_initial']
            stoichiometric_num = parameters['stoichiometric_num']
            reaction_num = parameters['reaction_num']
            reac_met_num = parameters['reaction_num']
            matrixinv = parameters['matrixinv']
            experiments = parameters['experiments']
            mdv_exp = numpy.array(parameters['mdv_exp'])
            mdv_use = parameters['mdv_use']
            covinv = parameters['covinv']
            lb = parameters['lb']
            ub = parameters['ub']
            lbi =parameters['lbi']
            ubi = parameters['ubi']
            df = parameters['df']
            rss = rss_original
            Rm_ind = Rm_ind_original[:]

            if isinstance(func, dict):
                calmdv = func["calmdv"]
                diffmdv = func["diffmdv"]
            else:
                locals_dic = locals()
                exec(func, globals(), locals_dic)
                calmdv = locals_dic["calmdv"]
                diffmdv = locals_dic["diffmdv"]

            Length_ind = len(Rm_ind)

            pdf = scipy.stats.chi2.pdf(x = rss, df = df)
            accept = 0
            acceptd = 0
            hamatta = 0

            while accept < 1000:
                if hamatta > 100:
                    return Rm_ind, rss, -1



                Rm_ind_next = Rm_ind[:]
                dame = 0
                for i in range(3):
                #for reac_num in range(Length_ind):
                    reac_num = numpy.random.randint(0, Length_ind)
                    for j in range(10):
                        Rm_ind_next_temp = Rm_ind_next[:]
                        perturbation =  (numpy.random.rand() - 0.5) * 2 * (ubi[reac_num]-lbi[reac_num])/100
                        if lbi[reac_num] < Rm_ind_next[reac_num] + perturbation < ubi[reac_num]:
                            Rm_ind_next_temp[reac_num] = Rm_ind_next[reac_num] + perturbation
                            Rm = numpy.array(list(Rm_initial))
                            Rm[stoichiometric_num: reaction_num] = list(Rm_ind_next_temp)
                            tmp_r = numpy.dot(matrixinv, Rm)
                            temp = sum([1 for x in numpy.array(lb) - tmp_r if x > 0]) + sum([1 for x in tmp_r - numpy.array(ub) if x > 0])
                            if temp > 0:
                                continue
                            Rm_ind_next[reac_num] = Rm_ind_next[reac_num] + perturbation
                            break
                    else:
                        dame = dame + 1
                if dame == 3:
                    #print("hamatta")
                    hamatta = hamatta + 1
                    continue



                Rm = numpy.array(list(Rm_initial))
                Rm[stoichiometric_num: reaction_num] = list(Rm_ind_next)
                tmp_r = numpy.dot(matrixinv, Rm)

                mdv_original = list(tmp_r)

                for experiment in sorted(experiments.keys()):
                    target_emu_list = experiments[experiment]['target_emu_list']
                    mdv_carbon_sources = experiments[experiment]['mdv_carbon_sources']
                    #
                    mdv_original_temp, mdv_hash = calmdv(list(tmp_r), target_emu_list, mdv_carbon_sources)

                    mdv_original.extend(mdv_original_temp)

                mdv = numpy.array([y for x, y in enumerate(mdv_original) if mdv_use[x] != 0])
                res = mdv_exp - mdv
                f = numpy.dot(res, numpy.dot(covinv, res))

                rss_next = f #+ sum


                pdf_next = scipy.stats.chi2.pdf(x = rss_next, df = df)
                #
                # If probabirity of next point is smaller than that of present point
                #
                if pdf_next < pdf:
                    if numpy.random.rand() > pdf_next/pdf:
                            acceptd = acceptd + 1
                            accept = accept + 1
                            continue
                Rm_ind[:] = Rm_ind_next[:]
                rss = rss_next
                pdf = pdf_next
                accept = accept + 1
                #print("accpeted", accept, acceptd, rss_next, pdf_next, rss, pdf)

            return (Rm_ind, rss, acceptd)

        #
        # Creation of header
        #

        recordini = []

        tmp_r = ['RSS']
        tmp_r.extend([id for (group, id) in self.vector["ids"]])
        recordini.append(tmp_r)
        record = recordini[:]


        rss = self.calc_rss(state_dict)
        df = 0
        for experiment in sorted(self.experiments.keys()):
            df += self.experiments[experiment]['number_of_measurement']


        sf = scipy.stats.chi2.sf(x = rss, df = df)
        pdf = scipy.stats.chi2.pdf(x = rss, df = df)
        if callbacklevel >= 1:
            print("RSS:",rss, "p-value of chisq test", sf, "pdf", pdf)
        if sf < 0.95:
            if callbacklevel >= 1:
                print("The metabolic model failed to overfit to isotopoer data")

        #
        # Preparation of parameters
        #"
        Rm_ind = [state_dict[group][id]["value"] for (group, id) in self.vector['independent_flux']]
        lbi = [state_dict[group][id]["lb"] for (group, id) in self.vector['independent_flux']]
        ubi = [state_dict[group][id]["ub"] for (group, id) in self.vector['independent_flux']]
        stoichiometric_num = self.numbers['independent_start']
        reaction_num=  self.numbers['total_number']
        Rm_initial= self.vector["Rm_initial"]
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
                if callbacklevel >= 1:
                    print("Error in ", i, std)
            covinv[i,i] = 1.0/(std**2)

        parameters ={"reaction_num": self.numbers['total_number'],
            "stoichiometric_num": self.numbers['independent_start'],
            "matrixinv":self.matrixinv,
            "experiments":self.experiments,
            "mdv_exp":mdv_exp,
            "mdv_use":mdv_use,
            "covinv":covinv,
            "Rm_initial": self.vector["Rm_initial"],
            "lb" : copy.copy(self.vector["lb"]),
            "ub" : copy.copy(self.vector["ub"]),
            "reac_met_number" : self.numbers['reac_met_number'],
            "lbi" : lbi,
            "ubi" : ubi,
            "df" : df
        }

        Inds = []
        RSSs = []
        for ii in range(numberofchains):
            Inds.append(Rm_ind)
            RSSs.append(rss)



        for ti in range(int(number_of_steps/1000)):
            Indsnext = []
            RSSsnext = []
            #
            #
            #
            r = Parallel(n_jobs=-1)([delayed(calcrsss)(Inds[i],RSSs[i], parameters, self.calmdv_text) for i in range(numberofchains)] )
            for results in r:
                Rm_ind, rss, acceptd = results
                if callbacklevel >= 1:
                    print(str(ti), "RSS:",rss, scipy.stats.chi2.pdf(x = rss, df = df), acceptd, df)
                if acceptd == -1:
                    t_num = numpy.random.randint(0,numberofchains)
                    Rm_ind = Inds[t_num]
                    rss = RSSs[t_num]
                    if callbacklevel >= 1:
                        print("Failed to find proposal flux>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                check, flux_vector, reaction_name = self.check_independents(Rm_ind)

                if check == False:
                    t_num = numpy.random.randint(0,numberofchains)
                    Rm_ind = Inds[t_num]
                    rss = RSSs[t_num]
                    if callbacklevel >= 1:
                        print("flux is out of range>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                Indsnext.append(Rm_ind)
                RSSsnext.append(rss)
                check, flux_vector, reaction_name = self.check_independents(Rm_ind)

                tmp_record = [rss]
                tmp_record.extend(flux_vector)

                record.append(tmp_record)
            Inds = Indsnext[:]
            RSSs = RSSsnext[:]

        return record

    def generate_ci_templete(self, targets = "normal"):
        return self.generate_ci_template(targets)


    def generate_ci_template(self, targets = "normal"):
        """Generator of a template dictionary to keep confidence interval search results.

        The templatedisctionary generated by this method is used for:
            * Setting reactions, metabolites, and reversible reactions for CI searching in self.search_ci method.
            * Store all record of CI searching by self.search_ci method.

        Args:
            targets (str):
                * 'normal' Free reversible and inreversible reactions
                * 'independent' Independent reactions
                * 'all' all reactions
                * 'without_reversible_reactions' inreversible reactions

        Returns:
            dict: template dictionary


        Examples:
            >>> ci = model.generate_ci_templete(targets = [("reaction","r29_g6pdh"),("reaction","r27_pc"),("reaction","r28_mae")])
            >>> ci = model.generate_ci_templete(targets = 'normal')

        See Also:
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



    def search_ci(self, ci, flux, method = 'grid', alpha = 0.05, dist = 'F_dist', outputthres = 2.0):
        """Method to estimate confidence intervals using information in a 'ci' dictionary.

        Used parameters in self.configuration().

            * 'grid_search_iterations' (default = 1): MNumber of trials for model fitting at each point in grid search. Example: model.set_configuration(grid_search_iterations = 1)

            * 'number_of_repeat'(default = 50): Number of initial metabolic flux disributions generated for each trial in grid search.Best initial metabolic flux disribution with smallest RSS is employed.

            * 'number_of_repeat': Number of repeated execution by 'deep' functions. Example: model.set_configuration(number_of_repeat = 2)

            * 'ncpus': Number of CPUs used in parallel processing. Example: model.set_configuration(ncpus = 16)


        Argss:
            ci (dict): dictionary for storeing confidence interval results. ci is generated by self.generate_ci_templete()

            flux (dict): dict of metabolic state (best fitted)

            method (str): 'grid' (Grid search method) is available

            alpha (float): (default 0.05) Confidence interval level. 0.05 means a 95% confidence interval level

            dist (str):
                * 'F-dist' (default): Using F-distribution
                * 'chai-dist':Using Chai-square distribution

            outputthres (float): The grid search results larger than outputthres * threshold rss were removed from output

        Returns:
            dict: dictionary of confidence interval results.

        Examples:
            >>> ci = model.generate_ci_templete(targets = 'independent')
            >>> ci = model.search_confidence_interval_parallel(ci, flux, method = 'grid')



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


        #if (callbacklevel >= 2):
        #    print("Number of active nodes is " + str(job_server.get_active_nodes()))
        ci['record']['flux'] = flux
        ci['record']['thres'] = thres
        ci['record']['rss'] = self.calc_rss(flux)



        #######################################################
        #
        # Grid search
        #
        ########################################################
        if method == "gridpp":
            #
            # Check Parallel python
            #
            try:
                import pp
                #job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
            except:
                print("This function requires Parallel Python!")
                return False

            starttime=time.perf_counter()
            step = 0.0
            data_tmp={}
            job_number = self.configuration["grid_search_iterations"]
            initial_search_repeats_in_grid_search = self.configuration["initial_search_repeats_in_grid_search"]
            #
            # confidence interval
            #
            rss_bestfit = self.calc_rss(flux)

            for group, rid in sorted(ci['data'].keys()):
                if ci['data'][(group, rid)]['use'] != 'on': continue
                flux_opt_rid = flux[group][rid]['value']


                data={(group, rid):{"flux_data":[flux_opt_rid],"rss_data":[rss_bestfit],"state":["Best fit"],"log":{},"raw_flux_data":[[flux["reaction"][i]["value"] for i in self.reaction_ids]]}}
                data_tmp.update(data)

                flux_lower =ci['data'][(group, rid)]["lower_boundary"]
                flux_upper =ci['data'][(group, rid)]["upper_boundary"]



                # Grid number is 20
                if (flux_upper-flux_lower)>=100:
                    n=15
                else:
                    n=10

                if callbacklevel >= 1:
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
                    if fixed_flux >= flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux <= flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux,flux_lower,flux_upper, "interation", i)

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallelpp")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallelpp")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag > 0:
                        counter_of_missed_initial_state = 0
                    else:
                        counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
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
                    if fixed_flux >= flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux <= flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux,flux_lower,flux_upper, "interation", i)

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallelpp")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallelpp")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag > 0:
                        counter_of_missed_initial_state = 0
                    else:
                        counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
                #
                # Pitch jobs to pp
                #
                jobs=[]
                job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in reversed(temp_array_initial_fluxes):
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    try:
                        rss = self.calc_rss(flux_opt)
                    except:
                        if callbacklevel >= 3:
                            print("Skipped:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)
                        continue

                    if callbacklevel >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)

                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                    jobs.append([(group, rid),fixed_flux, flux_opt, job_server.submit(optimize.fit_r_mdv_deep, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                    if callbacklevel >= 2:
                        print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)

                #
                # Return to original position
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()
                if callbacklevel >= 2:
                    print("Waiting for pp response", rid, "flux_opt ")
                #
                # Retrive results
                #
                for (group, rid), flux_value, flux_opt, job in jobs:
                    print(group, rid)
                    results = job()
                    if results == None:
                        #continue
                        #
                        # Large value is used when falied
                        #
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        #opt_flux = []
                        state = "Failed in optimization"
                    else:
                        state, rss, opt_flux, Rm_ind_sol = results
                        state =  "Finished successfully "
                    #
                    #  Large value is used when falied
                    #
                    if len(opt_flux) == 0:
                        #continue
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        opt_flux = []
                        state = "Failed in optimization"
                    #
                    #  Store data
                    #
                    data_tmp[(group, rid)]["flux_data"].append(flux_value)
                    data_tmp[(group, rid)]["rss_data"].append(rss)
                    data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)
                    data_tmp[(group, rid)]["state"].append(state) #★200517追加　落ちてないのかfialしたのか判断できないので固定値にした

                    if callbacklevel >= 3:
                        print("Finished:", group, rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, "flux_fixed", flux_value, "rss", rss, "state", state)
                job_server.destroy()


            #
            # Cals again around thres
            #
            for group, rid in sorted(ci['data'].keys()):
                if ci['data'][(group, rid)]['use'] != 'on': continue
                print(group, rid)

                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']

                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[0])
                #
                # Add right and left end points
                #
                interval = flux_upper - flux_lower # Added at 200517 to exactly detect edges
                data_sorted.append((flux_upper + interval * 0.0001, thres * 11.0))
                data_sorted.insert(0, (flux_lower - (interval * 0.0001), thres * 11.0))
                #
                # points under the threshold
                #
                below_threshold  = [x for x in data_sorted if x[1] < thres]
                # Detect just before the lower boundary
                flux_previous_left = below_threshold[0][0]
                # Detect just before the upper boundary
                flux_previous_right = below_threshold[-1][0]

                # Detect just next to the boundary
                #print(data_sorted, flux_previous_left, flux_previous_right)
                flux_next_left  = max([x[0] for x in data_sorted if x[0] < flux_previous_left])
                flux_next_right  = min([x[0] for x in data_sorted if x[0] > flux_previous_right])
                #print("flux_next_left ",flux_next_left,", flux_next_right ",flux_next_right)
                data_tmp[(group, rid)]["log"]["1st previous_left"] = flux_previous_left
                data_tmp[(group, rid)]["log"]["1st previous_right"] = flux_previous_right
                data_tmp[(group, rid)]["log"]["1st next_left"] = flux_next_left
                data_tmp[(group, rid)]["log"]["1st next_right"] = flux_next_right

                if callbacklevel >= 2:
                    print("flux_previous_left ",flux_previous_left,", flux_previous_right ",flux_previous_right)
                    print("flux_next_left ",flux_next_left,", flux_next_right ",flux_next_right)
                #
                # This part should be reconsiderd 200517
                #
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
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e)

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallelpp")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallelpp")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag == 0:
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
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
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e)

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallelpp")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallelpp")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 2:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag == 0:
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
                #
                # Pitch jobs to pp
                #
                jobs=[]
                job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in temp_array_initial_fluxes:
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()

                    if callbacklevel >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)

                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    functions = (optimize.calc_MDV_residue_scipy, optimize.fit_r_mdv_scipy,optimize.calc_MDV_residue_nlopt, optimize.fit_r_mdv_nlopt)
                    jobs.append([(group, rid),fixed_flux, flux_opt, job_server.submit(optimize.fit_r_mdv_deep, parameters, functions,("numpy","nlopt","scipy","scipy.integrate"))])
                    if callbacklevel >= 2:
                        print("New job was added to pp:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)
                #
                # Return to original position [need modificaton]
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()


                #
                # Retrieve results
                #
                for (group, rid), flux_value, flux_opt, job in jobs:
                    print(group, rid)
                    results = job()
                    if results == None:
                        #continue
                        #
                        # Large value is used when falied
                        #
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        #opt_flux = []
                        state = "Failed in optimization"
                    else:
                        state, rss, opt_flux, Rm_ind_sol = results
                        state =  "Finished successfully "
                    #
                    #  Large value is used when falied
                    #
                    if len(opt_flux) == 0:
                        #continue
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        opt_flux = []
                        state = "Failed in optimization"
                    #
                    #  Store data
                    #
                    data_tmp[(group, rid)]["flux_data"].append(flux_value)
                    data_tmp[(group, rid)]["rss_data"].append(rss)
                    data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)
                    data_tmp[(group, rid)]["state"].append(state) #★200517追加　落ちてないのかfialしたのか判断できないので固定値にした

                    if callbacklevel >= 3:
                        print("Finished:", group, rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, "flux_fixed", flux_value, "rss", rss, "state", state)
                job_server.destroy()

            if callbacklevel >= 1:
                print("Data collection was finished")

            for (group, rid) in ci['data'].keys():
                if ci['data'][(group, rid)]['use'] != 'on':
                    continue
                #
                # Retirieve results
                #

                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']
                state_array_tmp = data_tmp[(group, rid)]['state']

                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[0])

                interval = flux_upper - flux_lower # Added at 200517 to exactly detect edges
                data_sorted.append((flux_upper + interval * 0.0001, thres * 10.0))
                data_sorted.insert(0, (flux_lower - (interval * 0.0001), thres * 10.0))

                #data_sorted.append((flux_upper, thres * 10.0))
                #data_sorted.insert(0, (flux_lower, thres * 10.0))
                below_threshold  = [x for x in data_sorted if x[1] < thres]
                # Detect just before the lower boundary
                flux_previous_left = below_threshold[0][0]
                rss_previous_left = below_threshold[0][1]
                # Detect just before the upper boundary
                flux_previous_right = below_threshold[-1][0]
                rss_previous_right = below_threshold[-1][1]
                # Detect just sfter the lower boundary
                left_group  = sorted([x for x in data_sorted if x[0] < flux_previous_left], key = lambda s: s[1])
                right_group  = sorted([x for x in data_sorted if x[0] > flux_previous_right], key = lambda s: s[1])


                #rss_next_left = left_group[-1][1]
                #flux_next_left = left_group[-1][0]
                rss_next_left = left_group[0][1]
                flux_next_left = left_group[0][0]
                rss_next_right = right_group[0][1]
                flux_next_right = right_group[0][0]


                flux_lower = flux_previous_left - (flux_previous_left - flux_next_left) * (thres - rss_previous_left)/(rss_next_left - rss_previous_left)
                state_lower = "Determined"
                #
                # Check
                #
                if flux_next_left < ci['data'][(group, rid)]['lower_boundary']:
                    state_lower = "Not determined. Rearched to lower boundary"
                if flux_lower < ci['data'][(group, rid)]['lower_boundary']:
                    flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_upper = flux_previous_right + (flux_next_right - flux_previous_right) * (thres - rss_previous_right)/(rss_next_right - rss_previous_right)
                state_upper = "Determined"
                if flux_next_right > ci['data'][(group, rid)]['upper_boundary']:
                    state_upper = "Not determined. Rearched to upper boundary"
                if flux_upper > ci['data'][(group, rid)]['upper_boundary']:
                    flux_upper = ci['data'][(group, rid)]['upper_boundary']

                if callbacklevel >= 3:
                    print(flux_lower, flux_previous_left, rss_previous_left,flux_next_left,rss_next_left)
                    print(flux_upper, flux_previous_right, rss_previous_right,flux_next_right,rss_next_right)
                #
                # Detect lower bounary of confidence interval
                #
                data_tmp[(group, rid)]["log"]["2nd previous_left"] = flux_previous_left
                data_tmp[(group, rid)]["log"]["2nd previous_right"] = flux_previous_right
                data_tmp[(group, rid)]["log"]["2nd next_left"] = flux_next_left
                data_tmp[(group, rid)]["log"]["2nd next_right"] = flux_next_right


                if callbacklevel >= 1:
                    print(rid, 'Grid search maeda method. Lower boundary:', flux_lower, 'Upper boundery' ,flux_upper)

                ci['data'][(group, rid)]['upper_boundary'] = flux_upper
                ci['data'][(group, rid)]['lower_boundary'] = flux_lower
                ci['data'][(group, rid)]['upper_boundary_state'] = state_upper
                ci['data'][(group, rid)]['lower_boundary_state'] = state_lower

                data_tmp[(group, rid)]["log"]["upper_boundary"] = flux_upper
                data_tmp[(group, rid)]["log"]["lower_boundary"] = flux_lower

                #
                # Grids with too large rss are removed
                #
                folds = outputthres #
                flux_values_array = [flux_values_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < (thres * folds) ]
                raw_flux_array = [raw_flux_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < thres * folds ]
                rss_values_array = [rss_values_array_tmp[i] for (i, rsst)  in enumerate(rss_values_array_tmp) if rsst < thres * folds ]
                state_array = [state_array_tmp[i] for (i, rsst)  in enumerate(rss_values_array_tmp) if rsst < thres * folds ]

                ci['data'][(group, rid)]['flux_data'] = flux_values_array
                ci['data'][(group, rid)]['rss_data'] = rss_values_array
                ci['data'][(group, rid)]['raw_flux_data'] = raw_flux_array
                ci['data'][(group, rid)]['state'] = state_array
                #Add upper boundary point
                ci['data'][(group, rid)]['flux_data'].append(flux_upper)
                ci['data'][(group, rid)]['rss_data'].append(thres)
                ci['data'][(group, rid)]['raw_flux_data'].append({})
                ci['data'][(group, rid)]['state'].append(state_lower)
                #Add lower boundary point
                ci['data'][(group, rid)]['flux_data'].append(flux_lower)
                ci['data'][(group, rid)]['rss_data'].append(thres)
                ci['data'][(group, rid)]['raw_flux_data'].append({})
                ci['data'][(group, rid)]['state'].append(state_lower)

                ci['data'][(group, rid)]['log'] = data_tmp[(group, rid)]["log"]
            #
            # Show calc time
            #
            if callbacklevel >= 1:
                endtime=time.perf_counter()
                print("time elapsed ", endtime-starttime,"s")




        #######################################################
        #
        # Grid search joblib
        #
        ########################################################
        if method == "grid":
            #
            # Check Parallel python
            #
            try:
                from joblib import Parallel, delayed
                #job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
            except:
                if callbacklevel >= 1:
                    print("This function requires joblib!")
                return False

            starttime=time.perf_counter()
            step = 0.0
            data_tmp={}
            job_number = self.configuration["grid_search_iterations"]
            initial_search_repeats_in_grid_search = self.configuration["initial_search_repeats_in_grid_search"]
            #
            # confidence interval
            #
            rss_bestfit = self.calc_rss(flux)

            for group, rid in sorted(ci['data'].keys()):
                if ci['data'][(group, rid)]['use'] != 'on': continue
                flux_opt_rid = flux[group][rid]['value']


                data={(group, rid):{"flux_data":[flux_opt_rid],"rss_data":[rss_bestfit],"state":["Best fit"],"log":{},"raw_flux_data":[[flux["reaction"][i]["value"] for i in self.reaction_ids]]}}
                data_tmp.update(data)

                flux_lower =ci['data'][(group, rid)]["lower_boundary"]
                flux_upper =ci['data'][(group, rid)]["upper_boundary"]



                # Grid number is 20
                if (flux_upper-flux_lower)>=100:
                    n=15
                else:
                    n=10

                if callbacklevel >= 2:
                    print("Searching:", rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower," ub ",flux_upper, "n=", n)
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
                    if fixed_flux >= flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux <= flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux,flux_lower,flux_upper, "interation", i)

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallel")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag > 0:
                        counter_of_missed_initial_state = 0
                    else:
                        counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
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
                    if fixed_flux >= flux_upper:
                        fixed_flux = flux_upper
                    if fixed_flux <= flux_lower:
                        fixed_flux = flux_lower
                    #
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux,flux_lower,flux_upper, "interation", i)
                    self.update()

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallel")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", i, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag > 0:
                        counter_of_missed_initial_state = 0
                    else:
                        counter_of_missed_initial_state = counter_of_missed_initial_state + 1
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
                #
                # Pitch jobs to pp
                #
                jobs=[]
                jobsparameters=[]
                #job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in reversed(temp_array_initial_fluxes):
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    self.update()
                    try:
                        rss = self.calc_rss(flux_opt)
                    except:
                        if callbacklevel >= 3:
                            print("Skipped:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)
                        continue

                    if callbacklevel >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)

                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    jobs.append(parameters)
                    jobsparameters.append(((group, rid),fixed_flux, flux_opt))
                    if callbacklevel >= 3:
                        print("New job was added to joblib:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)

                if callbacklevel >= 3:
                    print("Waiting for joblib response")
                #
                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_deep)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) for (configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) in jobs])

                #
                # Return to original position
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()

                if callbacklevel >= 3:
                    print("Waiting for joblib response", rid, "flux_opt ")
                #
                # Retrive results
                #
                #for (group, rid), flux_value, flux_opt, job in jobs:
                for i in range(len(result)):
                    job = result[i]
                    ((group, rid), flux_value, flux_opt) = jobsparameters[i]

                    results = job
                    if results == None:
                        #continue
                        #
                        # Large value is used when falied
                        #
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        #opt_flux = []
                        state = "Failed in optimization"
                    else:
                        state, rss, opt_flux, Rm_ind_sol = results
                        state =  "Finished successfully "
                    #
                    #  Large value is used when falied
                    #
                    if len(opt_flux) == 0:
                        #continue
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        opt_flux = []
                        state = "Failed in optimization"
                    #
                    #  Store data
                    #
                    data_tmp[(group, rid)]["flux_data"].append(flux_value)
                    data_tmp[(group, rid)]["rss_data"].append(rss)
                    data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)
                    data_tmp[(group, rid)]["state"].append(state) #★200517追加　落ちてないのかfialしたのか判断できないので固定値にした

                    if callbacklevel >= 3:
                        print("Finished:", group, rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, "flux_fixed", flux_value, "rss", rss, "state", state)
            #
            # Cals again around thres
            #
            for group, rid in sorted(ci['data'].keys()):
                if ci['data'][(group, rid)]['use'] != 'on': continue


                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']

                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[0])
                #
                # Add right and left end points
                #
                interval = flux_upper - flux_lower # Added at 200517 to exactly detect edges
                data_sorted.append((flux_upper + interval * 0.0001, thres * 11.0))
                data_sorted.insert(0, (flux_lower - (interval * 0.0001), thres * 11.0))
                #
                # points under the threshold
                #
                below_threshold  = [x for x in data_sorted if x[1] < thres]
                # Detect just before the lower boundary
                flux_previous_left = below_threshold[0][0]
                # Detect just before the upper boundary
                flux_previous_right = below_threshold[-1][0]

                # Detect just next to the boundary
                flux_next_left  = max([x[0] for x in data_sorted if x[0] < flux_previous_left])
                flux_next_right  = min([x[0] for x in data_sorted if x[0] > flux_previous_right])
                data_tmp[(group, rid)]["log"]["1st previous_left"] = flux_previous_left
                data_tmp[(group, rid)]["log"]["1st previous_right"] = flux_previous_right
                data_tmp[(group, rid)]["log"]["1st next_left"] = flux_next_left
                data_tmp[(group, rid)]["log"]["1st next_right"] = flux_next_right

                #
                # This part should be reconsiderd 200517
                #
                if (flux_next_right-flux_next_left)>=40:
                    n=20
                elif (flux_next_right-flux_next_left)>=5:
                    n=10
                else:
                    n=5
                #
                #
                #
                if callbacklevel == 2:
                    print("Searching again:", rid, "lb ", flux_next_right," ub ",flux_next_left, "n=", n)
                if callbacklevel >= 3:
                    print("Searching again:", rid,"previous_left ",flux_previous_left,", previous_right ",flux_previous_right, "next_left ",flux_next_left,", flux_next_right ",flux_next_right, "n=", n)
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
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e)
                    self.update()

                    initial_state_flag = 0

                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallel")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag == 0:
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
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
                    # Fix reaction
                    #
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    if callbacklevel >= 3:
                        print("Fixed for initial search:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e)
                    self.update()
                    initial_state_flag = 0
                    for job_n in range(job_number):
                        state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, template = flux, method = "parallel")
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state using template:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            state, flux_opt = self.generate_initial_states(initial_search_repeats_in_grid_search, 1, method = "parallel")
                        #
                        # When initial flux could not be found
                        #
                        if len(flux_opt) == 0:
                            if callbacklevel >= 3:
                                print("Can't find initial state:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux, "interation", e, "job", job_n)
                            continue
                        temp_array_initial_fluxes.append((group, rid, fixed_flux, flux_opt))
                        initial_state_flag = initial_state_flag + 1
                    #
                    # if initial state was found
                    #
                    if initial_state_flag == 0:
                        #
                        # Set large data
                        #
                        data_tmp[(group, rid)]["flux_data"].append(fixed_flux)
                        data_tmp[(group, rid)]["rss_data"].append(thres * 10.0)
                        data_tmp[(group, rid)]["raw_flux_data"].append([])
                        data_tmp[(group, rid)]["state"].append("Failed to find initial state") #
                #
                # Pitch jobs to pp
                #
                jobs=[]
                jobsparameters=[]
                #job_server = pp.Server(ncpus = ncpus, ppservers=ppservers, restart=True, socket_timeout=7200)
                for (group_temp, rid_temp, fixed_flux, flux_opt) in temp_array_initial_fluxes:
                    self.set_constrain(group, rid, "fixed", value = fixed_flux, stdev = 1.0)
                    if callbacklevel >= 3:
                        print("Fixed for fitting:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)
                    self.update()
                    parameters = self.fitting_flux(method = 'deep', flux = flux_opt, output = 'for_parallel')
                    jobs.append(copy.deepcopy(parameters))
                    jobsparameters.append(((group, rid),fixed_flux, flux_opt))
                    if callbacklevel >= 3:
                        print("New job was added to joblib:", rid, "flux_opt ",flux_opt_rid, "value", fixed_flux)
                if callbacklevel >= 3:
                    print("Waiting for joblib response")

                result = Parallel(n_jobs=ncpus)([delayed(optimize.fit_r_mdv_deep)(configuration, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) for (configurationx, experiments, numbers, vectors, matrixinv, calmdv_text, flux_temp) in jobs])
                #
                # Return to original position [need modificaton]
                #
                self.set_constrain(group, rid, temp_type, temp_value, temp_stdev)
                self.update()
                #
                # Retrive results
                #
                #for (group, rid), flux_value, flux_opt, job in jobs:
                for i in range(len(result)):
                    results = result[i]
                    ((group, rid), flux_value, flux_opt) = jobsparameters[i]
                    if results == None:
                        #continue
                        #
                        # Large value is used when falied
                        #
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        #opt_flux = []
                        state = "Failed in optimization"
                    else:
                        state, rss, opt_flux, Rm_ind_sol = results
                        state =  "Finished successfully "
                    #
                    #  Large value is used when falied
                    #
                    if len(opt_flux) == 0:
                        #continue
                        rss = 1000000000000 #★200317追加　落ちてないのかfialしたのか判断できないので固定値にした
                        #rss = thres * 10.0 + abs(flux_value-flux_opt) * 0.00001
                        opt_flux = []
                        state = "Failed in optimization"
                    #
                    #  Store data
                    #
                    data_tmp[(group, rid)]["flux_data"].append(flux_value)
                    data_tmp[(group, rid)]["rss_data"].append(rss)
                    data_tmp[(group, rid)]["raw_flux_data"].append(opt_flux)
                    data_tmp[(group, rid)]["state"].append(state) #★200517追加　落ちてないのかfialしたのか判断できないので固定値にした

                    if callbacklevel >= 3:
                        print("Finished:", group, rid, "flux_opt ",flux_opt_rid, "lb ", flux_lower,", ub ",flux_upper, "flux_fixed", flux_value, "rss", rss, "state", state)

            for (group, rid) in ci['data'].keys():
                if ci['data'][(group, rid)]['use'] != 'on':
                    continue
                #
                # Retirieve results
                #

                flux_upper = ci['data'][(group, rid)]['upper_boundary']
                flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_values_array_tmp = data_tmp[(group, rid)]['flux_data']
                rss_values_array_tmp = data_tmp[(group, rid)]['rss_data']
                raw_flux_array_tmp = data_tmp[(group, rid)]['raw_flux_data']
                state_array_tmp = data_tmp[(group, rid)]['state']

                data_sorted = sorted([(x, rss_values_array_tmp[i]) for i, x in enumerate(flux_values_array_tmp)], key = lambda s: s[0])

                interval = flux_upper - flux_lower # Added at 200517 to exactly detect edges
                data_sorted.append((flux_upper + interval * 0.0001, thres * 10.0))
                data_sorted.insert(0, (flux_lower - (interval * 0.0001), thres * 10.0))

                below_threshold  = [x for x in data_sorted if x[1] < thres]
                # Detect just before the lower boundary
                flux_previous_left = below_threshold[0][0]
                rss_previous_left = below_threshold[0][1]
                # Detect just before the upper boundary
                flux_previous_right = below_threshold[-1][0]
                rss_previous_right = below_threshold[-1][1]
                # Detect just sfter the lower boundary
                left_group  = sorted([x for x in data_sorted if x[0] < flux_previous_left], key = lambda s: s[1])
                right_group  = sorted([x for x in data_sorted if x[0] > flux_previous_right], key = lambda s: s[1])

                rss_next_left = left_group[0][1]
                flux_next_left = left_group[0][0]
                rss_next_right = right_group[0][1]
                flux_next_right = right_group[0][0]


                flux_lower = flux_previous_left - (flux_previous_left - flux_next_left) * (thres - rss_previous_left)/(rss_next_left - rss_previous_left)
                state_lower = "Determined"
                #
                # Check
                #
                if flux_next_left < ci['data'][(group, rid)]['lower_boundary']:
                    state_lower = "Not determined. Rearched to lower boundary"
                if flux_lower < ci['data'][(group, rid)]['lower_boundary']:
                    flux_lower = ci['data'][(group, rid)]['lower_boundary']

                flux_upper = flux_previous_right + (flux_next_right - flux_previous_right) * (thres - rss_previous_right)/(rss_next_right - rss_previous_right)
                state_upper = "Determined"
                if flux_next_right > ci['data'][(group, rid)]['upper_boundary']:
                    state_upper = "Not determined. Rearched to upper boundary"
                if flux_upper > ci['data'][(group, rid)]['upper_boundary']:
                    flux_upper = ci['data'][(group, rid)]['upper_boundary']
                #
                # Detect lower bounary of confidence interval
                #
                data_tmp[(group, rid)]["log"]["2nd previous_left"] = flux_previous_left
                data_tmp[(group, rid)]["log"]["2nd previous_right"] = flux_previous_right
                data_tmp[(group, rid)]["log"]["2nd next_left"] = flux_next_left
                data_tmp[(group, rid)]["log"]["2nd next_right"] = flux_next_right


                if callbacklevel >= 2:
                    print("Finished:", rid, 'Lower boundary:', flux_lower, 'Upper boundery' ,flux_upper)

                ci['data'][(group, rid)]['upper_boundary'] = flux_upper
                ci['data'][(group, rid)]['lower_boundary'] = flux_lower
                ci['data'][(group, rid)]['upper_boundary_state'] = state_upper
                ci['data'][(group, rid)]['lower_boundary_state'] = state_lower

                data_tmp[(group, rid)]["log"]["upper_boundary"] = flux_upper
                data_tmp[(group, rid)]["log"]["lower_boundary"] = flux_lower

                #
                # Grids with too large rss are removed
                #
                folds = outputthres #
                flux_values_array = [flux_values_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < (thres * folds) ]
                raw_flux_array = [raw_flux_array_tmp[i] for (i, rsst) in enumerate(rss_values_array_tmp) if rsst < thres * folds ]
                rss_values_array = [rss_values_array_tmp[i] for (i, rsst)  in enumerate(rss_values_array_tmp) if rsst < thres * folds ]
                state_array = [state_array_tmp[i] for (i, rsst)  in enumerate(rss_values_array_tmp) if rsst < thres * folds ]

                ci['data'][(group, rid)]['flux_data'] = flux_values_array
                ci['data'][(group, rid)]['rss_data'] = rss_values_array
                ci['data'][(group, rid)]['raw_flux_data'] = raw_flux_array
                ci['data'][(group, rid)]['state'] = state_array
                #Add upper boundary point
                ci['data'][(group, rid)]['flux_data'].append(flux_upper)
                ci['data'][(group, rid)]['rss_data'].append(thres)
                ci['data'][(group, rid)]['raw_flux_data'].append({})
                ci['data'][(group, rid)]['state'].append(state_lower)
                #Add lower boundary point
                ci['data'][(group, rid)]['flux_data'].append(flux_lower)
                ci['data'][(group, rid)]['rss_data'].append(thres)
                ci['data'][(group, rid)]['raw_flux_data'].append({})
                ci['data'][(group, rid)]['state'].append(state_lower)

                ci['data'][(group, rid)]['log'] = data_tmp[(group, rid)]["log"]
            #
            # Show calc time
            #
            if callbacklevel >= 2:
                endtime=time.perf_counter()
                print("time elapsed ", endtime-starttime,"s")

        return ci


    def load_states(self, filename, format = 'text'):
        """Load a text/csv file with 'reacton type' information to generate new states dict.

        Args:
            filename (str): filename of flux data with following format::

                type	Id	type	flux_calue	flux_std	lb	ub
                reaction	v1	fixed	100	1	0.001	100
                reaction	v2	free	100	1	0.001	100
                reaction	v3	free	50	1	0.001	100
                metabolites	AcCoA	fixed	1	1	0.001	100
                metabolites	OAC	fixed	1	1	0.001	100
                metabolites	OACs	fixed	1	1	0.001	100
                metabolites	OACx	fixed	1	1	0.001	100
                reversible	FUM	free	1	1	0.001	100

            format: 'csv' CSV or 'text'tab-deliminated text.

        Returns:
            dict: Dictionary of states dict.

        Examples:
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
                    print("This row was ignored", row)
                    continue

                state, rid, type, value, stdev, lb, ub, *over = row
                if state == "reaction":
                    dict = flux_dict
                elif state == "metabolite":
                    dict = conc_dict
                elif state == "reversible":
                    dict = reversible_dict
                else:
                    print("This row was ignored", row)
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
        """Save state dict to a text/csv file with 'type' information.

        Args:
            filename (str): filename of MDV data with the format.
            dict (dict) : dictinary of metabolic state (flux)
            format: 'csv' CSV or 'text'tab-deliminated text.

        Returns:
            Boolean: True/False

        Examples:
            >>> model.save_states(flux_dict, "test.csv", format = 'csv')

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
        """Load MDV data from the text file. This function generate new instance of mfapy.mdv class from a instance of mfapy.model.

        Args:
            filename (str) : filename of MDV data with following format::

                Name	Spectrum	Select	MDV	Std
                Ala57	m0	1	0.566990778	0.000774686
                Ala57	m1	1	0.148623963	0.000774686
                Ala57	m2	1	0.039467636	0.000774686
                Ala57	m3	1	0.244917622	0.000774686

            format (str) : "text" (defalut) or "csv"
            output (str) : "normal" (defalut) or "debug"

        Returns:
            Boolean: True/False

        Examples:
            >>> mdv = load_mdv_data('filename')



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

