# -*- coding: utf-8 -*-

#from .context import mfapy
import mfapy
import unittest


class Testmfapy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # procedures before tests are started. This code block is executed only once
        pass

    def setUp(self):
        self.state_filename = "./tests/test_metabolicmodel_state.csv"

        reactions = {'e1': {'atommap': 'nd',
                'lb': 0.1,
                'order': 9,
                'reaction': 'nd',
                'stoichiometry': 'Glu-->Gluex',
                'ub': 300.0,
                'use': 'use'},
         'v1': {'atommap': 'AB+CDEF-->FEDBAC',
                'lb': 0.1,
                'order': 0,
                'reaction': 'AcCoA+OAC-->Cit',
                'stoichiometry': 'AcCoA+OAC-->Cit',
                'ub': 300.0,
                'use': 'use'},
         'v10': {'atommap': 'ABCD-->ABCD',
                 'lb': 0.1,
                 'order': 10,
                 'reaction': 'OAC-->OACx',
                 'stoichiometry': 'OAC-->OACx',
                 'ub': 300.0,
                 'use': 'use'},
         'v11': {'atommap': 'ABCD-->ABCD',
                 'lb': 0.001,
                 'order': 11,
                 'reaction': 'OACs-->OACx',
                 'stoichiometry': 'OACs-->OACx',
                 'ub': 300.0,
                 'use': 'use'},
         'v12': {'atommap': 'nd',
                 'lb': 0.1,
                 'order': 12,
                 'reaction': 'nd',
                 'stoichiometry': 'OACx-->OACe',
                 'ub': 300.0,
                 'use': 'use'},
         'v13': {'atommap': 'ABCD+EFGH-->EFAB',
                 'lb': 0.1,
                 'order': 13,
                 'reaction': 'OAC+OAC-->Fumx',
                 'stoichiometry': 'OAC+OAC-->Fumx',
                 'ub': 300.0,
                 'use': 'use'},
         'v14': {'atommap': 'nd',
                 'lb': 0.1,
                 'order': 14,
                 'reaction': 'nd',
                 'stoichiometry': 'OA-->Fume',
                 'ub': 300.0,
                 'use': 'use'},
         'v15': {'atommap': 'ABCD+EFGH+IJKL-->CDHL',
                 'lb': 0.1,
                 'order': 15,
                 'reaction': 'OAC+Fum+Suc-->Tx',
                 'stoichiometry': 'OAC+Fum+Suc-->Tx',
                 'ub': 300.0,
                 'use': 'use'},
         'v16': {'atommap': 'nd',
                 'lb': 0.1,
                 'order': 16,
                 'reaction': 'nd',
                 'stoichiometry': 'Tx-->Te',
                 'ub': 300.0,
                 'use': 'use'},
         'v17': {'atommap': 'ABCD-->AB+CD',
                 'lb': 0.1,
                 'order': 17,
                 'reaction': 'Fumx-->OA+OA',
                 'stoichiometry': 'Fumx-->{2}OA',
                 'ub': 300.0,
                 'use': 'use'},
         'v2': {'atommap': 'ABCDEF-->ABCDE+F',
                'lb': 0.1,
                'order': 1,
                'reaction': 'Cit-->AKG+CO2ex',
                'stoichiometry': 'Cit-->AKG+CO2ex',
                'ub': 300.0,
                'use': 'use'},
         'v3': {'atommap': 'ABCDE-->ABCDE',
                'lb': 0.1,
                'order': 2,
                'reaction': 'AKG-->Glu',
                'stoichiometry': 'AKG-->Glu',
                'ub': 300.0,
                'use': 'use'},
         'v4': {'atommap': 'ABCDE-->BCDE+A',
                'lb': 0.1,
                'order': 3,
                'reaction': 'AKG-->Suc+CO2ex',
                'stoichiometry': 'AKG-->Suc+CO2ex',
                'ub': 300.0,
                'use': 'use'},
         'v5': {'atommap': 'ABCD-->ABCD',
                'lb': 0.1,
                'order': 4,
                'reaction': 'Suc-->Fum',
                'stoichiometry': 'Suc-->Fum',
                'ub': 300.0,
                'use': 'use'},
         'v6': {'atommap': 'ABCD-->ABCD',
                'lb': 0.1,
                'order': 5,
                'reaction': 'Fum-->OAC',
                'stoichiometry': 'Fum-->OAC',
                'ub': 300.0,
                'use': 'use'},
         'v7': {'atommap': 'ABCD-->ABCD',
                'lb': 0.1,
                'order': 6,
                'reaction': 'OAC-->Fum',
                'stoichiometry': 'OAC-->Fum',
                'ub': 300.0,
                'use': 'use'},
         'v8': {'atommap': 'ABCD-->ABCD',
                'lb': 0.1,
                'order': 7,
                'reaction': 'Asp-->OAC',
                'stoichiometry': 'Asp-->OAC',
                'ub': 300.0,
                'use': 'use'},
         'v9': {'atommap': 'nd',
                'lb': 0.1,
                'order': 8,
                'reaction': 'nd',
                'stoichiometry': 'Suc-->Sucex',
                'ub': 300.0,
                'use': 'use'}}
        reversible =  {'FUM': {'forward': 'v6', 'order': 0, 'reverse': 'v7', 'type': 'free'},
         'Input': {'forward': 'v9+e1',
                   'order': 1,
                   'reverse': 'nothing',
                   'type': 'free'}}
        metabolites = {'AKG': {'C_number': 5,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 12,
             'symmetry': 'no'},
        'AcCoA': {'C_number': 2,
               'carbonsource': 'carbonsource',
               'excreted': 'no',
               'order': 1,
               'symmetry': 'no'},
        'Asp': {'C_number': 4,
             'carbonsource': 'carbonsource',
             'excreted': 'no',
             'order': 18,
             'symmetry': 'no'},
        'CO2ex': {'C_number': 1,
               'carbonsource': 'no',
               'excreted': 'excreted',
               'order': 0,
               'symmetry': 'no'},
        'Cit': {'C_number': 6,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 11,
             'symmetry': 'no'},
        'Fum': {'C_number': 4,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 14,
             'symmetry': 'symmetry'},
        'Fume': {'C_number': 4,
              'carbonsource': 'no',
              'excreted': 'excreted',
              'order': 8,
              'symmetry': 'no'},
        'Fumx': {'C_number': 4,
              'carbonsource': 'no',
              'excreted': 'no',
              'order': 6,
              'symmetry': 'no'},
        'Glu': {'C_number': 5,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 15,
             'symmetry': 'no'},
        'Gluex': {'C_number': 5,
               'carbonsource': 'no',
               'excreted': 'excreted',
               'order': 16,
               'symmetry': 'no'},
        'OA': {'C_number': 2,
            'carbonsource': 'no',
            'excreted': 'no',
            'order': 7,
            'symmetry': 'no'},
        'OAC': {'C_number': 4,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 2,
             'symmetry': 'no'},
        'OACe': {'C_number': 4,
              'carbonsource': 'no',
              'excreted': 'excreted',
              'order': 5,
              'symmetry': 'no'},
        'OACs': {'C_number': 4,
              'carbonsource': 'carbonsource',
              'excreted': 'no',
              'order': 3,
              'symmetry': 'no'},
        'OACx': {'C_number': 4,
              'carbonsource': 'no',
              'excreted': 'no',
              'order': 4,
              'symmetry': 'no'},
        'Suc': {'C_number': 4,
             'carbonsource': 'no',
             'excreted': 'no',
             'order': 13,
             'symmetry': 'symmetry'},
        'Sucex': {'C_number': 4,
               'carbonsource': 'no',
               'excreted': 'excreted',
               'order': 17,
               'symmetry': 'no'},
        'Te': {'C_number': 4,
            'carbonsource': 'no',
            'excreted': 'excreted',
            'order': 10,
            'symmetry': 'no'},
        'Tx': {'C_number': 4,
            'carbonsource': 'no',
            'excreted': 'no',
            'order': 9,
            'symmetry': 'no'}}
        target_fragments = {'AKGc': {'atommap': 'AKG_1:2+AKG_3:4:5',
                  'order': 11,
                  'type': 'gcms',
                  'use': 'no'},
         'AKGe': {'atommap': 'AKG_1:2:3:4:5', 'order': 1, 'type': 'gcms', 'use': 'use'},
         'AKGms': {'atommap': 'AKG_1:2:3:4:5', 'order': 10, 'type': 'gcms', 'use': 'no'},
         'AKGmsms': {'atommap': 'AKG_1:2:3:4:5+AKG_1+AKG_2:3:4:5',
                     'order': 9,
                     'type': 'msms',
                     'use': 'no'},
         'Glue': {'atommap': 'Glu_1:2:3:4:5', 'order': 0, 'type': 'gcms', 'use': 'use'},
         'Gluf': {'atommap': 'Glu_1:2+Glu_3:4:5', 'order': 8, 'type': 'gcms', 'use': 'no'},
         'OACc': {'atommap': 'OAC_1:2+OAC_1:2', 'order': 5, 'type': 'gcms', 'use': 'use'},
         'OACi': {'atommap': 'OAC_1:2:3:4', 'order': 3, 'type': 'gcms', 'use': 'use'},
         'OACo': {'atommap': 'OACx_1:2:3:4', 'order': 2, 'type': 'gcms', 'use': 'use'},
         'OACt': {'atommap': 'Fumx_1:2:3:4', 'order': 4, 'type': 'gcms', 'use': 'use'},
         'Txc': {'atommap': 'OAC_3:4+Fum_4+Suc_4',
                 'order': 7,
                 'type': 'gcms',
                 'use': 'use'},
         'Txt': {'atommap': 'Tx_1:2:3:4', 'order': 6, 'type': 'gcms', 'use': 'use'}}

        self.model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
        # procedures before every tests are started. This code block is executed every time
        pass

    def test_load_metabolic_model_construction(self):
        expected = 'v1'
        actual = self.model.reaction_ids[0]
        self.assertEqual(expected, actual)


    def test_set_constrain(self):
        self.model.set_constrain('reaction', 'v1', "fixed", 100.0,1)
        expected = 'fixed'
        actual = self.model.reactions["v1"]["type"]
        self.assertEqual(expected, actual)

    def test_set_boundary(self):
        self.model.set_boundary('reaction', 'v12', 0.001, 1.0)
        expected = 1.0
        actual = self.model.reactions["v12"]["ub"]
        self.assertEqual(expected, actual)

    def test_set_configuration(self):
        self.model.set_configuration(iteration_max = 10001)
        expected = 10001
        actual = self.model.configuration["iteration_max"]
        self.assertEqual(expected, actual)

    def test_update(self):
        self.model.set_constrain('reaction', 'v1', "fixed", 100.0,1)
        expected = True
        actual = self.model.update()
        self.assertEqual(expected, actual)

    def test_load_states(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        expected = 100.0
        actual = state_dic["reversible"]["FUM"]['ub']
        self.assertEqual(expected, actual)
        expected = 500.0
        actual = state_dic["reaction"]["v1"]['ub']
        self.assertEqual(expected, actual)
        expected = 100.0
        actual = state_dic['metabolite']['AKG']['ub']
        self.assertEqual(expected, actual)
        flux_list = [state_dic["reaction"][id]['value'] for id in self.model.reaction_ids]
        expected = [100.0, 100.0, 50.0, 50.0, 50.0, 125.0, 75.0, 50.0, 0.0, 50.0, 1.0, 0.0, 100.0, 1.0, 1.0, 1.0, 1.0, 0.0]
        actual = flux_list
        self.assertEqual(expected, actual)

    def test_set_constraints_from_state_dict(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        expected = 100.0
        actual = self.model.reversible["FUM"]['ub']
        self.assertEqual(expected, actual)
        expected = 500.0
        actual = self.model.reactions["v1"]['ub']
        self.assertEqual(expected, actual)
        expected = 100.0
        actual = self.model.metabolites['AKG']['ub']
        self.assertEqual(expected, actual)

    def test_generate_flux_distribution(self):
        flux_opt, state = self.model.generate_flux_distribution()
        expected = "Determined"
        actual = state
        self.assertEqual(expected, actual)

    def test_matrix_inv(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.update()
        tmp_r= [state_dic[group][id]["value"] for (group, id) in self.model.vector["ids"]]
        ind_r = [state_dic[group][id]["value"] for (group, id) in self.model.vector["independent_flux"]]
        import numpy

        Rm_initial = list(self.model.vector["Rm_initial"])
        Rm_initial[self.model.numbers["independent_start"]:self.model.numbers["independent_end"]] = ind_r
        flux_list = [state_dic["reaction"][id]['value'] for id in self.model.reaction_ids]
        expected =  [100.0, 100.0, 50.0, 50.0, 50.0, 125.0, 75.0, 50.0, 0.0, 50.0, 100.0, 0.0, 100.0, 0.5, 1.0, 1.0, 1.0, 0.5, 10.0, 10.0, 10.0, 10.0, 10.0, 100.0, 100.0, 100.0, 100.0, 100.0, 50.0, 50.0]
        actual = list(numpy.dot(self.model.matrixinv, Rm_initial))
        self.assertEqual(expected, actual)


    def test_calmdv(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        cs.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.25, 0.25])
        cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = 'no')
        flux_list = [state_dic["reaction"][id]['value'] for id in self.model.reaction_ids]
        glu, glut = self.model.func["calmdv"](flux_list, ['Glue'], cs.generate_dict())
        expected =  [0.34635416666666663, 0.26953124999999994, 0.27083333333333326, 0.08072916666666664, 0.028645833333333325, 0.0039062499999999987]
        actual = list(glu)
        self.assertAlmostEqual(expected[0], actual[0])
        expected = 1.0
        actual = sum(list(glu))
        self.assertAlmostEqual(expected, actual)
    """
    def test_save_states(self):
        flux_opt, state = self.model.generate_flux_distribution()
        actual = self.model.save_states(flux_opt, "Example_1_test.csv", format = 'csv')
        expected = True
        self.assertEqual(expected, actual)
    """

    def test_set_experiment(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)

        self.model.set_experiment('ex1', mdv1, cs)
        actual = self.model.experiments["ex1"]["mode"]
        expected = "ST"
        self.assertEqual(expected, actual)

    def test_calc_rss(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)
        rss = self.model.calc_rss(flux_opt)
        actual = rss
        expected = 0.0
        self.assertAlmostEqual(expected, actual)

    def test_generate_initial_states(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.set_constrain('reaction','v8','free')
        self.model.set_constrain('reversible','FUM','free')
        self.model.update()
        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)

        state, flux_initial1 = self.model.generate_initial_states(50, 1)
        actual = state
        expected = 50
        self.assertEqual(expected, actual)

    def test_fitting_flux(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.set_constrain('reaction','v8','free')
        self.model.set_constrain('reversible','FUM','free')
        self.model.update()
        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)

        state, flux_initial1 = self.model.generate_initial_states(50, 1)
        # call scipy
        state, RSS_bestfit, flux_opt_slsqp = self.model.fitting_flux(method = 'SLSQP', flux = flux_initial1)
        actual = state
        expected = "Optimization terminated successfully."
        self.assertEqual(expected, actual)
        # call nlopt
        state, RSS_bestfit, flux_opt_slsqp = self.model.fitting_flux(method = 'LN_PRAXIS', flux = flux_initial1)
        actual = state["text"]
        expected = 'Function was called'
        self.assertEqual(expected, actual)
    """
    def test_fitting_flux_parallel(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.set_constrain('reaction', 'v9', "fitting", 0.0,1)
        self.model.set_constrain('reaction','v8','free')
        self.model.set_constrain('reversible','FUM','free')
        self.model.update()
        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)

        state, flux_initial1 = self.model.generate_initial_states(50, 2)
        # call scipy
        state, RSS_bestfit, flux_opt_slsqp = self.model.fitting_flux(method = 'SLSQP', flux = flux_initial1)
        actual = state[0]
        expected = "Optimization terminated successfully."
        self.assertEqual(expected, actual)
        # call nlopt
        state, RSS_bestfit, flux_opt_slsqp = self.model.fitting_flux(method = 'LN_PRAXIS', flux = flux_initial1)
        actual = state[0]["text"]
        expected = 'Function was called'
        self.assertEqual(expected, actual)
    """


    def test_get_thres_confidence_interval(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)
        thres, number_of_measurements, degree_of_freedom  = self.model.get_thres_confidence_interval(flux_opt, alpha = 0.05)
        actual = degree_of_freedom
        expected = 2
        self.assertEqual(expected, actual)
        actual = number_of_measurements
        expected = 34
        self.assertEqual(expected, actual)
        actual = thres
        expected = 0.0
        self.assertAlmostEqual(expected, actual)

    def test_get_number_of_independent_measurements(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)
        number_of_independent_measurements  = self.model.get_number_of_independent_measurements()
        actual = number_of_independent_measurements
        expected = 34
        self.assertEqual(expected, actual)
    def test_get_degree_of_freedom(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)
        result  = self.model.get_degree_of_freedom()
        #print(result)
        actual = result
        expected = 2
        self.assertEqual(expected, actual)
    def test_goodness_of_fit(self):
        state_dic = self.model.load_states(self.state_filename, format = 'csv')
        self.model.set_constraints_from_state_dict(state_dic)
        self.model.update()

        cs = self.model.generate_carbon_source_templete()
        flux_opt, state = self.model.generate_flux_distribution()
        mdv1 = self.model.generate_mdv(flux_opt, cs)
        self.model.set_experiment('ex1', mdv1, cs)
        pvalue, rss_thres= self.model.goodness_of_fit(flux_opt, alpha =0.05)
        #print(rss_thres)
        actual = rss_thres
        expected = 46.19425952027847
        self.assertAlmostEqual(expected, actual)

if __name__ == '__main__':
    unittest.main()