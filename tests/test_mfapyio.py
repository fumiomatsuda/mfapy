# -*- coding: utf-8 -*-

#from .context import mfapy
import mfapy
import unittest


class Testmfapy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
        # procedures before tests are started. This code block is executed only once



    def setUp(self):
        # procedures before every tests are started. This code block is executed every time
        self.filename = "./tests/test_mfapyio_model.txt"




    def test_load_metabolic_model_reactions(self):

        reactions = mfapy.mfapyio.load_metabolic_model_reactions(self.filename, format="text", output = "normal")
        expected = float(300.0)
        actual = reactions["v17"]['ub']
        self.assertEqual(expected, actual)
    	#{'atommap': 'ABCD-->AB+CD',
    	# 'lb': 0.1,
    	# 'order': 17,
    	# 'reaction': 'Fumx-->OA+OA',
    	# 'stoichiometry': 'Fumx-->{2}OA',
        # 'ub': 300.0,
        # 'use': 'use'}

    def test_load_metabolic_model_reversibles(self):
        reversible_reactions = mfapy.mfapyio.load_metabolic_model_reversibles(self.filename, format="text", output = "normal")
        expected = 'free'
        actual = reversible_reactions["Input"]['type']
        self.assertEqual(expected, actual)
    	#>>> reversible["Input"]
        #{'forward': 'v9+e1', 'order': 1, 'reverse': 'nothing', 'type': 'free'}

    def test_load_metabolic_model_metabolites(self):
        metabolites = mfapy.mfapyio.load_metabolic_model_metabolites(self.filename, format="text", output = "normal")
        expected = 'no'
        actual = metabolites["Asp"]['excreted']
        self.assertEqual(expected, actual)
    	#>>> metabolites["Asp"]
        #{'C_number': 4,'carbonsource': 'carbonsource', 'excreted': 'no', 'order': 18, 'symmetry': 'no'}

    def test_load_metabolic_model_fragments(self):
        target_fragments = mfapy.mfapyio.load_metabolic_model_fragments(self.filename, format="text", output = "normal")
        expected = 'no'
        actual = target_fragments["AKGc"]['use']
        self.assertEqual(expected, actual)
    	#>>> target_fragments["AKGc"]
        #{'atommap': 'AKG_12+AKG_345', 'order': 11, 'type': 'gcms', 'use': 'no'}

    def test_load_metabolic_model(self):
        reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model(self.filename, format = "text")


        expected = 'no'
        actual = target_fragments["AKGc"]['use']
        self.assertEqual(expected, actual)

        expected = 'no'
        actual = metabolites["Asp"]['excreted']
        self.assertEqual(expected, actual)

        expected = 'free'
        actual = reversible["Input"]['type']
        self.assertEqual(expected, actual)

        expected = float(300.0)
        actual = reactions["v17"]['ub']
        self.assertEqual(expected, actual)



if __name__ == '__main__':
    unittest.main()