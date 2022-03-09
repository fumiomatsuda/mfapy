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
        carbon_sources =  {'AcCoA': {'IDV': [0.5, 0.0, 0.25, 0.25], 'size': 2},
         'Asp': {'IDV': [1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'size': 4},
         'OACs': {'IDV': [1.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0],
                  'size': 4}}
        carbon_source_emu = {'Asp_1:2:3:4': {'metabolite_name': 'Asp', 'position_list': ['1', '2', '3', '4']},'Asp_1:2:3': {'metabolite_name': 'Asp',\
        'position_list': ['1', '2', '3']},'Asp_1:2': {'metabolite_name': 'Asp', 'position_list': ['1', '2']},'Asp_1': {'metabolite_name': 'Asp', 'position_list': ['1']},\
        'AcCoA_1': {'metabolite_name': 'AcCoA', 'position_list': ['1']},'AcCoA_1:2': {'metabolite_name': 'AcCoA', 'position_list': ['1','2']}}
        self.cs = mfapy.carbonsource.CarbonSource(carbon_sources, carbon_source_emu)


    def test_set_all_isotopomers(self):
        boolean = self.cs.set_all_isotopomers('AcCoA', [0.3, 0.3, 0.3, 0.1])
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
        expected = [0.6, 0.4]
        actual = self.cs.mdv_carbon_sources['AcCoA_1']
        #print(actual)
        self.assertEqual(expected, actual)
    def test_set_each_isotopomer(self):
        boolean = self.cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = 'yes')
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
        actual = self.cs.mdv_carbon_sources['Asp_1'][0]
        expected = 0.9893000129677051
        self.assertEqual(expected, actual)
        actual = self.cs.mdv_carbon_sources['Asp_1:2'][0]
        expected = 0.9787145028289506
        self.assertEqual(expected, actual)
        actual = self.cs.mdv_carbon_sources['Asp_1:2:3'][0]
        expected = 0.9682422576486809
        self.assertEqual(expected, actual)
        actual = self.cs.mdv_carbon_sources['Asp_1:2:3:4'][0]
        expected = 0.95788206549184
        self.assertEqual(expected, actual)
if __name__ == '__main__':
    unittest.main()