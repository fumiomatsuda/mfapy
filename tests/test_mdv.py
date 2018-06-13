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
        target_fragment =  {
        'AKGe': {'atommap': 'AKG_12345','number': 6,'order': 0,'type': 'gcms','use': 'use'},
        'AKGmsms': {'atommap': 'AKG_12345+AKG_1+AKG_2345', 'number': 10,'order': 1, 'type': 'msms', 'use': 'use'}}
        self.mdv = mfapy.mdv.MdvData(target_fragment)
        #print(self.mdv.mdv)
        self.mdv.mdv['AKGe'][0]["ratio"] = 0.5
        self.mdv.mdv['AKGe'][1]["ratio"] = 0.1
        self.mdv.mdv['AKGe'][2]["ratio"] = 0.1
        self.mdv.mdv['AKGe'][3]["ratio"] = 0.1
        self.mdv.mdv['AKGe'][4]["ratio"] = 0.1
        self.mdv.mdv['AKGe'][5]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][0]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][1]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][2]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][3]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][4]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][5]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][6]["ratio"] = 0.5
        self.mdv.mdv['AKGmsms'][7]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][8]["ratio"] = 0.1
        self.mdv.mdv['AKGmsms'][9]["ratio"] = 0.1


    def test_has_data(self):
        boolean = self.mdv.has_data('AKGe', 5)
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
    def test_set_data(self):
        boolean = self.mdv.set_data('AKGe', 5, 0.2, 0.01, "no")
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
        expected = 0.2
        actual = self.mdv.mdv['AKGe'][5]["ratio"]
        self.assertEqual(expected, actual)
    def test_get_data(self):
        boolean = self.mdv.set_data('AKGe', 5, 0.2, 0.01, "no")
        ratio, std, use =  self.mdv.get_data('AKGe', 5)
        expected = 0.2
        actual = self.mdv.mdv['AKGe'][5]["ratio"]
        self.assertEqual(expected, actual)
        expected = 0.01
        actual = self.mdv.mdv['AKGe'][5]["std"]
        self.assertEqual(expected, actual)
        expected = "no"
        actual = self.mdv.mdv['AKGe'][5]["use"]
        self.assertEqual(expected, actual)
    def test_set_mdv_for_ignore(self):
        boolean = self.mdv.set_mdv_for_ignore('AKGe', 5)
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
        expected = "no"
        actual = self.mdv.mdv['AKGe'][5]["use"]
        self.assertEqual(expected, actual)
    def test_set_mdv_for_comparison(self):
        boolean = self.mdv.set_mdv_for_ignore('AKGe', 5)
        boolean = self.mdv.set_mdv_for_comparison('AKGe', 5)
        expected = True
        actual = boolean
        self.assertEqual(expected, actual)
        expected = "use"
        actual = self.mdv.mdv['AKGe'][5]["use"]
        self.assertEqual(expected, actual)
    def test_set_std(self):
        self.mdv.set_data('AKGe', 5, 1.0, 0.01, "no")
        self.mdv.set_std(0.01, method = 'absolute')
        expected = 0.01
        actual = self.mdv.mdv['AKGe'][5]["std"]
        self.assertEqual(expected, actual)
        self.mdv.set_std(0.01)
        expected = 0.01
        actual = self.mdv.mdv['AKGe'][5]["std"]
        self.assertEqual(expected, actual)
    def test_set_mdvs_for_comparison(self):
        self.mdv.set_data('AKGe', 5, 0.005, 0.01, "use")
        self.mdv.set_mdvs_for_comparison(0.01)
        expected = "no"
        actual = self.mdv.mdv['AKGe'][5]["use"]
        self.assertEqual(expected, actual)
        expected = "use"
        actual = self.mdv.mdv['AKGe'][4]["use"]
        self.assertEqual(expected, actual)
    def test_check(self):
        self.mdv.set_data('AKGe', 5, 1.01, 0.01, "use")
        expected = False
        actual = self.mdv.check(output = "normal")
        self.assertEqual(expected, actual)
        self.mdv.set_data('AKGe', 5, 0.00, 0.01, "use")
        expected = False
        actual = self.mdv.check(output = "normal")
        self.assertEqual(expected, actual)
        self.mdv.set_data('AKGe', 5, 0.02, 0.01, "use")
        expected = True
        actual = self.mdv.check(output = "normal")
        self.assertEqual(expected, actual)

    def test_get_fragments_for_mdv_calculation(self):
        list = self.mdv.get_fragments_for_mdv_calculation()
        expected = ['AKGe', 'AKGmsms']
        actual = list
        self.assertEqual(expected, actual)
    def test_set_observed_fragments(self):
        self.mdv.set_observed_fragments(['AKGe'])
        list = self.mdv.get_fragments_for_mdv_calculation()
        expected = ['AKGe']
        actual = self.mdv.observed_fragments
        self.assertEqual(expected, actual)
        expected = 'no'
        actual = self.mdv.mdv['AKGmsms'][4]["use"]
        self.assertEqual(expected, actual)
    def test_add_gaussian_noise(self):
        self.mdv.add_gaussian_noise(0.01, 1, method = 'absolute')

        expected = 0.0
        actual = self.mdv.mdv['AKGmsms'][4]["std"]
        self.assertEqual(expected, actual)
    def test_generate_observed_mdv(self):
        list, ratio, *remains = self.mdv.generate_observed_mdv()
        expected = 0.5
        actual = ratio[0]
        self.assertEqual(expected, actual)
    def test_get_number_of_measurement(self):
        expected = 14
        actual = self.mdv.get_number_of_measurement()
        self.assertEqual(expected, actual)
        self.mdv.set_mdv_for_ignore('AKGe', 5)
        expected = 13
        actual = self.mdv.get_number_of_measurement()
        self.assertEqual(expected, actual)

if __name__ == '__main__':
    unittest.main()