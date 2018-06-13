#-------------------------------------------------------------------------------
# Name:        mdv.py
# Purpose:     MdvData class in mfapy
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------

import numpy as numpy
from scipy.stats import chi2, f
import itertools
import re
import math


class MdvData:
    """
    Class of mass distribution vector (MDV) information. A instance of the MdvData class could be generate from
    an instance of MetabolicModel class using model.generate_mdv_templete()
    """
    def __init__(self, target_fragments):
        """
        Constructer of MdvData instance from target_fragments information of the metabolic model

        Parameters
        ----------
        target_fragments : dictionary of target fragments.

        Examples
        --------
        >>> mdv = MdvData(target_fragments)
        Mdv data can be accessed by:
        mdv.mdv['fragment_name']['0'][ratio]

        See Also
        --------

        """
        self.mdv = {}
        self.mode = "single"
        #
        # List of fragments now considered
        #
        self.observed_fragments_set = set()
        #
        # List of fragments described in model definition file
        #
        self.fragments_for_mdv_calculation_set = set()
        for fragment in target_fragments:
            #
            #ignore unused data
            #
            if target_fragments[fragment]['use'] != 'use':
                continue
            #
            #discard space for data clean up
            #
            self.observed_fragments_set.update([fragment])
            self.fragments_for_mdv_calculation_set.update([fragment])
            self.number_of_replicate = 1
            self.mdv[fragment] = {}
            for i in range(target_fragments[fragment]['number']):
                id = fragment+"_"+str(i)
                self.mdv[fragment][i] = {"id":id, "ratio":float(0.0), "std":float(0.0), "use": "use", "data":[]}
        self.fragments_for_mdv_calculation = sorted(list(self.fragments_for_mdv_calculation_set))
        self.observed_fragments = sorted(list(self.observed_fragments_set))

    def has_data(self, fragment, number):
        """
        Checker whether the MdvData instance has data of [number] th isotope of [fragment] mass spectra.

        Parameters
        ----------
        fragment: name of target_fragment
        number: number of isotope

        Reterns
        ----------
        Boolean: True\False

        Examples
        --------
        >>> if mdv.has_data('Leu85', 2): brabra

        See Also
        --------

        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                return True
        return False

    def set_data(self, fragment, number, ratio, std, use):
        """
        Setter of single MDV data

        Parameters
        ----------
        fragment: name of target_fragment
        number: Number of isotope
        ratio: Relative abundance of mass spectra data
        std: Standard deviation of measured MDV data
        use: Data to be used for MDV comparison. 'use' or 'no'

        Reterns
        ----------
        Boolean: True\False

        Examples
        --------
        >>> mdv.set_data('Leu85', 2, 0.763, 0.0054, 'use')

        See Also
        --------

        """
        if self.has_data(fragment, number):
            self.mdv[fragment][number]['ratio'] = ratio
            self.mdv[fragment][number]['std'] = std
            self.mdv[fragment][number]['use'] = use
            return True
        return False

    def get_data(self, fragment, number):
        """
        Getter of single MDV data

        Parameters
        ----------
        fragment: name of target_fragment
        number: Number of isotope

        Reterns
        ----------
        ratio: Relative abundance of mass spectra data
        std: Standard deviation of measured MDV data
        use: Data to be used for MDV comparison. 'use' or 'no'

        Examples
        --------
        >>> ratio, std, use = mdv.get_data('Leu85', 2)

        See Also
        --------

        """
        if self.has_data(fragment, number):
            return self.mdv[fragment][number]['ratio'], self.mdv[fragment][number]['std'], self.mdv[fragment][number]['use']
        return False, False, False


    def set_observed_fragments(self, fragments):
        """
        Setter of observed fragment. When new MdvData instance is generated 'observed_fragment' == 'target_fragment'.

        Parameters
        ----------
        fragments: List of fragment names calculated in calmdv fucntion

        Examples
        --------
        >>> mdv.set_observed_fragments(['Phe_85', 'Ala_57'])

        See Also
        --------

        """
        # Empty set
        self.observed_fragments_set = set()
        # Add new fragments
        for fragment in fragments:
            # Check fragments
            if fragment in self.fragments_for_mdv_calculation:
                self.observed_fragments_set.update([fragment])
        self.observed_fragments = sorted(self.observed_fragments_set)

        #
        # Set "no" to other fragments
        #

        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                if fragment in self.observed_fragments:
                    self.set_mdv_for_comparison(fragment, number)
                else:
                    self.set_mdv_for_ignore(fragment, number)
        return True


    def get_fragments_for_mdv_calculation(self):
        """
        Getter of fragment list whose MDVs are calculated in calmdv fucntion

        Reterns
        ----------
        fragments: List of fragment names to be used for MDV comparison


        Examples
        --------
        >>> list_of_fragment = mdv.get_fragments_for_mdv_calculation()

        See Also
        --------

        """
        return list(self.fragments_for_mdv_calculation)

    def set_mdv_for_comparison(self, fragment, number):
        """
        Setter of single mass data to be used for MDV comparison

        Reterns
        ----------
        fragments and number: fragment names and its mass data  to be used for MDV comparison

        Examples
        --------
        >>> mdv.set_mdv_for_comparison(fragment, number)

        See Also
        --------

        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                self.mdv[fragment][number]['use'] = 'use'
                return True
        return False

    def set_mdv_for_ignore(self, fragment, number):
        """
        Setter of single mass data to be ignored for MDV comparison

        Reterns
        ----------
        fragments and number: fragment names and its mass data  to be ignored for MDV comparison

        Examples
        --------
        >>> mdv.set_mdv_for_ignored(fragment, number)

        See Also
        --------
        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                self.mdv[fragment][number]['use'] = 'no'
                return True
        return False

    def set_std(self, value, method = 'absolute'):
        """
        Set standard deviation levels from mass sepctral intensity

        Parameters
        ----------
        method: Method to calculate stdev levels:
            'relative' (detault): Levels are set by stdev = [value] * signal intensity
            'abusolute': Levels are set by stdev = [value]

        Examples
        --------
        >>> mdv.set_std(0.05, method = 'absolute')

        See Also
        --------
        """
        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                if method == 'absolute':
                    self.mdv[fragment][number]['std'] = value * 1.0
                else:
                    self.mdv[fragment][number]['std'] = self.mdv[fragment][number]['ratio'] * value

    def check(self, output = "debug"):
        """
        Function to check missing values in the MDV data before MDV comparison

        Parameters
        ----------
        output: Show details when "debug" mode

        Returns
        --------
        boolean: True/False

        Examples
        --------
        >>> mdv.check()

        See Also
        --------
        """
        counter = 0;
        for fragment in self.mdv:
            if fragment in self.observed_fragments:
                for number in self.mdv[fragment]:
                    if self.mdv[fragment][number]['use'] == 'use':
                        if self.mdv[fragment][number]['ratio'] > 1.0:
                            counter = counter + 1
                            if output == "debug":
                                print(fragment, number, self.mdv[fragment][number]['ratio'],"> 1.0")
                        if self.mdv[fragment][number]['ratio'] <= 0.0:
                            counter = counter + 1
                            if output == "debug":
                                print(fragment, number, self.mdv[fragment][number]['ratio'],"<= 0.0")
            else:
                if output == "debug":
                    print('The MDV object does not include', fragment)
                counter = counter + 1
        if counter == 0:
            return True
        return False

    def set_mdvs_for_comparison(self, threshold_ratio, threshold_std = 1.0):
        """
        Setters of multiple data to be used for MDV comparison

        Parameters
        --------
        threshold_ratios: Mass spectral signals with normalized intensities above threshold_ratios are used to be used for MDV comparison
        threshold_std: Threshold for std data

        Returns
        --------
        boolean: True/False

        Examples
        --------
        >>> mdv.set_mdvs_for_comparison(0.05, 0.005)

        See Also
        --------
        """
        for fragment in self.fragments_for_mdv_calculation:
            for number in self.mdv[fragment]:
                self.set_mdv_for_ignore(fragment, number)
                if fragment in self.observed_fragments:
                    ratio, std, use = self.get_data(fragment, number)
                    if (ratio > threshold_ratio) and (std < threshold_std):
                        self.set_mdv_for_comparison(fragment, number)

    def add_gaussian_noise(self, stdev, iteration, method = 'absolute'):
        """
        Add gausian noise to MDV data

        Parameters
        --------
        ratio: Intensity of noise;
        iteration: Number of sampleing
        method:
            relative: intentsity = (randn() * stdev + 1) * original_intensity
            absolute: intentsity = randn() * stdev) + original_intensity

        Returns
        --------
        boolean: True/False

        Examples
        --------
        >>> mdv.add_gaussian_noise(0.05, 5)

        See Also
        --------
        """

        self.number_of_replicate = iteration
        for fragment in self.fragments_for_mdv_calculation:
            number_of_mass_data = max(self.mdv[fragment].keys()) + 1
            noise = numpy.zeros((number_of_mass_data, iteration))
            for i in range(iteration):
                for number in range(number_of_mass_data):
                    if method == 'relative':
                        noise[number, i] = (numpy.random.randn() * stdev + 1) * self.mdv[fragment][number]['ratio']
                    else:
                        noise[number, i] = (numpy.random.randn() * stdev) + self.mdv[fragment][number]['ratio']
                    if noise[number, i] < 0.0:
                        noise[number, i] = 0.0
                #各フラグメント毎に総和を１にする。
                sumvalue = sum(noise[:,i])
                noise[:,i] = noise[:,i] / sumvalue
            for number in range(number_of_mass_data):
                #各フラグメント毎に総和を１にする。
                self.mdv[fragment][number]['ratio']= sum(noise[number,:])/iteration
                self.mdv[fragment][number]['std'] = numpy.std(noise[number,:])
                self.mdv[fragment][number]['data'] = numpy.array(noise[number,:])

    def generate_observed_mdv(self):
        """
        Generator of MDV related inforamtion in list format.
        This function is used in model.set_experiment()

        Returns
        --------
        id_array: List of ids
        ratio_array: List of ratio data
        std_array:List of std data
        use_array: List of use data
        observed_fragments: List of observed fragment
        data:

        Examples
        --------
        >>> id_array, ratio_array, std_array, use_array, observed_fragments, data = mdv.generate_observed_mdv()

        See Also
        --------
        """
        id_array = []
        ratio_array = []
        std_array = []
        use_array = []
        data = numpy.zeros(self.number_of_replicate)
        for fragment in sorted(self.observed_fragments):
            for number in sorted(self.mdv[fragment].keys()):
                ratio, std, use = self.get_data(fragment, number)
                id_array.append(self.mdv[fragment][number]['id'])#これださい
                ratio_array.append(ratio)
                std_array.append(std)
                if self.number_of_replicate >= 3:
                    data = numpy.vstack((data, self.mdv[fragment][number]['data']))
                if use == 'use':
                    use_array.append(1)
                else:
                    use_array.append(0)
        if self.number_of_replicate >= 3:
            data = data[1:,:]
        return id_array, numpy.array(ratio_array), numpy.array(std_array), use_array, self.observed_fragments, data

    def get_number_of_measurement(self):
        """
        Getter of number of measurement of the MDV data set.

        Returns
        --------
        number_of_measurement : Number of measurement

        Examples
        --------
        >>> number_of_measurement = mdv.get_number_of_measurement()

        See Also
        --------
        """
        used_fragments = set()
        counter = 0
        for fragment in self.observed_fragments:
            for i in self.mdv[fragment]:
                if self.mdv[fragment][i]['use'] == 'use':
                    used_fragments.add(fragment)
                    counter = counter + 1
        return counter-len(used_fragments)

    def compare_mdv(self, mdv):
        """
        Comparison of two MDV data

        Parameters
        --------
        mdv: an instanse of MdvData class for comparison

        Examples
        --------
        >>> mdv.compare_mdv(mdv_fitted)

        See Also
        --------
        """
        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                my =self.mdv[fragment][number]['ratio']
                you = mdv.mdv[fragment][number]['ratio']
                print("{0:8.8s} {1:1d} {2:6.5f} {3:6.5f} {4:-6.5f} {5:6.5f} {6:5.5s}".format(fragment, number, my ,you, my-you,self.mdv[fragment][number]['std'] ,self.mdv[fragment][number]['use']))
    def save(self, filename, format = "text"):
        """
        Save MDV data to the text/csv file

        Parameters
        ----------
        mdv : a instance of mfapy.mdv for templete
        filename : filename of MDV data with the format.

        Name	Spectrum	Select	MDV	Std
        Ala57	m0	1	0.566990778	0.000774686
        Ala57	m1	1	0.148623963	0.000774686
        Ala57	m2	1	0.039467636	0.000774686
        Ala57	m3	1	0.244917622	0.000774686

        format:
            'csv' : CSV.
            'text' : tab-deliminated text.

        Reterns
        ----------
        Boolean: True/False

        Examples
        --------
        >>> self.save_data(mdv, 'filename')

        See Also
        --------

        """
        #
        # preparation of data
        #
        Data = []
        Data.append(['Name','Spectrum','Select','MDV','Std'])
        for fragment in self.mdv.keys():
            for isotope_number in self.mdv[fragment].keys():
                if self.mdv[fragment][isotope_number]['use'] == 'use':
                    use_or_no = 1
                else:
                    use_or_no = 0
                Data.append([fragment,
                            str(isotope_number),
                            str(use_or_no),
                            str(self.mdv[fragment][isotope_number]['ratio']),
                            str(self.mdv[fragment][isotope_number]['std'])])
        try:
            with open(filename, 'w', newline='') as f:
                import csv
                if format == "text":
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerows(Data)
                if format == "csv":
                    writer = csv.writer(f, dialect='excel')
                    writer.writerows(Data)
        except:
            return False

        return True

    def load(self, filename, format = 'text',output = "normal"):
        """
        Load MDV data from the text file.

        Parameters
        ----------
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
        >>> mdv.load('filename', format = 'text',output = "normal")

        See Also
        --------

        """
        #
        #
        observed_fragments_set = set()
        with open(filename, 'r') as f:
            import csv
            if format == "text":
                reader = csv.reader(f, delimiter='\t')
            elif format == "csv":
                reader = csv.reader(f, dialect='excel')
            else:
                print("Unknown format!")
                return False

            for i, row in enumerate(reader):
                if output == "debug":
                    print(row)
                if i == 0: continue
                if len(row) != 5:

                    continue
                fragment_temp, number, use, ratio, std, *over = row
                if ratio == "":
                    ratio = 0
                if len(fragment_temp) > 0:
                    fragment = fragment_temp
                if len(fragment_temp) == 0:
                    continue
                if fragment not in self.mdv:
                    continue
                number = int(number.replace("m", ""))
                self.mdv[fragment][number]['ratio'] = float(ratio)
                if len(std) <= 0:
                    std = 1.0
                self.mdv[fragment][number]['std'] = float(std)
                if str(use) == "1":
                    observed_fragments_set.update([fragment])
                    self.mdv[fragment][number]['use'] = 'use'
                else:
                    self.mdv[fragment][number]['use'] = 'no'

        self.set_observed_fragments(observed_fragments_set)

        return True

class MdvTimeCourseData:
    def __init__(self):
        """
        Constructer of MdvTimeCourseData instance of the metabolic model

        Parameters
        ----------
        target_fragments : dictionary of target fragments.

        Examples
        --------
        >>> mdv_timecourse = MdvTimeCourseData(model)
        Mdv data can be accessed by:
        mdv_timecourse.mdv['fragment_name']['0'][ratio]

        See Also
        --------

        """
        #MdvData. __init__(self,model.target_fragments)
        #self.mdv = {}
        self.mdvtc ={}
        self.mode = "timecourse"

    def has_data(self, time, fragment, number):
        """
        Checker whether the MdvData instance has data of [number] th isotope of [fragment] mass spectra.

        Parameters
        ----------
        fragment: name of target_fragment
        number: number of isotope

        Reterns
        ----------
        Boolean: True\False

        Examples
        --------
        >>> if mdv.has_data(0, 'Leu85', 2): brabra

        See Also
        --------

        """
        if time in self.mdvtc.keys():
            if fragment in self.mdvtc[time].mdv.keys():
                if number in self.mdvtc[time].mdv[fragment].keys():
                    return True
        return False

    def set_data(self, time, fragment, number, ratio, std, use):
        """
        Setter of single MDV data

        Parameters
        ----------
        time
        fragment: name of target_fragment
        number: Number of isotope
        ratio: Relative abundance of mass spectra data
        std: Standard deviation of measured MDV data
        use: Data to be used for MDV comparison. 'use' or 'no'

        Reterns
        ----------
        Boolean: True\False

        Examples
        --------
        >>> mdv.set_data(0, 'Leu85', 2, 0.763, 0.0054, 'use')

        See Also
        --------

        """
        result = False
        if time in self.mdvtc.keys():
            result = self.mdvtc[time].set_data(fragment, number, ratio, std, use)
        return result

    def get_data(self, time, fragment, number):
        """
        Getter of single MDV data

        Parameters
        ----------
        fragment: name of target_fragment
        number: Number of isotope

        Reterns
        ----------
        ratio: Relative abundance of mass spectra data
        std: Standard deviation of measured MDV data
        use: Data to be used for MDV comparison. 'use' or 'no'

        Examples
        --------
        >>> ratio, std, use = mdv.get_data('Leu85', 2)

        See Also
        --------

        """
        if time in self.mdvtc.keys():
            if self.mdvtc[time].has_data(fragment, number):
                return self.mdvtc[time].get_data(fragment, number)
        return False, False, False

    def get_timecourse_data(self, fragment, number):
        """
        Getter of time course MDV data

        Examples
        --------
        >>> tc = mdv.get_timecourse_data('Leu85', 2)

        Parameters
        ----------
        fragment: name of target_fragment
        number: Number of isotope

        Reterns
        ----------
        tc: Array of time course of Mass isotopomer


        See Also
        --------

        """
        timepoints = self.get_timepoints()
        timecourse = []
        for time in timepoints:
            if self.mdvtc[time].has_data(fragment, number):
                ratio, stdev, use = self.mdvtc[time].get_data(fragment, number)
                timecourse.append(ratio)
        return timecourse

    def get_timepoints(self):
        """
        Getter of Time points

        Parameters
        ----------


        Reterns
        ----------
        timepoints: list of time points
        Examples
        --------
        >>> timepoints = mdv.get_timepoints()
        >>> timepoints
        >>> [0., 1., 5,]

        See Also
        --------

        """
        return list(sorted(self.mdvtc.keys()))

    def set_observed_fragments(self, fragments):
        """
        Setter of observed fragment. When new MdvData instance is generated 'observed_fragment' == 'target_fragment'.

        Parameters
        ----------
        fragments: List of fragment names to be used for MDV comparison

        Examples
        --------
        >>> mdv.set_observed_fragments(['Phe_85', 'Ala_57'])

        See Also
        --------

        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].set_observed_fragments(fragments)

    def add_time_point(self,time, mdv_instance):
        self.mdvtc[time] = mdv_instance


    def add_gaussian_noise(self, stdev, iteration, method = 'absolute'):
        """
        Add gausian noise to MdvTimeCourseData

        Parameters
        --------
        ratio: Intensity of noise;
        iteration: Number of sampleing
        method:
            relative: intentsity = (randn() * stdev + 1) * original_intensity
            absolute: intentsity = randn() * stdev) + original_intensity

        Returns
        --------
        boolean: True/False

        Examples
        --------
        >>> mdv.add_gaussian_noise(0.05, 5)

        See Also
        --------
        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].add_gaussian_noise(stdev, iteration, method)

    def set_std(self, value, method = 'absolute'):
        """
        Set standard deviation levels from mass sepctral intensity

        Parameters
        ----------
        method: Method to calculate stdev levels:
            'relative' (detault): Levels are set by stdev = [value] * signal intensity
            'abusolute': Levels are set by stdev = [value]

        Examples
        --------
        >>> mdv.set_std(0.05, method = 'absolute')

        See Also
        --------
        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].set_std(value, method)

    def set_mdvs_for_comparison(self, threshold_ratio, threshold_std = 1.0):
        """
        Setters of multiple data to be used for MDV comparison

        Parameters
        --------
        threshold_ratios: Mass spectral signals with normalized intensities above threshold_ratios are used to be used for MDV comparison
        threshold_std: Threshold for std data

        Returns
        --------
        boolean: True/False

        Examples
        --------
        >>> mdv.set_mdvs_for_comparison(0.05, 0.005)

        See Also
        --------
        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].set_mdvs_for_comparison(threshold_ratio, threshold_std)



    def generate_observed_mdv(self):
        """
        Generator of MDV related inforamtion in list format.
        This function is used in model.set_experiment()

        Returns
        --------
        id_array: List of ids
        ratio_array: List of ratio data
        std_array:List of std data
        use_array: List of use data
        observed_fragments: List of observed fragment
        data:

        Examples
        --------
        >>> id_array, ratio_array, std_array, use_array, observed_fragments, data = mdv.generate_observed_mdv()

        See Also
        --------
        """
        id_array = []
        ratio_array = []
        std_array = []
        use_array = []
        raw_data_array = []
        for time in sorted(self.mdvtc.keys()):
            ids, mdv_exp_original, mdv_std_original, mdv_use, target_emu_list, rawdata = self.mdvtc[time].generate_observed_mdv()
            id_array.extend(ids)
            ratio_array.extend(mdv_exp_original)
            std_array.extend(mdv_std_original)
            use_array.extend(mdv_use)
            raw_data_array.extend(rawdata)


        return id_array, numpy.array(ratio_array), numpy.array(std_array), use_array, target_emu_list, raw_data_array

    def get_number_of_measurement(self):
        """
        Getter of number of measurement of the MDV data set.

        Returns
        --------
        number_of_measurement : Number of measurement

        Examples
        --------
        >>> number_of_measurement = mdv.get_number_of_measurement()

        See Also
        --------
        # Revised 180527
        """
        num_of_meas = 0
        for time in self.mdvtc.keys():
            num_of_meas = num_of_meas + self.mdvtc[time].get_number_of_measurement()
        #
        return num_of_meas

    def save(self, filename, format = "text"):
        """
        Save MDV data to the text/csv file

        Parameters
        ----------
        mdv : a instance of mfapy.mdv for templete
        filename : filename of MDV data with the format.

        Name	Spectrum	Select	MDV	Std
        Ala57	m0	1	0.566990778	0.000774686
        Ala57	m1	1	0.148623963	0.000774686
        Ala57	m2	1	0.039467636	0.000774686
        Ala57	m3	1	0.244917622	0.000774686

        format:
            'csv' : CSV.
            'text' : tab-deliminated text.

        Reterns
        ----------
        Boolean: True/False

        Examples
        --------
        >>> self.save_data(mdv, 'filename')

        See Also
        --------

        """
        for time in self.mdvtc.keys():
            if format == "csv":
                save_filename = filename + str(int(time)) + "sec.csv"
            elif format == "text":
                save_filename = filename + str(int(time)) + "sec.txt"
            else:
                save_filename = filename + str(int(time)) + "sec.txt"
            self.mdvtc[time].save(save_filename, format)




def binomial(n,k):
    a=math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
    return a

def INV_correcting(num_N,num_Si,num_O,num_H,num_C,num_S):

    row = num_C + num_H + 2*num_O + num_N + num_S + 2*num_Si + 1
    column =num_C + num_H + 2*num_O + num_N +num_S + 2*num_Si + 1

    DerivC=numpy.zeros((row,column))
    DerivH=numpy.zeros((row,column))
    DerivO=numpy.zeros((row,column))
    DerivN=numpy.zeros((row,column))
    DerivS=numpy.zeros((row,column))
    DerivSi=numpy.zeros((row,column))
    #
    #   Cの補正
    #
    for i in range(row):
        for j in range(num_C+1):
            if j+i >= row :
                break
            DerivC[j+i,i] = binomial(num_C,j)*0.9893**(num_C - j)*0.0107**j
    #
    #
    #   H
    #
    for i in range(row):
        for j in range(num_H +1):
            if j+i >=row:
                break
            DerivH[j+i,i]=binomial(num_H,j)*0.99985**(num_H-j)*0.00015**j
    #
    #
    #   N
    #
    for i in range(row):
        for j in range(num_N+1):
            if j+i >= row :
                break
            DerivN[j+i,i] = binomial(num_N,j)*0.99632**(num_N - j)*0.00368**j
    #
    #
    #   O
    #

    #iは16Oの数、jは17Oの数

    for i in range(row):

        if num_O == 0:
            DerivO=numpy.identity(column)
            break

        for j in range(2*num_O +1):
            if j+i>= row:
                break

            tmp=j
            if tmp > num_O:
                tmp = num_O

            for s in range(tmp+1):
                for k in range(tmp - s +1):
                    if s + 2*k == j:
                        DerivO[j+i,i] = DerivO[j+i,i] +binomial(num_O,s)*binomial((num_O-s),k)*0.99757**(num_O-(s+k))*0.00038**s*0.00205**k
    #
    #
    #   Si
    #
    for i  in range(row):
        if num_Si == 0:
            DerivSi =numpy.identity(column)
            break

        for j in range(2*num_Si +1):
            if j+i >= row:
                break

            tmp=j
            if tmp > num_Si:
                tmp = num_Si


            for s in range(tmp+1):
                for k in range(tmp-s+1):
                    if s + 2*k ==j:
                        DerivSi[j+i,i] =DerivSi[j+i,i] + binomial(num_Si,s)*binomial((num_Si-s),k)*0.922297**(num_Si-(s+k))*0.046832**s*0.030871**k
    #
    #
    #   S
    #

    #iは16Oの数、jは17Oの数
    for i in range(row):
        if num_S==0:
            DerivS=numpy.identity(column)
            break

        for j in range(2*num_S+1):
            if j+i >= row:
                break

            tmp=j
            if tmp > num_S:
                tmp = num_S

            for s in range(tmp+1):
                for k in range(tmp-s +1):
                    if s+2*k == j:
                        DerivS[j+i,i] = DerivS[j+i,i] + binomial(num_S,s)*binomial((num_S-s),k)*0.9493**(num_S-(s+k))*0.0076**s*0.0429**k



    N =numpy.linalg.inv(DerivN)
    Si=numpy.linalg.inv(DerivSi)
    O=numpy.linalg.inv(DerivO)
    H=numpy.linalg.inv(DerivH)
    C=numpy.linalg.inv(DerivC)
    S=numpy.linalg.inv(DerivS)

    tmp_INV = C.dot(H).dot(O).dot(Si).dot(N).dot(S)

    return tmp_INV

def transition_matrix(num_N,num_Si,num_O,num_H,num_C,num_S):

    row = num_C + num_H + 2*num_O + num_N + num_S + 2*num_Si + 1
    column =num_C + num_H + 2*num_O + num_N +num_S + 2*num_Si + 1

    DerivC=numpy.zeros((row,column))
    DerivH=numpy.zeros((row,column))
    DerivO=numpy.zeros((row,column))
    DerivN=numpy.zeros((row,column))
    DerivS=numpy.zeros((row,column))
    DerivSi=numpy.zeros((row,column))
    #
    #   Cの補正
    #
    for i in range(row):
        for j in range(num_C+1):
            if j+i >= row :
                break
            DerivC[j+i,i] = binomial(num_C,j)*0.9893**(num_C - j)*0.0107**j
    #
    #
    #   H
    #
    for i in range(row):
        for j in range(num_H +1):
            if j+i >=row:
                break
            DerivH[j+i,i]=binomial(num_H,j)*0.99985**(num_H-j)*0.00015**j
    #
    #
    #   N
    #
    for i in range(row):
        for j in range(num_N+1):
            if j+i >= row :
                break
            DerivN[j+i,i] = binomial(num_N,j)*0.99632**(num_N - j)*0.00368**j
    #
    #
    #   O
    #

    #iは16Oの数、jは17Oの数

    for i in range(row):

        if num_O == 0:
            DerivO=numpy.identity(column)
            break

        for j in range(2*num_O +1):
            if j+i>= row:
                break

            tmp=j
            if tmp > num_O:
                tmp = num_O

            for s in range(tmp+1):
                for k in range(tmp - s +1):
                    if s + 2*k == j:
                        DerivO[j+i,i] = DerivO[j+i,i] +binomial(num_O,s)*binomial((num_O-s),k)*0.99757**(num_O-(s+k))*0.00038**s*0.00205**k
    #
    #
    #   Si
    #
    for i  in range(row):
        if num_Si == 0:
            DerivSi =numpy.identity(column)
            break

        for j in range(2*num_Si +1):
            if j+i >= row:
                break

            tmp=j
            if tmp > num_Si:
                tmp = num_Si


            for s in range(tmp+1):
                for k in range(tmp-s+1):
                    if s + 2*k ==j:
                        DerivSi[j+i,i] =DerivSi[j+i,i] + binomial(num_Si,s)*binomial((num_Si-s),k)*0.922297**(num_Si-(s+k))*0.046832**s*0.030871**k
    #
    #
    #   S
    #

    #iは16Oの数、jは17Oの数
    for i in range(row):
        if num_S==0:
            DerivS=numpy.identity(column)
            break

        for j in range(2*num_S+1):
            if j+i >= row:
                break

            tmp=j
            if tmp > num_S:
                tmp = num_S

            for s in range(tmp+1):
                for k in range(tmp-s +1):
                    if s+2*k == j:
                        DerivS[j+i,i] = DerivS[j+i,i] + binomial(num_S,s)*binomial((num_S-s),k)*0.9493**(num_S-(s+k))*0.0076**s*0.0429**k

    tmp_r = DerivC.dot(DerivH).dot(DerivO).dot(DerivSi).dot(DerivN).dot(DerivS)

    return tmp_r



