#!/usr/bin/python
# -*- coding: utf-8 -*-

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
"""mdv.py:MdvData class in mfapy

The module includes::

    MdvData class
    MdvTimeCourseData class

Todo:
    * I/O

"""

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
        self.formula = {}
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
            self.formula[fragment] = target_fragments[fragment]["formula"]
            for i in range(target_fragments[fragment]['number']):
                id = fragment+"_"+str(i)
                self.mdv[fragment][i] = {"id":id, "ratio":float(0.0), "std":float(0.0), "use": "use", "data":[]}

        self.fragments_for_mdv_calculation = sorted(list(self.fragments_for_mdv_calculation_set))
        self.observed_fragments = sorted(list(self.observed_fragments_set))

    def has_data(self, fragment, number):
        """Checker whether the MdvData instance has data of [number] th isotopomer in the MDV of [fragment].

        Args:

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            Boolean: True\False

        Examples:
            >>> if mdv.has_data('Leu85', 2): brabra

        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                return True
        return False

    def set_data(self, fragment, number, ratio, std, use):
        """Setter of [number] th isotopomer in the MDV of [fragment].

        Args:

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

            ratio (float): Relative abundance of mass spectra data

            std (float): Standard deviation of measured MDV data

            use (str): Data to be used for MDV comparison. 'use' or 'no'

        Returns:
            Boolean: True\False

        Examples:
            >>> mdv.set_data('Leu85', 2, 0.763, 0.0054, 'use')


        """
        if self.has_data(fragment, number):
            self.mdv[fragment][number]['ratio'] = ratio
            self.mdv[fragment][number]['std'] = std
            self.mdv[fragment][number]['use'] = use
            return True
        return False

    def get_data(self, fragment, number):
        """Getter of [number] th isotopomer in the MDV of [fragment].

        Args:

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            ratio (float) Relative abundance of mass spectra data

            std (float) Standard deviation of measured MDV data

            use (str) Data to be used for MDV comparison. 'use' or 'no'

        Examples:
            >>> ratio, std, use = mdv.get_data('Leu85', 2)


        """
        if self.has_data(fragment, number):
            return self.mdv[fragment][number]['ratio'], self.mdv[fragment][number]['std'], self.mdv[fragment][number]['use']
        return False, False, False

    def get_fragment_mdv(self, fragment):
        """Getter of  MDV data of [fragment].


        Args:

            fragment (str): name of target_fragment

        Returns:
            mdv_pattern (array) Array of MDV of fragment


        Examples:
            >>> mdv_pattern = mdv.get_fragment_mdv('Leu85')


        """
        mdv_pattern = []
        if fragment in self.mdv:
            idvlist = sorted(self.mdv[fragment].keys())
            for number in idvlist:
                mdv_pattern.append(self.mdv[fragment][number]['ratio'])
        return mdv_pattern

    def correct_natural_isotope(self, mode = "normal"):
        """Subtraction of natural isotope effect from fragments

        Molecular formula information of "target fragment" in the model definition file is used.

        Args:
            mode (str): "normal"(defalut) and "correction" (forced removal of negatives value and sum(ratio) is set at 1.0)

        Returns:
            Nothing

        Examples:
            >>> mdv.correct_natural_isotope()

        """

        for fragment in self.mdv:
            if not self.formula[fragment] == "":
                pattern = self.get_fragment_mdv(fragment)
                num = len(pattern)
                matrix = INV_correcting(self.formula[fragment])[0:num,0:num]
                transformed = numpy.dot(matrix, pattern)
                if mode == "correction":
                    total = sum([x for x in transformed if x > 0])
                    transformed_temp = []
                    for x in transformed:
                        if x < 0:
                            transformed_temp.append(0)
                        else:
                            transformed_temp.append(x/total)
                    transformed = transformed_temp
                for i, ratio in enumerate(transformed):
                    self.mdv[fragment][i]["ratio"] =  ratio
        return


    def add_natural_isotope(self):
        """Addition of natural isotope effect to fragments

        Molecular formula information in "target fragment" of model file is used.

        Args:
            Not required.

        Returns:
            Nothing.


        Example:
            >>> mdv.add_natural_isotope()
        """

        for fragment in self.mdv:
            if not self.formula[fragment] == "":
                pattern = self.get_fragment_mdv(fragment)
                num = len(pattern)
                matrix = transition_matrix(self.formula[fragment])[0:num,0:num]
                transformed = numpy.dot(matrix, pattern)
                for i, ratio in enumerate(transformed):
                    self.mdv[fragment][i]["ratio"] =  ratio
        return

    def set_observed_fragments(self, fragments):
        """Setter of observed fragment.

        When new MdvData instance is generated 'observed_fragment' == 'target_fragment'.

        Args:
            fragments (Array): List of fragment names calculated in calmdv fucntion

        Returns:
            Nothing.

        Examples:
            >>> mdv.set_observed_fragments(['Phe_85', 'Ala_57'])

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
        """
        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                if fragment in self.observed_fragments:
                    self.set_mdv_for_comparison(fragment, number)
                else:
                    self.set_mdv_for_ignore(fragment, number)
        """

        return True


    def get_fragments_for_mdv_calculation(self):
        """Getter of fragment list whose MDVs are calculated in calmdv fucntion

        Args:
            Not required.

        Returns:
            array: List of fragment names calculated in calmdv fucntion


        Examples:
            >>> list_of_fragment = mdv.get_fragments_for_mdv_calculation()


        """
        return list(self.fragments_for_mdv_calculation)

    def set_mdv_for_comparison(self, fragment, number):
        """Setter of mass data to be used for MDV comparison

        Args:
            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            Nothing.

        Examples:
            >>> mdv.set_mdv_for_comparison(fragment, number)

        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                self.mdv[fragment][number]['use'] = 'use'
                return True
        return False

    def set_unused_mdv_for_comparison(self, fragment, number):
        """Setter of mass data to be ignored for MDV comparison

        Args:
            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            Nothing.


        Examples:
            >>> mdv.set_unused_mdv_for_comparison(fragment, number)


        """
        if fragment in self.mdv:
            if number in self.mdv[fragment]:
                self.mdv[fragment][number]['use'] = 'no'
                return True
        return False

    def set_std(self, value, method = 'absolute'):
        """Setter of standard deviation level

        Args:
            value (float): level of standard deviation.

            method (str): Method to calculate stdev levels.

                * 'relative': Levels are set by stdev = [value] * signal intensity.

                * 'absolute' (detault) : Levels are set by stdev = [value].

        Returns:
            Nothing.

        Examples:
            >>> mdv.set_std(0.05, method = 'absolute')


        """
        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                if method == 'absolute':
                    self.mdv[fragment][number]['std'] = value * 1.0
                else:
                    self.mdv[fragment][number]['std'] = self.mdv[fragment][number]['ratio'] * value

    def check(self, output = "debug"):
        """Function to check missing values in the MDV data before MDV comparison

        Args:
            output (str): Show details when "debug" mode

        Returns:
            boolean: True/False

        Examples:
            >>> mdv.check()

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
        """Setters of multiple data to be used for MDV comparison

        Args:
            threshold_ratios (float): MDV data whose ratio (0.0-1.0) is above threshold_ratios are used for MDV comparison
            threshold_std (float): MDV data whose std is below threshold_std are used for MDV comparison

        Returns:
            boolean: True/False

        Examples:
            >>> mdv.set_mdvs_for_comparison(0.05, 0.005)

        """
        for fragment in self.fragments_for_mdv_calculation:
            for number in self.mdv[fragment]:
                self.set_unused_mdv_for_comparison(fragment, number)
                if fragment in self.observed_fragments:
                    ratio, std, use = self.get_data(fragment, number)
                    if (ratio > threshold_ratio) and (std < threshold_std):
                        self.set_mdv_for_comparison(fragment, number)

    def add_gaussian_noise(self, stdev, iteration = 1, method = 'absolute', normalize = "on"):
        """Addition of gausian noise to MDV data

        Args:
            stdev (float): noise level (standard deviation of normal distribution)

            iteration (int): number of sampleing (experimental)

            method (str): method to add gausian noise

                * "relative"  intentsity = (randn() * stdev + 1) * original_intensity

                * "absolute" (default) intentsity = randn() * stdev) + original_intensity

            normalize (str): on (default)/off sum of ratio is set to 1.0

        Returns:
            boolean: True/False

        Examples:
            >>> mdv.add_gaussian_noise(0.01)

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
                if normalize == "on":
                    sumvalue = sum(noise[:,i])
                    noise[:,i] = noise[:,i] / sumvalue
            for number in range(number_of_mass_data):
                self.mdv[fragment][number]['ratio']= sum(noise[number,:])/iteration
                #self.mdv[fragment][number]['std'] = numpy.std(noise[number,:])
                self.mdv[fragment][number]['data'] = numpy.array(noise[number,:])

    def generate_observed_mdv(self):
        """Generator of MDV related inforamtion in list format.

        This function is onely used in model.set_experiment()

        Args:
            Nor required.

        Returns:
            * id_array (array) List of ids

            * ratio_array (array) List of ratio data

            * std_array (array) List of std data

            * use_array (array) List of use data

            * observed_fragments (array) List of observed fragments

            * data (array) List of other information

        Examples:
            >>> id_array, ratio_array, std_array, use_array, observed_fragments, data = mdv.generate_observed_mdv()

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
        """Getter of number of measurement of the MDV data set.

        Args:
            Nor required.

        Returns:
            int: Number of measurement

        Examples:
            >>> number_of_measurement = mdv.get_number_of_measurement()

        """
        used_fragments = set()
        counter = 0
        for fragment in self.observed_fragments:
            num_of_isotope = 0
            used_counter = 0
            for i in self.mdv[fragment]:
                num_of_isotope = num_of_isotope + 1
                if self.mdv[fragment][i]['use'] == 'use':

                    counter = counter + 1
                    used_counter = used_counter + 1
            if num_of_isotope == used_counter:
                used_fragments.add(fragment)
        return counter-len(used_fragments)

    def comparison_with_another_mdv(self, mdv):
        """Comparison of two MDV data

        Args:
            mdv: instanse of another MdvData class for comparison

        Returns:
            Nothing

        Examples:
            >>> mdv.compare_mdv(mdv_fitted)

        """
        for fragment in self.mdv:
            for number in self.mdv[fragment]:
                my =self.mdv[fragment][number]['ratio']
                you = mdv.mdv[fragment][number]['ratio']
                print("{0:8.8s} {1:1d} {2:6.5f} {3:6.5f} {4:-6.5f} {5:6.5f} {6:5.5s}".format(fragment, number, my ,you, my-you,self.mdv[fragment][number]['std'] ,self.mdv[fragment][number]['use']))
    def save(self, filename, format = "text"):
        """Method to save MDV data in text/csv file

        Args:
            filename (str): filename of MDV data with the format::

                Name	Spectrum	Select	MDV	Std
                Ala57	m0	1	0.566990778	0.000774686
                Ala57	m1	1	0.148623963	0.000774686
                Ala57	m2	1	0.039467636	0.000774686
                Ala57	m3	1	0.244917622	0.000774686

            format (str): file format:

                * 'csv' : CSV.
                * 'text' (defalut) : tab-deliminated text.

        Returns:
            Boolean: True/False

        Examples:
            >>> mdv.save('filename', format = "csv")


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
        """Method to load MDV data from text/csv file

        Args:
            filename (str): filename of MDV data with the format::

                Name	Spectrum	Select	MDV	Std
                Ala57	m0	1	0.566990778	0.000774686
                Ala57	m1	1	0.148623963	0.000774686
                Ala57	m2	1	0.039467636	0.000774686
                Ala57	m3	1	0.244917622	0.000774686

            format (str): file format:

                * 'csv' : CSV.
                * 'text' (defalut) : tab-deliminated text.

            output : "normal" (defalut) or "debug"


        Returns:
            Boolean: True/False

        Examples:
            >>> mdv.load('filename', format = 'text',output = "normal")

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
                if number not in self.mdv[fragment]:
                    continue

                self.mdv[fragment][number]['ratio'] = float(ratio)
                if len(std) <= 0:
                    std = 1.0
                self.mdv[fragment][number]['std'] = float(std)
                if str(use) == '1':

                    observed_fragments_set.update([fragment])
                    self.mdv[fragment][number]['use'] = "use"


                else:
                    self.mdv[fragment][number]['use'] = "no"


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
        """Checker whether the MdvData instance has data of [number] th isotopomer in the MDV of [fragment] at [time].

        Args:
            time (float): time point

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            Boolean: True\False

        Examples:
            >>> if mdv.has_data(0, 'Leu85', 2): brabra

        """
        if time in self.mdvtc.keys():
            if fragment in self.mdvtc[time].mdv.keys():
                if number in self.mdvtc[time].mdv[fragment].keys():
                    return True
        return False

    def set_data(self, time, fragment, number, ratio, std, use):
        """Setter of [number] th isotopomer in the MDV of [fragment] at [time].

        Args:
            time (float): time point

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

            ratio (float): Relative abundance of mass spectra data

            std (float): Standard deviation of measured MDV data

            use (str): Data to be used for MDV comparison. 'use' or 'no'

        Returns:
            Boolean: True\False

        Examples:
            >>> mdv.set_data(0, 'Leu85', 2, 0.763, 0.0054, 'use')
        """
        result = False
        if time in self.mdvtc.keys():
            result = self.mdvtc[time].set_data(fragment, number, ratio, std, use)
        return result

    def get_data(self, time, fragment, number):
        """Getter of [number] th isotopomer in the MDV of [fragment] at [time].

        Args:
            time (float): time point

            fragment (str): name of target_fragment

            number (int): [number] th isotopomer

        Returns:
            ratio (float) Relative abundance of mass spectra data

            std (float) Standard deviation of measured MDV data

            use (str) Data to be used for MDV comparison. 'use' or 'no'

        Examples:
            >>> ratio, std, use = mdv.get_data('Leu85', 2)


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
        """Setter of observed fragment.

        When new MdvData instance is generated 'observed_fragment' == 'target_fragment'.

        Args:
            fragments (Array): List of fragment names calculated in calmdv fucntion

        Returns:
            Nothing.

        Examples:
            >>> mdv.set_observed_fragments(['Phe_85', 'Ala_57'])

        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].set_observed_fragments(fragments)

    def add_time_point(self,time, mdv_instance):
        """Addition of new time point

        Args:
            time (float): time point

            mdv_instance (MdvData): MdvData instance

        Returns:
            nothing

        Examples:
            >>> mdv.add_time_point(3.0, mdv3)

        """

        self.mdvtc[time] = mdv_instance


    def add_gaussian_noise(self, stdev, iteration, method = 'absolute', normalize = "on"):
        """Addition of gausian noise to MDV data

        Args:
            stdev (float): noise level (standard deviation of normal distribution)

            iteration (int): number of sampleing (experimental)

            method (str): method to add gausian noise

                * "relative"  intentsity = (randn() * stdev + 1) * original_intensity

                * "absolute" (default) intentsity = randn() * stdev) + original_intensity

            normalize (str): on (default)/off sum of ratio is set to 1.0

        Returns:
            boolean: True/False

        Examples:
            >>> mdv.add_gaussian_noise(0.01)

        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].add_gaussian_noise(stdev, iteration, method, normalize)

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
        """Setters of multiple data to be used for MDV comparison

        Args:
            threshold_ratios (float): MDV data whose ratio (0.0-1.0) is above threshold_ratios are used for MDV comparison
            threshold_std (float): MDV data whose std is below threshold_std are used for MDV comparison

        Returns:
            boolean: True/False

        Examples:
            >>> mdv.set_mdvs_for_comparison(0.05, 0.005)

        """
        for time in self.mdvtc.keys():
            self.mdvtc[time].set_mdvs_for_comparison(threshold_ratio, threshold_std)



    def generate_observed_mdv(self):
        """Generator of MDV related inforamtion in list format.

        This function is onely used in model.set_experiment()

        Args:
            Nor required.

        Returns:
            * id_array (array) Array of list of ids at each time point

            * ratio_array (array) Array of list of ratio data at each time point

            * std_array (array) Array of list of std data at each time point

            * use_array (array) Array of list of use data at each time point

            * observed_fragments (array) Array of list of observed fragments at each time point

            * data (array) Array of list of other information at each time point

        Examples:
            >>> id_array, ratio_array, std_array, use_array, observed_fragments, data = mdv.generate_observed_mdv()

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
        """Getter of number of measurement of the MDV data set.

        Args:
            Nor required.

        Returns:
            int: Number of measurement

        Examples:
            >>> number_of_measurement = mdv.get_number_of_measurement()

        History
            Revised at 27/5/2018
        """
        num_of_meas = 0
        for time in self.mdvtc.keys():
            num_of_meas = num_of_meas + self.mdvtc[time].get_number_of_measurement()
        #
        return num_of_meas

    def save(self, filename, format = "text"):
        """Method to save MDV data in text/csv file

        Files of each timepoint are generated with names filename+"timepoint".txt

        Args:
            filename (str): filename of MDV data with the format::

                Name	Spectrum	Select	MDV	Std
                Ala57	m0	1	0.566990778	0.000774686
                Ala57	m1	1	0.148623963	0.000774686
                Ala57	m2	1	0.039467636	0.000774686
                Ala57	m3	1	0.244917622	0.000774686

            format (str): file format:

                * 'csv' : CSV.
                * 'text' (defalut) : tab-deliminated text.

        Returns:
            Boolean: True/False

        Examples:
            >>> mdv.save('filename', format = "csv")


        """
        #
        for time in self.mdvtc.keys():
            if format == "csv":
                save_filename = filename + str(int(time)) + ".csv"
            elif format == "text":
                save_filename = filename + str(int(time)) + ".txt"
            else:
                save_filename = filename + str(int(time)) + ".txt"
            self.mdvtc[time].save(save_filename, format)




def binomial(n,k):

    """Method to calc binomial coefficients.

        This function calculates nCk.

    Args:
        n (int): n
        k (int): k


    Returns:
        a (int): math.factorial(n)/(math.factorial(k)*math.factorial(n-k))



    """

    a=math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
    return a

def INV_correcting(formula):
    """Method to create a matrix to deduct natural isotope effect of given formula.

    Args:
        formula (str): chemical formula
            * Now this function supports "C", "H", "N", "O", "P", "S", "Si".

    Returns:
        Numpy matrix


    """

    atomdic = count_atom_number(formula)
    num_C = atomdic["C"]
    num_H = atomdic["H"]
    num_N = atomdic["N"]
    num_O = atomdic["O"]
    num_P = atomdic["P"]
    num_S = atomdic["S"]
    num_Si = atomdic["Si"]



    row = num_C + num_H + 2*num_O + num_N + num_S + num_P + 2*num_Si + 1
    column =num_C + num_H + 2*num_O + num_N +num_S + num_P + 2*num_Si + 1

    DerivC=numpy.zeros((row,column))
    DerivH=numpy.zeros((row,column))
    DerivO=numpy.zeros((row,column))
    DerivN=numpy.zeros((row,column))
    DerivS=numpy.zeros((row,column))
    DerivSi=numpy.zeros((row,column))
    DerivP=numpy.zeros((row,column))
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
            DerivH[j+i,i]=binomial(num_H,j)*0.999885**(num_H-j)*0.000115**j
    #
    #
    #   N
    #
    for i in range(row):
        for j in range(num_N+1):
            if j+i >= row :
                break
            DerivN[j+i,i] = binomial(num_N,j)*0.99636**(num_N - j) * 0.00364**j
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
    #   P
    #
    for i in range(row):
        for j in range(num_P+1):
            if j+i >= row :
                break
            DerivP[j+i,i] = binomial(num_P,j)*1.000**(num_P - j)*0.0000**j
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
    P=numpy.linalg.inv(DerivP)

    tmp_INV = C.dot(H).dot(O).dot(Si).dot(N).dot(S)

    return tmp_INV

def transition_matrix(formula):
    """Method to create a matrix to add natural isotope effect of given formula.

    Args:
        formula (str): chemical formula
            * Now this function supports "C", "H", "N", "O", "P", "S", "Si".

    Returns:
        Numpy matrix


    """

    atomdic = count_atom_number(formula)
    num_C = atomdic["C"]
    num_H = atomdic["H"]
    num_N = atomdic["N"]
    num_O = atomdic["O"]
    num_P = atomdic["P"]
    num_S = atomdic["S"]
    num_Si = atomdic["Si"]


    row = num_C + num_H + 2*num_O + num_N + num_S + num_P + 2*num_Si + 1
    column =num_C + num_H + 2*num_O + num_N +num_S + num_P + 2*num_Si + 1

    DerivC=numpy.zeros((row,column))
    DerivH=numpy.zeros((row,column))
    DerivO=numpy.zeros((row,column))
    DerivN=numpy.zeros((row,column))
    DerivS=numpy.zeros((row,column))
    DerivSi=numpy.zeros((row,column))
    DerivP=numpy.zeros((row,column))
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
            DerivH[j+i,i]=binomial(num_H,j)*0.999885**(num_H-j)*0.000115**j
    #
    #
    #   N
    #
    for i in range(row):
        for j in range(num_N+1):
            if j+i >= row :
                break
            DerivN[j+i,i] = binomial(num_N,j)*0.99636**(num_N - j) * 0.00364**j
    #
    #
    #   P
    #
    for i in range(row):
        for j in range(num_P+1):
            if j+i >= row :
                break
            DerivP[j+i,i] = binomial(num_P,j)*1.000**(num_P - j)*0.0000**j
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

    tmp_r = DerivC.dot(DerivH).dot(DerivO).dot(DerivSi).dot(DerivN).dot(DerivS).dot(DerivP)

    return tmp_r

def count_atom_number(formula):
    """Method to extract atom numbers from given moleculr formula

    Args:
        formula (str): chemical formula
            * Now this function supports "C", "H", "N", "O", "P", "S", "Si".

    Returns:
        Dict

    Examples:
        >>> print(mdv.count_atom_number("C6H12O6"))
        >>> {"C": 6, "H": 12, "O": 6}


    """


    import re
    atoms = {}
    for atom in ["C", "H", "N", "O", "P", "S", "Si"]:
        arrey = re.findall(atom+'([0-9]*)', formula)
        cnum = 0
        for i in arrey:
            if  i == "":
                cnum = cnum+1
            else:
                cnum = cnum + int(i)
        atoms[atom] = cnum
    return atoms



