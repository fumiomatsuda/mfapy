#!/usr/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        carbonsource.py
# Purpose:     CarbonSource class in mfapy
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------
"""carbonsource.py:CarbonSource class in mfapy

The module includes::

    CarbonSource class

Todo:
    * Subtraction of natural isotope

"""

import numpy as numpy
import itertools

class CarbonSource:
    """Class for carbon source information.

    Instance of this class is generated in the MetabolicModel instance

    """
    def __init__(self, carbon_sources, emu_dict):
        """Constructer of CarbonSource instance from target_fragments information of the metabolic model

        Args:
            carbon_sources (dict): dictionary of carbon sources information in model.carbon_sources.

        Returns:
            instance : CarbonSource instance

        Examples:
            >>> cs = CarbonSource(arbon_sources)
            CarbonSource data can be accessed by:
            cs.cs['fragment_name']['IDV']


        """
        self.cs = {}
        #for each line in mass data
        for compound in carbon_sources:
            self.cs[compound] = {
            'IDV': carbon_sources[compound]['IDV'][:],
            'size': carbon_sources[compound]['size']
            }
        self.emu_dict = {}
        for compound, position_list in emu_dict.items():
            self.emu_dict[compound] = position_list
        self.mdv_carbon_sources = {}
        #Set MDVs of EMUs of each carbon sources
        self.generate_carbonsource_MDV()

    def show(self):
        """Method to display contents of CarbonSource instance

        Args:
            Not required.

        Returns:
            Nothing.

        Examples:

            >>> cs.show()
            Name: Asp
            Carbon number: 4
            #0000	1.0
            Name: AcCoA
            Carbon number: 2
            #00	0.5
            #10	0.25
            #11	0.25


        """
        for compound in self.cs:
            IDV = self.cs[compound]['IDV'][:]
            size = self.cs[compound]['size']
            print('Name: ' +  compound)
            print('Carbon number: ' +  str(size))
            for i in range(0, 2**size):
                if IDV[i] > 0.0:
                    strbin = format(i, 'b')# strbin = format(10,'b') #'1010'
                    strbin ="0" * (size - len(strbin)) + strbin# strbin = "0" * (6-4) + '1010' = '001010'
                    print("#"+strbin[::-1]+"\t", IDV[i])# #010100 reverse 0/1
            for emu in sorted(self.mdv_carbon_sources.keys()):
                if compound in emu:
                    print(emu+"\t", self.mdv_carbon_sources[emu], "sum "+str(sum(self.mdv_carbon_sources[emu])))

    def generate_dict(self):
        """Generator of a dictionary of MDVs of all EMUs of carbon sources

        Args:
            Not required.

        Returns:
            dict : Dictionary of MDVs of all EMUs of carbon sources

        Examples:
            >>> mdvs = cs.generate_dict()



        """
        mdv_cs = {}
        #for each line in mass data
        for label in self.mdv_carbon_sources:
            mdv_cs[label] = self.mdv_carbon_sources[label][:]
        return mdv_cs

    def set_all_isotopomers(self, compound, list, correction = 'no'):
        """Setter of IDV data by list of all mass isotopomer distribution

        Args:
            compound (str): Name of carbon source.

            list (array): list of mass isotopomer distribution (from 000, 001, 010, to 111)

            correction (str): (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns:
            Booleans: True/False


        Examples:
            >>> cs.set_all_isotopomers('AcCoA', [0.3, 0.3, 0.3, 0.1])



        """
        if compound not in self.cs:
            print("False compound name")
            return False
        if len(list) != 2**self.cs[compound]['size']:
            print("False list length")
            return False
        if abs(1.0-sum(list)) > 0.00001:
            print("Sum of list is not 1.0")
            return False
        self.cs[compound]['IDV'] = list[:]
        self.generate_carbonsource_MDV(carbonsource = [compound], correction = correction)
        return True

    def set_carbonsources(self, filename, correction = 'no', format = 'text',output = "normal"):
        """Setter of isotopomer data of multiple carbon sourses from text file.

        Args:
            filename (str): filename of MDV data with following format::

                Name	Isotopomer	Ratio
                Asp	#0000	0.5
                Asp	#1111	0.5
                AcCoA	#00	0.5
                AcCoA	#11	0.5

            correction (str) : (yes/no) Correction of isotopomer distribution considering natural 13C occurence

            format (str) : "text" (defalut) or "csv"

            output (str) : "normal" (defalut) or "debug"

        Returns:
            Boolean: True/False

        Examples:
            >>> cs2.set_isotopomers_from_file('Example_1_carbonsource2.txt', correction = "yes")



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
            dict = {}

            for i, row in enumerate(reader):
                if output == "debug":
                    print(row)
                if i == 0: continue
                if len(row) != 3:
                    continue
                fragment, isotopomer, ratio, *over = row
                if ratio == "":
                    ratio = 0

                if fragment not in self.cs:
                    continue
                if fragment not in dict:
                    dict[fragment] = {}
                dict[fragment][isotopomer] = float(ratio)
        for fragment in dict:
            self.set_each_isotopomer(fragment, dict[fragment], correction = correction)
        return True

    def set_each_isotopomer(self, compound, dict, correction = 'no'):
        """Setter of IDV data of selected mass isotopomers

        Args:
            compound (str): An ID of the carbon source metabolite. The ID must be defined in the model definition file.

            dict (dict): A dictionary of isotopomers and their ratios. For example, for SubsGlc with six carbons, the keys '#000000', '#100000', and '#111111' represent non-labeled glucose, [1-13C]glucose, and [U-13C]glucose, respectively. The values (e.g., 0.02, 0.7, and 0.28) denote their relative abundances. Ex. {'#000000': 0.02, '#100000': 0.7, '#111111': 0.28 }
		
            correction (str): (yes/no) If 'yes', the composition of the isotopomer is automatically adjusted considering the occurrence of natural 13C.

        Returns:
            Boolean: True/False

        Examples:
            >>> cs.set_each_isotopomer('SubsGlc',{'#000000': 0.02, '#100000': 0.7, '#111111': 0.28 }, correction = 'yes')




        """

        #check compound name
        #
        if compound not in self.cs:
            print("False compound name")
            return False

        size = self.cs[compound]['size']
        #
        #check compound name
        #
        sum = 0.0
        for isotopomer in dict:
            sum = sum + dict[isotopomer]
        if sum > 1.0:
            print("Sum isctopomers is over 1.0")
            return False

        for isotopomer in dict:
            sum = 0
            for letter in isotopomer:
                if letter == "0": sum = sum + 1
                if letter == "1": sum = sum + 1
            if sum != size:
                print("False isotopomer" + isotopomer)
                return False

        #
        # generate IDV
        #
        self.cs[compound]['IDV'] = [0] * 2**size

        for isotopomer in dict:
            str_bin = ""
            for letter in isotopomer[::-1]:
                if letter == "0": str_bin = str_bin + "0"
                if letter == "1": str_bin = str_bin + "1"

            self.cs[compound]['IDV'][int(str_bin, 2)] = dict[isotopomer]
        self.generate_carbonsource_MDV(carbonsource = [compound], correction = correction)
        return True

    def set_labeled_compounds(self, compound, dict, correction = 'no'):
        """Setter of IDV data by distribution of selected mass isotopomers

        Following symbols are available::

            "[13C]CO2": "#1",
            "[12C]CO2": "#0",
            "[13C]THF": "#1",
            "[12C]THF": "#0",
            "[1-13C]glucose": "#100000",
            "[2-13C]glucose": "#010000",
            "[1,2-13C]glucose": "#110000",
            "[U-13C]glucose": "#111111",
            "[1-13C]glutamine": "#10000",
            "[2-13C]glutamine": "#01000",
            "[5-13C]glutamine": "#00001",
            "[U-13C]glutamine": "#11111",

        Args:
            compound (str): Name of carbon source.

            dict (dict): Dictionary of mass isotopomer and its relative abundance {'#111': 0.5, '#001': 0.5}

            correction (str): (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns:
            Boolean: True/False

        Examples:
            >>> cs.set_labeled_compounds('Glc',{'[1_13C]glucose': 0.5, '[U_13C]glucose':0.5}, correction = 'yes')


        """
        #
        # please enrich this list
        #
        labeledcompound_dict = {
        "[13C]CO2": "#1",
        "[12C]CO2": "#0",
        "[13C]THF": "#1",
        "[12C]THF": "#0",
        "[1-13C]glucose": "#100000",
        "[2-13C]glucose": "#010000",
        "[1,2-13C]glucose": "#110000",
        "[U-13C]glucose": "#111111",
        "[1-13C]glutamine": "#10000",
        "[2-13C]glutamine": "#01000",
        "[5-13C]glutamine": "#00001",
        "[U-13C]glutamine": "#11111",
        }


        #check compound name
        #
        if compound not in self.cs:
            print("False compound name")
            return False

        size = self.cs[compound]['size']
        #
        #check compound name
        #
        sum = 0.0
        for isotopomer in dict:
            sum = sum + dict[isotopomer]
        if sum > 1.0:
            print("Sum isotopomers is over 1.0")
            return False

        for labeledcompound in dict:
            if labeledcompound not in labeledcompound_dict:
                print("False lableled compound name "+labeledcompound)
                return False
        for labeledcompound in dict:
            sum = 0
            for letter in labeledcompound_dict[labeledcompound]:
                if letter == "0": sum = sum + 1
                if letter == "1": sum = sum + 1
            if sum != size:
                print("False isotopomer" + isotopomer)
                return False

        #
        # IDV????
        #
        self.cs[compound]['IDV'] = [0] * 2**size

        for labeledcompound in dict:
            str_bin = ""
            for letter in labeledcompound_dict[labeledcompound]:
                if letter == "0": str_bin = str_bin + "0"
                if letter == "1": str_bin = str_bin + "1"

            self.cs[compound]['IDV'][int(str_bin, 2)] = dict[isotopomer]
        self.generate_carbonsource_MDV(carbonsource = [compound], correction = correction)
        return True


    def generate_carbonsource_MDV(self, carbonsource = [], correction = 'no'):
        """Generator of MDVs of all EMUs of each carbon source.

        This function is called in the mfapy.metabolicmodel.MetabolicModel

        Args:
            carbonsource (array): List of carbon source metabolite (optional)

            correction (str): (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns:
            dict: Dictionary of mdv data of carbon source.

        Examples:

            >>> cs.generate_carbonsource_MDV()

        History:

            Correcton of IDV by natural 13C method was improved.

            labeled by 2 13C and 3 13C was taken into consideration.

        """

        #
        # Initialization
        #
        stable_isotope_ratio = 0.0107
        if len(carbonsource)==0:
            carbonsource = self.cs.keys()
        #
        # For each carbon source
        #
        for compound in carbonsource:

            if compound not in self.cs:
                continue
            #
            #  correction by natural isotope abundance
            #
            size = self.cs[compound]["size"]
            if correction != 'no':
                #
                # Preparation of IDV array for the natural isotope correaction.
                # Asp, {'IDV': [1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'size': 4}
                #
                IDV_array = numpy.zeros((2 ** int(size)))
                #
                # For each idv..
                #
                for t in range((2 ** int(size))):
                    #
                    # Ignore when carbon source does not contain the IDV
                    #  'IDV': [1.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                    if self.cs[compound]["IDV"][t] == 0:
                        continue
                    #
                    # Generate label patterns of IDV
                    #
                    # ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
                    #['0', '0', '0', '0']
                    #bindata = list('{0:012b}'.format(t))[12-size:12]# Max 11
                    #bindata.reverse()
                    #
                    # New in ver 0.6.1
                    bindata = list('{:b}'.format(t))
                    bindata.reverse()
                    diffsize = size-len(bindata)
                    bindata.extend(['0'] * diffsize)

                    #
                    #['0', '0', '0', '0'] <= Label pattern = #0000
                    #
                    # Count nonlabled carbon (0)
                    #
                    zero_count = float(sum([1 for x in bindata if x == '0']))
                    # zero_count = 4
                    for i in range(0, len(bindata)):
                        # scan ['0', '0', '0', '0']
                        #
                        # If there is one natural 13C among 12C
                        #
                        isotope = bindata[i]
                        if isotope == "0":
                            bin_temp = list(bindata) #copy
                            bin_temp[i] = '1' # Change to 13C ['1', '0', '0', '0']
                            bin_temp.reverse()  # reverse ['0', '0', '0', '1']
                            number_temp = int("".join(bin_temp),2) # convert to number of natural isotope 3
                            IDV_array[number_temp] = IDV_array[number_temp] + \
                            self.cs[compound]["IDV"][t] * stable_isotope_ratio * ((1.0- stable_isotope_ratio) ** (zero_count-1.0))
                            # 0 + 1 * 0.0107 + (1- 0.0107) ** (1-1)
                            for i2 in range(len(bin_temp)):
                                # scan ['0', '0', '0', '1']
                                # If there are two natural 13C among 12C
                                #
                                isotope2 = bin_temp[i2]
                                if isotope2 == "0":
                                    bin_temp2 = list(bin_temp)
                                    bin_temp2[i2] = '1'
                                    bin_temp2.reverse()
                                    number_temp2 = int("".join(bin_temp2),2)
                                    IDV_array[number_temp2] = IDV_array[number_temp2] + \
                                    0.5 * self.cs[compound]["IDV"][t] * (stable_isotope_ratio ** 2.0) * ((1.0- stable_isotope_ratio) ** (zero_count-2.0))
                                    for i3 in range(len(bin_temp2)):
                                        #
                                        # If there are three natural 13C among 12C
                                        #
                                        isotope3 = bin_temp2[i3]
                                        if isotope3 == "0":
                                            bin_temp3 = list(bin_temp2)
                                            bin_temp3[i3] = '1'
                                            bin_temp3.reverse()
                                            number_temp3 = int("".join(bin_temp3),2)
                                            IDV_array[number_temp3] = IDV_array[number_temp3] +  self.cs[compound]["IDV"][t] / 6.0 * (stable_isotope_ratio ** 3.0) * ((1.0- stable_isotope_ratio) ** (zero_count-3.0))
                    IDV_array[t] = IDV_array[t] + self.cs[compound]["IDV"][t] * ((1.0- stable_isotope_ratio) ** zero_count)
            else:
                IDV_array = list(self.cs[compound]["IDV"])
            #
            # EMU generation New in ver 0.6.1
            #
            for emu in self.emu_dict:
                compound_emu = self.emu_dict[emu]["metabolite_name"]
                if compound != compound_emu:
                    continue
                position_list = self.emu_dict[emu]["position_list"]
                emu_size = len(position_list)

                MID = [0] * (int(emu_size) + 1)#Mass isotopomer distribution
                filtered = [0]  * int(size);# Carbons for EMU
                #
                # Name of EMUs.
                #
                #numbers =  ":".join(map(str,c))
                #emu = compound + "_" + numbers
                #
                # Set carbon positions
                #
                for i in position_list:
                    filtered[int(i)-1] = 1
                for t in range((2 ** int(size))):
                    #bindata = list(format (t, '012b'))#for all isotopomers
                    #bindata.reverse()
                    #
                    #
                    bindata = list('{:b}'.format(t))
                    bindata.reverse()
                    diffsize = size-len(bindata)
                    bindata.extend(['0'] * diffsize)
                    #Cals number of 13C in EMU
                    #print(bin, filter)
                    number = sum([int(bindata[x]) for x in range(len(filtered)) if filtered[x] == 1])
                    #Integraton
                    MID[number] += IDV_array[t]
                #Normalised to 1.0
                sum_MID = sum(MID)
                self.mdv_carbon_sources[emu] = [x/sum_MID for x in MID]
        return(self.mdv_carbon_sources)

def main():
    pass

if __name__ == '__main__':
    main()
