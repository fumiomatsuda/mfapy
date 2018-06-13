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


import numpy as numpy
import itertools

class CarbonSource:
    """
    Class for carbon source information.
    """
    def __init__(self, carbon_sources):
        """
        Constructer of CarbonSource instance from target_fragments information of the metabolic model

        Parameters
        ----------
        carbon_sources : dictionary of carbon sources information produced in model.carbon_sources

        nothing.

        Returns
        -------
        cs : CarbonSource instance

        Examples
        --------
        >>> cs = CarbonSource(arbon_sources)
        CarbonSource data can be accessed by:
        cs.cs['fragment_name']['IDV']

        See Also
        --------


        """
        self.cs = {}
        #for each line in mass data
        for compound in carbon_sources:
            self.cs[compound] = {
            'IDV': carbon_sources[compound]['IDV'][:],
            'size': carbon_sources[compound]['size']
            }
        self.mdv_carbon_sources = {}
        #Set MDVs of EMUs of each carbon sources
        self.generate_carbonsource_MDV()

    def show(self):
        """
        Show concents of CarbonSource instance

        Parameters
        ----------
        nothing.

        Returns
        -------
        nothing.

        Examples
        --------
        >>> cs.show()
        Name: Asp
        Carbon number: 4
        #0000	1.0
        Name: AcCoA
        Carbon number: 2
        #00	0.5
        #10	0.25
        #11	0.25

        See Also
        --------

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
        """
        Generator of a dictionary of MDVs of all EMUs

        Parameters
        ----------
        nothing

        Returns
        -------
        mdvs : Dictionary of MDVs of all EMUs

        Examples
        --------
        >>> mdvs = cs.generate_dict()


        See Also
        --------

        """
        mdv_cs = {}
        #for each line in mass data
        for label in self.mdv_carbon_sources:
            mdv_cs[label] = self.mdv_carbon_sources[label][:]
        return mdv_cs

    def set_all_isotopomers(self, compound, list, correction = 'no'):
        """
        Setter of IDV data by list of all mass isotopomer distribution

        Parameters
        ----------
        compound: Name of carbon source
        list: list of mass isotopomer distribution (from 000, 001, 010, to 111)
        correction: (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns
        -------
        Booleans.


        Examples
        --------
        >>> cs.set_all_isotopomers('AcCoA', [0.3, 0.3, 0.3, 0.1])


        See Also
        --------

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

    def set_each_isotopomer(self, compound, dict, correction = 'no'):
        """
        Setter of IDV data by distribution of selected mass isotopomers

        Parameters
        ----------
        compound: Name of carbon source
        dict: Dictionary of mass isotopomer and its relative abundance {'#111': 0.5, '#001': 0.5}
        correction: (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns
        -------
        Booleans.

        Examples
        --------
        >>> cs.set_each_isotopomers('AcCoA', {'#11':0.5, '#10':0.25}, correction = 'yes')


        See Also
        --------

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
        """
        Setter of IDV data by distribution of selected mass isotopomers

        Parameters
        ----------
        compound: Name of carbon source
        dict: Dictionary of lableled compound name and its relative abundance {'[1_13C]glucose': 0.5, '[U_13C]glucose':0.5}
        correction: (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns
        -------
        Booleans.

        Examples
        --------
        >>> cs.set_labeled_compounds('Glc',{'[1_13C]glucose': 0.5, '[U_13C]glucose':0.5}, correction = 'yes')

        Availables
        --------
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
        """
        Generate MDVs of all EMUs of each carbon source.
        This function is called in the mfapy.metabolicmodel.MetabolicModel

        Parameters
        ----------
        carbonsource: List of carbon source metabolite (optional)
        correction: (yes/no) Correction of isotopomer distribution considering natural 13C occurence

        Returns
        -------
        Dictionary of mdv data of carbon source.

        Examples
        --------
        >>> cs.generate_carbonsource_MDV()


        See Also
        --------

        History
        -------
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
                    bin = list('{0:010b}'.format(t))[10-size:10]# Max 10
                    bin.reverse()
                    #
                    #['0', '0', '0', '0'] <= Label pattern = #0000
                    #
                    # Count nonlabled carbon (0)
                    #
                    zero_count = float(sum([1 for x in bin if x == '0']))
                    # zero_count = 4
                    for i in range(0, len(bin)):
                        # scan ['0', '0', '0', '0']
                        #
                        # If there is one natural 13C among 12C
                        #
                        isotope = bin[i]
                        if isotope == "0":
                            bin_temp = list(bin) #copy
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
            # EMU generation
            #
            for emu_size in range(1, size+1):
                #
                # Generate all EMUs by itertools.combinations
                #
                for c in itertools.combinations(range(1,size+1),emu_size):
                    MID = [0] * (int(emu_size) + 1)#Mass isotopomer distribution
                    filter = [0]  * int(size);# Carbons for EMU
                    #
                    # Name of EMUs.
                    #
                    numbers =  "".join(map(str,c))
                    emu = compound + "_" + numbers
                    #
                    # Set carbon positions
                    #
                    for i in numbers:
                        filter[int(i)-1] = 1
                    for t in range((2 ** int(size))):
                        bin = list(format (t, '08b'))#for all isotopomers
                        bin.reverse()
                        #Cals number of 13C in EMU
                        number = sum([int(bin[x]) for x in range(len(filter)) if filter[x] == 1])
                        #Integraton。
                        MID[number] += IDV_array[t]
                    #Normalised to 1.0
                    sum_MID = sum(MID)
                    self.mdv_carbon_sources[emu] = [x/sum_MID for x in MID]
        return(self.mdv_carbon_sources)

def main():
    pass

if __name__ == '__main__':
    main()
