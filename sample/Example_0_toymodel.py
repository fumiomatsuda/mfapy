#-------------------------------------------------------------------------------
# Name:        mfapy Example_0_toymodel.py
#              Test code for mass isotopomer distribution vector (MDV) calculation related functions.
#              Toy model used in this code in deribed from Antoniewicz et al Metab. Eng. 2007, 9, 68-86.
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2021
# Licence:     MIT license
#-------------------------------------------------------------------------------
import mfapy

if __name__ == '__main__':
    #
    # Load metabolic model from txt file
    #
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("example_0_toymodel_model.txt", format = "text")
    #
    # Construction of MetabolicModel instance
    #
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    #
    # Configurations
    #
    model.set_configuration(callbacklevel = 0) #
    #
    # Modifiction of constraint "type" of reaction "v1".
    #
    model.set_constrain('reaction', 'v1', "fixed", 100.0, 1)
    #
    # self.update() is required after modification of any "type" of constraint
    #
    model.update()
    #
    # Load metabolic state from file
    #
    flux = model.load_states("Example_0_toymodel_status.csv", format = 'csv')
    #
    # Generation of instances of CarbonSource class from model
    #
    cs = model.generate_carbon_source_templete()
    cs2 = model.generate_carbon_source_templete()
    #
    # Isotope labelling information of carbon sources is set by ratio of all isotopomer (#00, #01, #10, #11)
    #
    cs.set_all_isotopomers('AcCoA', [0.5, 0.0, 0.25, 0.25])#a mixture of 25% [2-13C]AcCoA and 25%[1,2-13C]AcCoA
    cs2.set_all_isotopomers('AcCoA', [1.0, 0.0, 0.0, 0.0])#a mixture of 25% [2-13C]AcCoA and 25%[1,2-13C]AcCoA
    #
    # Isotope labelling information of carbon sources is set by specific isotopomers
    #
    cs.set_each_isotopomer('Asp', {'#0000':1.0}, correction = "no")
    cs2.set_each_isotopomer('Asp', {'#0000':1.0}, correction = "no")
    #
    # Generation of MDV instances from metabolic state and carbon source.
    #
    mdv = model.generate_mdv(flux, cs)
    data = mdv.get_fragment_mdv("Glue")
    print("Check MDV of Glu by the method demonstated in Antoniewiez et al 2007")
    print("Calced MDV data of              : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Answer in Antoniewiez et al 2007: [0.3464 0.2695 0.2708 0.0807 0.0286 0.0039]")
    mdv = model.generate_mdv(flux, cs2)
    data = mdv.get_fragment_mdv("Glue")
    print("Non labled carobon soureses without addition of isotope effects")
    print("Without natural isotope         : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Answer                          : [1.0000 0.0000 0.0000 0.0000 0.0000 0.0000]")
    mdv.add_natural_isotope()
    data = mdv.get_fragment_mdv("Glue")
    print("Addition of isotope effect by mdv.add_natural_isotope()")
    print("Correction by natural isotope   : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("C5H10N2O3                       : [0.9328 0.0594 0.0074 0.0004 0.0000 0.0000]")
    mdv.correct_natural_isotope()
    data = mdv.get_fragment_mdv("Glue")
    print("Cancellation of isotope effect by mdv.correct_natural_isotope()")
    print("Calced MDV data of              : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Return                          : [1.0000 0.0000 0.0000 0.0000 0.0000 0.0000]")
    #
    # Addition of isotope effect during MDV calculation
    #
    model.set_configuration(add_naturalisotope_in_calmdv = "yes") #
    #
    # Requires model.reconstruct()
    #
    model.reconstruct()
    #
    # Generate MDVs from metabolic state and carbon source.
    #
    mdv = model.generate_mdv(flux, cs2)
    data = mdv.get_fragment_mdv("Glue")
    print("Addition of isotope effect during MDV calculation by model.set_configuration(add_naturalisotope_in_calmdv = 'yes')")
    print("Without natural isotope         : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("C5H10N2O3                       : [0.9328 0.0594 0.0074 0.0004 0.0000 0.0000]")
    print("Cancellation of isotope effect by mdv.correct_natural_isotope()")
    mdv.correct_natural_isotope()
    data = mdv.get_fragment_mdv("Glue")
    print("Calced MDV data of              : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Return                          : [1.0000 0.0000 0.0000 0.0000 0.0000 0.0000]")
    mdv = model.generate_mdv(flux, cs2)
    data = mdv.get_fragment_mdv("Glugc")
    print("Addition of isotope effect during MDV calculation by model.set_configuration(add_naturalisotope_in_calmdv = 'yes')")
    print("Without natural isotope         : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("C5H10N2O3                       : [0.9328 0.0594 0.0074 0.0004 0.0000 0.0000]")
    print("Cancellation of isotope effect by mdv.correct_natural_isotope()")
    mdv.correct_natural_isotope()
    data = mdv.get_fragment_mdv("Glugc")
    print("Calced MDV data of              : [{0:5.4f} {1:5.4f} {2:5.4f} {3:5.4f} {4:5.4f} {5:5.4f}]".format(data[0], data[1], data[2], data[3], data[4], data[5]))
    print("Return                          : [1.0000 0.0000 0.0000 0.0000 0.0000 0.0000]")
