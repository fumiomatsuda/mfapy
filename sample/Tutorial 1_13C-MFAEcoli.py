#-------------------------------------------------------------------------------
# Name:        Tutiral 1 13C-MFA of metabolically engineered E. coli
#
#
# Author:      Fumio_Matsuda
#
# Created:     20/02/2024
# Copyright:   (c) Fumio_Matsuda 2024
# Licence:     MIT license
#-------------------------------------------------------------------------------
"""
############################################################
Part 1: Importing the mfapy Module
############################################################

mfapy is a Python toolbox for 13C-metabolic flux analysis.
This tutorial demonstrates the steps to perform 13C metabolic
flux analysis (13CMFA) of E. coli using the mfapy package.
The metabolic model and data used in this tutorial are derived
from Okahashi et al., Biotechnol. Bioeng., 2017, 114, 2782-2793.

The mfapy package encompasses all the necessary functions
for 13CMFA and can be downloaded from https://github.com/fumiomatsuda/mfapy.
Detailed installation instructions can be found on the GitHub page.
Additionally, comprehensive documentation for every function of
mfapy is available https://fumiomatsuda.github.io/mfapy-document/.

To begin, import the mfapy package:

"""

import mfapy

"""
############################################################
Part 2: Loading the Metabolic Model Definition File
############################################################

The mfapy package consists of five modules: mfapyio, metabolicmodel,
mdv, carbonsource, and optimize.
The optimize module provides low-level functions for various
optimization problems but is not directly used by users.

The metabolicmodel module defines a class called MetabolicModel.
An instance of this class (named "model" in this script) contains
a series of functions to perform 13C-MFA.


To generate a "model," information about an E. coli metabolic
model is loaded from a model definition file named
"Tutorial 1_13C-MFAEcoli_model.txt" using the mfapy.mfapyio.load_metabolic_model
function.
This information is stored in four dictionaries: reactions, reversible,
metabolites, and target_fragments.
These dictionaries contain information described in the "Metabolites,"
"Reversible Reactions," and "Target Fragments" sections of the model
definition file.
Users can modify the metabolic model by editing these dictionaries.

It's recommended to use "debug" mode if you encounter errors during the loading process.



"""

reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("Tutorial 1_13C-MFAEcoli_model.txt", format = "text", mode = "normal")

# mfapy.mfapyio.load_metabolic_model
# Function to load metabolic model information from a text file
# CAUTION: This function has limited error-checking functions.
# Args:
#     filename (str): Metabolic model file
#     format (str): "text" (default, tab-delimited) or "csv"
#     mode (str): "normal" (default) or "debug" (to show loaded metabolic file data)
# Returns:
#     - reactions (dict): Dictionary describing metabolite reactions
#     - reversible (dict): Dictionary defining reversible reactions
#     - metabolites (dict): Dictionary including metabolite information
#     - target_fragments (dict): Dictionary of target fragments


"""
############################################################
Model Definition File Structure:
############################################################

The model definition file, "Tutorial 1_13C-MFAEcoli_model.txt,"
comprises four sections identified by the labels: //Reactions,
//Metabolites, //Reversible_reactions, and //Target_fragments.
Each section ends with the //End label.

//Reactions Section:

The  //Reactions part in "Tutorial 1_13C-MFAEcoli_model.txt"
consists of 103 metabolic reactions in total.
The //Reactions section in "Tutorial 1_13C-MFAEcoli_model.txt"
contains 103 metabolic reactions.
Each line defines the stoichiometry and atom mapping of one reaction.
All reactions proceed in one direction from the left side to the right side.

Example:
//Reactions
# ID    For stoichiometry matrix    For atom mapping    Atom mapping    External ID Lower boundary    Upper boundary
r1    SubsGlc --> G6P    SubsGlc --> G6P    ABCDEF --> ABCDEF    (kegg:R02848)    0    1000
r2    G6P --> F6P    G6P --> F6P    ABCDEF --> ABCDEF    (kegg:R00771f)    0    1000
r3    F6P --> G6P    F6P --> G6P    ABCDEF --> ABCDEF    (kegg:R00771r)    0    1000

Explanation:
1. id, text
    Please use alphabet (+ number) as ids. Do not use space and
    special characters to avoid unexpected errors.
2. Stoichiometry for the construction of stoichiometry matrix, text
    All IDs of metabolites must be defined in //Metabolites.
    Only "+" and "-->" are acceptable separators.
    Here, expression of reaction using coefficient such as
    "Fumx --> {2}OA" is allowed.
    "nd" means this reaction is a "pseudo" reaction and ignored
    in the flux calculation of substrate by stoichiometry matrix.
3. Stoichiometry for the construction of atom mapping, text
    Only "+" and "-->" are acceptable separators.
    Here, expression of reaction without using coefficient is
    needed such as  "Fumx --> OA+OA".
    "nd" means this reaction is an excretion reaction and ignored
    in the MDV calculation
4. Atom mapping, text
    Only "+" and "-->" are acceptable separators. Please use A-Z to
    express carbon atoms.
    "nd" means this reaction is an excretion reaction and is ignored
    in the MDV calculation
5. Links to external IDs (Optional)
    Links to external IDs can be described here
6. Lower boundary of metabolic flux level, float (Optional)
    Value must exceed zero.
7. Upper boundary of metabolic flux level, float (Optional)
    Value must exceed zero.


//Reversible_reactions Section:

The //Reversible_reactions section in "Tutorial 1_13C-MFAEcoli_model.txt"
defines 13 reversible reactions.
Each line defines a pair of reaction IDs from the //Reactions section.

Reversible reactions can be used for various purposes such as to output
metabolic flux distribution data as net flux values and to constrain
metabolic flux.


Example:
//Reversible_reactions
#Name    Forward reaction ID    Reverse reaction ID    External ID (Experimental)    Lower boundary    Upper boundary
PGI    r2    r3    (kegg:R00771)    -1000    1000
FBA    r4    r5    (kegg:R01068)    -1000    1000


Explanations:
1. Id, text
    Use alphanumeric characters, avoiding spaces and special characters
    to avoid unexpected errors.
2. Dorward reaction, text
    Reaction ID
    Composition of reaction IDs separated by "+" such as "v9+e1" is also supported.
3. Reverse reaction, text
    Reacton id
    Use "nothing" if there is no reversible data.
4. Links to external IDs (Optional)
    Links to external IDs can be described here
5. Lower boundary of metabolic flux level, float (Optional)
    Value must be less than zero
6. Upper boundary of metabolic flux level, float (Optional)
    Value must exceed zero.


/Metabolites Section:

The //Metabolites section in "Tutorial 1_13C-MFAEcoli_model.txt" defines
69 metabolites.
Each line specifies the properties of one metabolite, including the number
of carbon atoms and the type of metabolite.



Example:
//Metabolites
# ID    Number of atom    Symmetry    Carbon source    Excreted metabolite    External ID (Experimental)    Lower boundary    Upper boundary
m6PG    6    no    no    no    (kegg:C00345)    0    300
SubsGlc    6    no    carbonsource    no    (dummy)    0    300
Suc    4    symmetry    no    no    (kegg:C00042)    0    300
CO2ex    1    no    no    excreted    (dummy)    0

Explanations:
1. ID, text
    Please use alphanumeric characters. Avoid using spaces, underscores,
    and other special characters to prevent unexpected errors.
    The initial character must be an alphabet. The id "6PG" was not allowed.
    Please convert to "m6PG" or other names
2. Number of atom, int
    Max number of carbon in "carbonsource" metabolites are 8.
3. Symmetry, text
    Set "symmetry" metabolites such as fumarate and succinate
4. Carbon source, text
    Set "carbonsource" metabolites
5. Excretion, text
    Set "excreted" metabolites
6. Links to external IDs (Optional)
    Links to external IDs can be described here
7. Lower boundary of metabolite pool size, float (Optional)
    This is required for INST-MFA mode and must be greater than zero.
8. Upper boundary of metabolite pool size, float (Optional)
    This is required for INST-MFA mode and must be greater than zero.

//Target_fragments Section:

The //Target_fragments section in "Tutorial 1_13C-MFAEcoli_model.txt"
defines 24 fragments.
In 13C-MFA, the term "fragments" refers to ions observed with a mass
spectrometer.
This is because fragments generated by GC-EI-MS are widely used to measure
mass distribution vectors (MDV) of amino acid derivatives.
Molecular-related ions observed by LC-ESI-MS are also applicable for 13C-MFA.

The example below illustrates a fragment, "Ala[M-57]," which includes
carbons from the Ala_1:2:3 (1st, 2nd, and 3rd) carbons of Ala.
Here, "Ala[M-57]" is the ID of a fragment, and "Ala" is the name of a
metabolite described in the //Metabolites section. Each line defines one fragment.

Example:
//Target_fragments
#ID     Type    Corresponding metabolite and atoms Usage     Formula (Experimental)
Ala[M-57]    gcms    Ala_1:2:3    use    C11H26N2O1Si2
Ala[M-85]    gcms    Ala_2:3    use    C10H26N1O1Si2
Asp[M-57]    gcms    Asp_1:2:3:4    use    C18H40N4O1Si3

Explanations:
1. ID, text
    Please use alphanumeric characters. Avoid using spaces, underscores,
    and other special characters to prevent unexpected errors.
2. Type, text
    Only "gcms" is available. Please use "gcms" even if using LCMS or CEMS data.
3. Corresponding metabolite and atoms, text
    "Glu_1:2:3:4:5" indicates that the target compound has five carbons
    derived from the 1st to 5th carbons of metabolite ID "Glu."
In "gcms" mode, "OAC_1:2+OAC_1:2" indicates that the target compound has
    four carbons derived from the 1st and 2nd carbons of metabolite ID "OAC."
4. Usage, text
    "use" indicates that this target fragment is used in the analysis.
5. Formula, text (experimental)
    Chemical formula of the target compound (Experimental).


############################################################
Part 3: Generation of an Instance of the mfapy.metabolicmodel.MetabolicModel Class
############################################################

An instance of the mfapy.metabolicmodel.MetabolicModel class,
named "model," is generated and initialized using information
from the four dictionaries: reactions, reversible, metabolites,
and target_fragments.

During initialization, a stoichiometry matrix and a function
are produced.
The matrix of the metabolic network allows the generation of
a metabolic flux distribution from an independent flux vector.
Additionally, the function allows the calculation of the predicted
Mass Distribution Vector (MDV) of fragments from a metabolic flux
distribution, using isotope labeling information from carbon sources.

Consider using "debug" mode for additional information if you
encounter errors during this process.


"""

model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments, mode = "normal")

#mfapy.metabolicmodel.MetabolicModel
#Generator of new instance
#       Args:
#            reactions (dict): Dictionary describing metabolite reactions
#           reversible (dict): Dictionary for defining reversible reactions
#            metabolites (dict): Dictionary including metabolite information
#            target_fragments (dict): Dictionary of target_fragments
#
#            mode (str): "normal" or "debug" (to show loaded metabolic file data)


"""
############################################################
Part 4.  Addition of constraints and structure of status file
############################################################

The metabolic model variables, including flux levels of metabolic
reactions and net flux levels of reversible reactions, can be constrained
by measured values in mfapy. The metabolite concentrations are also constrained in the case of INST-MFA as additional variables.

There are two types of constraints:

* Lower and upper boundaries of metabolic flux/metabolite concentration
 (INST-MFA only) levels.
* Measured metabolic flux/metabolite concentration values with their
standard deviations (stdev).

All reactions and reversible reactions require lower and upper boundaries.
Reactions are classified into four types:
* free: Variable within lower and upper boundaries. The "value" and "stdev"
are ignored.
* fixed: Fixed to "value". The "stdev" is ignored.
* fitting: Variable within lower and upper boundaries. Also included in the
Residual Sum of Squares (RSS) calculation using "value" and "stdev".
* pseudo: Fixed to "value". Used for Mass Distribution Vector (MDV)
calculations but ignored in the stoichiometry matrix. The "stdev" is
ignored. Only used for Reactions.
In this tutorial, all constraint information is described in the status
file, Tutorial 1_13C-MFAEcoli_status.csv.


Example:
State    ID    type    value    std    lb    ub
reaction    r1    fitting    100    1.59    0.0001    1000
reaction    r2    free    1    1    0.0001    1000

Explanation:
State: Type of entities. "reaction", "metabolite", and "reversible" are available
Id: IDs of reactions, metabolites, and reversible reactions.
    IDs must be unique and identical to those in the model definition file.
type: Type of constraints including free, fitting, fixed, and pseudo
value: Measured values of metabolic flux and metabolite concentration.
Used for "fixed" and "fitting" reactions/metabolites, and “pseudo” reaction
stdev: Standard deviation of metabolic flux and metabolite concentration measurement.
       Used only for "fitting" reactions/metabolites.
lb: lower boundary of metabolic flux and metabolite concentration
ub: upper boundary of metabolic flux and metabolite concentration

For most reactions in "Tutorial 1_13C-MFAEcoli_status.csv", lower
and upper boundaries are set at 0.0001 and 10000, respectively.
Lower and upper boundary levels of -300 and 300 are employed for
reversible reactions.

In this example, reactions for substrate uptake (e.g., r1) and
product excretion (e.g., r47, r48, and r50) have "value" and
"stdev" values.
These are specific rates calculated from measured time-course data
of medium component analysis. They are constrained as "fitting" types.


reaction    r1    fitting    100    1.59    0.0001    1000    (SubsGlc --> G6P)
reaction    r47    fitting    1.48    0.96    0.0001    1000    (Acetone --> Acetoneex)
reaction    r48    fitting    13.34    1.50    0    1000    (IPA --> IPAex)
reaction    r50    fitting    36.81    2.81    0.0001    1000    (Acetate --> Acetateex)

Artificial reactions for biomass synthesis from precursors are
constrained as "fixed" types:

reaction    r72    fixed    0.95    1    0.0001    1000    G6P --> Biomass
reaction    r73    fixed    0.32    1    0.0001    1000    F6P --> Biomass

In this example, all metabolite concentrations are constrained
using dummy information (e.g., 10) because these values were not
used in 13C-MFA. Measured values are required for INST-MFA.

metabolite    m6PG    fixed    10    1    0.001    100
metabolite    AcAc    fixed    10    1    0.001    100

The "pseudo" type is used to consider the carry-over natural isotopes
in the MDV data. In this 13C-MFA analysis, proteinous Ala is a target
fragment for the MDV measurement. It includes a small and unknown
part of non-labeled Ala derived from the protein fraction in the
starter cell biomass.

In E. coli, Ala is synthesized from Pyr. The amount of Pyr used
for Ala synthesis is included in r80.

reaction    r80    fixed    13.13    1    0.0001    1000    Pyr --> Biomass

The "Tutorial 1_13C-MFAEcoli_status.csv" file defines r58 as
a "pseudo" reaction:

reaction    r58    pseudo    1    1    0.0001    1000    Pyr --> Ala
reaction    r95    free    1    1    0.0001    1000    SubsAla --> Ala

It means that r58 is ignored in the stoichiometry matrix and the
calculation of metabolic flux distribution.
On the other hand, the calculation of MDV uses r58 considering
that 13C-labeled Ala is virtually synthesized from 13C-labled Pyr
at a fixed rate of 1.0.
r95 indicates that non-labeled Ala is synthesized from non-labeled
SubsAla.
Since r95 is a “free” reaction, a mixing rate of 13C-labeled Ala
and non-labeled Ala can be determined by fitting the measured MDV
of fragments of Ala using the non-linear optimization.


The ’model’ instance has a method, ’load_states’, to read a status
file such as "Tutorial 1_13C-MFAEcoli_status.csv" file, generating a
dictionary including the constraint information.
The ’set_constraints_from_state_dict’ method sets the constraints
all at once.
After modification of constraints, the "model" instance needs to be
updated by the "update" method.
The stoichiometry matrix of a metabolic network and the function to
calculate a predicted MDV of fragments are reconstructed by this procedure.

"""

state_dic = model.load_states("Tutorial 1_13C-MFAEcoli_status.csv", format = 'csv')
    #Load a text/csv file with 'reacton type' information to generate new states dict.
    # Args:
    #     filename (str): filename of flux data with following format::
    #   format: 'csv' CSV or 'text'tab-deliminated text.
    #Returns:
    #   dict: Dictionary of states dict.

model.set_constraints_from_state_dict(state_dic)
    # Reaction types, value, stdev, lb, and ub are set from the state dict data at once
    #Args:
    #    dict (dict): Dictionary of flux. Data in 'value', 'stdev', 'lb', 'ub' and 'type' fields are used.


model.update()
    #Method to generate stoichiometry matrix.
    #A metabolic model has to be updated when any types (fixed, free, fitting, pseudo) of reactions, metabolites, and reversible reactions are modified.
    #Args:
    #    Not required.

"""
############################################################
Part 5: Configuration of Isotope Labeling of Carbon Sources
############################################################

In mfapy, the CarbonSource class in the ‘carbonsource’ module
handles all information regarding the Isotope labeling of carbon sources.
The "model" instance has a method, generate_carbon_source_template,
to create an instance of the CarbonSource class.
In this example, the instance named "carbon_source1" is generated.

"""

carbon_source1 = model.generate_carbon_source_template()

"""
The instance "carbon_source1" contains information on metabolites
defined as "carbonsource" in the model definition file.
In this example, there are 11 carbon sources in total.
Initially, all carbon sources are set to a non-labeled form.
The "carbon_source1" instance has the method set_each_isotopomer
 to configure the isotope labeling.

The set_each_isotopomer method accepts three arguments:

1) An ID of the carbon source metabolite.
    The ID must be defined in the model definition file.
2) A dictionary of isotopomers and their ratios.
    For example, for SubsGlc with six carbons, the keys '#000000',
    '#100000', and '#111111' represent non-labeled glucose, [1-13C]glucose,
    and [U-13C]glucose, respectively.
    The values (e.g., 0.02, 0.7, and 0.28) denote their relative abundances.
3) Correction by natural isotope.
    If correction = 'yes', the composition of the isotopomer is automatically
    adjusted considering the occurrence of natural 13C.

"""

carbon_source1.set_each_isotopomer('SubsGlc',{'#000000': 0.02, '#100000': 0.7, '#111111': 0.28 }, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsCO2',{'#0': 1.0, '#1': 0.0}, correction = 'no')
carbon_source1.set_each_isotopomer('SubsAla',{'#000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsAsp',{'#0000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsThr',{'#0000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsGlu',{'#00000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsPhe',{'#000000000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsTyr',{'#000000000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsIle',{'#000000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsLeu',{'#000000': 1.0}, correction = 'yes')
carbon_source1.set_each_isotopomer('SubsVal',{'#00000': 1.0}, correction = 'yes')

"""
############################################################
Part 6: Loading of measured mass isotopomer distribution vector
(MDV) data.
############################################################

In mfapy, the MdvData class in the ‘mdv’ module offers a series
of functions to handle MDV (Mass Distribution Vector) data.
The "model" instance provides a method, load_mdv_data, to read
an MDVfile and generate an instance of the MdvData class.
In this example, the instance named "mdv_observed1" is created
from 'Tutorial 1_13C-MFAEcoli_mdv.txt'.

The MDV file contains MDV data of measured fragments with the
following format

Explanation.
1 Name: ID of the target fragment. The ID must be specified in
the model description file.
2 Spectrum: The values 0, 1, ... correspond to m0, m1, ... signals
 in the mass spec data of the target fragment.
3 Selection: 0 or 1.  A data point with 0/1 is used/ignored for
 the calculation of MDV.
4 MDV: Relative abundance of mass signals
5 Std: Standard deviation level of MDV measurement

Example:
Name    Spectrum    Selection    MDV    Std
Ala[M-57]    0    1    0.31367837    0.01
Ala[M-57]    1    1    0.400545455    0.01
Ala[M-57]    2    1    0.024866841    0.01
Ala[M-57]    3    1    0.260909333    0.01
Ala[M-85]    0    1    0.698099329    0.01
Ala[M-85]    1    1    0.033802571    0.01
Ala[M-85]    2    1    0.2680981    0.01
Asp[M-57]    0    1    0.374679462    0.01


"""

mdv_observed1 = model.load_mdv_data('Tutorial 1_13C-MFAEcoli_mdv.txt')
       #def load_mdv_data(self, filename, format = 'text',output = "normal"):
       #Load MDV data from the text file. This function generates a new instance of mfapy.mdv class from an instance of mfapy.model.
       # Args:
       #     filename (str) : filename of MDV data file
       #     format (str) : "text" (defalut) or "csv"
       #     output (str) : "normal" (defalut) or "debug"
        #Returns:
        #    Boolean: True/False


"""

All data points of MDV must have a standard deviation (std) level
of measurement to calculate the residual sum of squares.
For 13C-MFA of microbes using GC-MS data, the std level of MDV
measurement is typically set to ‘0.01’.
The mdv_observed1 instance provides a method set_std by which std
values of all data points can be set at once.

"""

mdv_observed1.set_std(0.01, method = 'absolute')
        #set_std(value, method):
        # Setter of standard deviation level
        #Args:
        #    value (float): level of standard deviation.
        #    method (str): Method to calculate stdev levels.
        #        * 'relative': Levels are set by stdev = [value] * signal intensity.
        #        * 'absolute' (detault) : Levels are set by stdev = [value].

"""
############################################################
Part 7: Registration of isotope labeling experiment
############################################################

You can utilize the set_experiment() method to register a labeling
experiment.
Here, a pair of carbon source information (‘carbon_source1’)and
MDV data (‘mdv_observed1’) is registered to the ‘model’ as a labeling
experiment named 'ex1'.

"""
model.set_experiment('ex1', mdv_observed1, carbon_source1)
        #def set_experiment(self, name, mdv, carbon_sources, startidv = []):
        #Setter of an 'experiment' to metabolic model.
        #Here an 'experiment' indicated a set of a carbon source labeling pattern and a measured MDV data.
        #Parallel labeling experiment can be performed by setting multiple sets of 'experiment'.
        #Args:
        #    name (str): name of the experiment (unique)
        #    mdv (instance): MDVdata instance including measured MDV data of target fragments.
         #   carbon_sources (instance): Instance of carbon source object
         #   startidv (array): array of idv as starting isotope distribution for INST-MFA

"""
############################################################
Part 8: Configuration of model parameters.
############################################################

The parameters of the “model” instance are set using the set_configuration method.
The parameters ”callbacklevel”, “ncpus”, and “iteration_max” are often configured by users,

            'callbacklevel': default = 0,
                'callbacklevel' determines a frequency level of callbacks from the metabolic model.
                # 0: No callback
                # 1: Error report
                # 2: Simple progress report for grid search
                # 3: Detailed progress report for grid search
                # 4: Simple progress report for flux fitting
                # 5: Detailed progress report for flux fitting
                # 7: Detailed report for model construction
                # 8: Debug
            'ncpus':default = 3,
                Number of CPUs used in parallel processing
            'iteration_max': default = 1000,
                Maximal number of iterations (steps) in each fitting task.
All parameters are listed and explained on the mfapy documentation page.

"""

model.set_configuration(callbacklevel = 0) # A frequency level of callbacks
model.set_configuration(ncpus = 4) #Number of CPUs used in parallel processing
model.set_configuration(iteration_max = 10000) #Maximal number of iterations.

"""
############################################################
Part 9: Generation of Initial Flux Distribution
############################################################

In mfapy, non-linear optimization is used to find a metabolic
state (flux distribution) that aligns best with the measured MDV data.
A metabolic "state" includes a metabolic flux distribution (13C-MFA).
For INST-MFA a metabolite concentration profile is also included.
The functions for non-linear optimization necessitate an initial
or seed state. To prevent local optima, the non-linear optimization
task should be performed iteratively using randomly generated initial states.

The method generate_initial_states is used to generate multiple
initial states.
In this tutorial, 10 initial states are generated and stored in
“flux_opt1” as a list of 10 initial states.

"""

print("Generating 10 random initial states using parallel processing")
state, flux_opt1 = model.generate_initial_states(100, 10, method ="parallel")
    #def generate_initial_states(self, iterations = 100, initial_states = 1, method = 'normal'):
    # Generator of random initial metabolic states for model fitting
    # This method randomly generates initial metabolic states for a specified number of iterations. From these iterations, better states with lower RSS (Residual Sum of Squares) are selected.
    # Args:
    #     iterations (int): Number of trials to generate initial metabolic state. A trial may fail due to a bad optimization state.
    #     initial_states (int): Number of initial metabolic state(s) with lower RSS to be generated.
    #     method (str): Computational method
    #         * 'normal' :  Generation of flux distributions by a single thread.
    #         * 'parallel' : Parallel generation of flux distributions using parallel processing.
    #  Returns:
    #     * state  (str or array) : stop condition
    #         * [number_of_initial_states > 1 ] Str of stop condition.
    #         * [number_of_initial_states == 1 ] Array of Str of stop condition
    #     * flux (array or dict) : Dictionary of the metabolic state including metabolic flux and metabolite pool size levels.
    #        * [number_of_initial_states > 1 ] Array of dictionary of metabolic state including metabolic flux and metabolite pool size levels.
    #        * [number_of_initial_states == 1 ] Dictionary of the metabolic state including metabolic flux and metabolite pool size levels.



"""
############################################################
Part 10: Model Fitting
############################################################


“fitting_flux” is the method to perform non-linear optimization
using initial states.
mfapy employs NLopt and SciPy as a library of nonlinear optimization
methods. In this example, two nonlinear optimization methods,
"GN_CRS2_LM (a global optimizer)" and SLSQP (a local optimizer)"
are iteratively used to avoid local optima.
The “fitting_flux” method uses the joblib package for parallel
computing If a list of multiple initial states is given.
In this example,  the “fitting_flux” method returns a list of
optimized states, flux_opt1.
The best-fitted state with the smallest RSS is available as “
flux_opt1[0]” because the list is ordered by RSS values.
The goodness_of_fit method is used to calculate the p-value of
fitting and a threshold level of RSS.



"""
print("Fitting by non-linear optimizations")
for method in ["GN_CRS2_LM", "SLSQP","GN_CRS2_LM", "SLSQP"]:
    state, RSS_bestfit, flux_opt1 = model.fitting_flux(method = method, flux = flux_opt1)
            # def fitting_flux(self, method = 'SLSQP', flux = []):
            # Method for fitting a metabolic model to multiple experiments
            # Args:
            #     method (str): Algorism used for optimization.
            #         * 'SLSQP': Sequential Least SQuares Programming algorithm implemented in opt
            #         * "GN_CRS2_LM": Controlled random search implemented for global optimization by nlopt
            #      flux (dict): Initial flux distribution generated by self.generate_initial_states():
            #          * When flux is a dict of one state, single fitting trial is executed.
            #          * When flux is an array of multiple dicts of states, multiple fitting trials are executed by using parallel processing.
            #  Returns:
            #     state  (str or array) : stop condition
            #          * [number_of_initial_states == 1 ] Str of stop condition.
            #          * [number_of_initial_states > 1 ] Array of str of the stop condition.
            #     flux_opt (array or dict): Initial flux distribution generated by self.generate_initial_states()
            #          * [number_of_initial_states == 1 ] dict of optimized metabolic state.
            #          * [number_of_initial_states > 1 ] Array of dict of optimized metabolic state. Sorted by ascending order of their RSS levels
            #     RSS_bestfit (array or float) : Dictionary of the metabolic state including metabolic flux and metabolite pool size levels.
            #          * [number_of_initial_states == 1 ] Dictionary of the metabolic state including metabolic flux and metabolite pool size levels.
            #          * [number_of_initial_states > 1 ] Array of dictionary of the metabolic state including metabolic flux and metabolite pool size levels.

    pvalue, rss_thres = model.goodness_of_fit(flux_opt1[0], alpha = 0.05)

            #def goodness_of_fit(self, flux, alpha = 0.05):
            # Getter to calculate goodness-of-fit of a given flux distribution

            #   Args:
            #       flux (dict): Dictionary of metabolic flux distribution
            #       alpha (float): Confidence level. Default is 0.05.
            #   Returns:
            #       * pvalue (float) p-value determined by chi-squere test.
            #       * rss_thres (float) Threshold RSS at the alpha levels.

    print("Method; {}, p-value: {:.1f}, RSS_bestfit: {:.1f}, RSS_thres: {:.1f}".format(method, pvalue, RSS_bestfit[0], rss_thres))


"""
############################################################
Part 11: Output
############################################################

The method show_results facilitates the display of one or more fitting results.
One or more metabolic states are provided as a list of tuples, each containing
pairs of "name" and a dictionary of metabolic states, like so:

[('name1', state1), ('name2', state2), ('name3', state3)]

By default, show_results displays the results on the console. However, it can
also generate a file if the filename argument is specified

"""
results = [(method, flux_opt1[0])]
model.show_results(results)
model.show_results(results, filename = "Tutorial 1_13C-MFAEcoli_output.csv", format = "csv")
       #Method to output metabolic state data.
       # Args:
       #     input (tuple): Tuple of ("name", dic of metabolic state) [('name1', state1),('name2', state2),('name3', state3),]
       #     filename (str): Results are saved when the file name is not ""
       #     format (str): "csv" or "text"



