mfapy: A Python toolbox for 13C-metabolic flux analysis
================================================

mfapy is a Python toolbox for 13C-metabolic flux analysis developed by Matsuda and Shimizu lab group in Osaka university, Japan.

mfapy supports:

1. The elementary metabolite unit (EMU) framework developed by Antoniewicz et al 2007.
2. Isotopically nonstationary 13C metabolic flux analysis (INST-MFA)
3. Parallel labeling experiment without limitation of number of experiments.
4. Metabolic flux, net metabolic flux of reversible reaction, and metabolite concentration could be considered as fitting parameters.
5. G-value by introducing the "pseudo" reaction
6. Goodness-of-fit analysis
7. Parallel performing non-linear optimization jobs using joblib package.
8. Automated model construction from model description file

Documentation
----------------------------------------
Please check 'mfapy-document' repository. https://fumiomatsuda.github.io/mfapy-document/

License
----------------------------------------
This software is released under the MIT License, see LICENSE.txt.

Requirement
----------------------------------------
Python 3.6

mfapy requires the packages:

NumPy and SciPy (van der Walt, Colbert & Varoquaux, 2011)

nlopt https://nlopt.readthedocs.io/en/latest/

joblib

mkl-service

mfapy was developed and tested using the PyScripter IDE, and was tested in 64 bit version of Anaconda3 distribution in Windows 10.


Install
----------------------------------------
1. Install 64 bit version of Anaconda3
2. Create virtual environment such as named "mfapy" and install required packages
~~~
    > conda create -n mfapy python=3.9 numpy scipy matplotlib joblib
    > conda activate mfapy
    > conda install -c conda-forge nlopt
    > conda install -c anaconda mkl-service
~~~

2. Install mfapy
Unzip and set as current folder
~~~
    > python setup.py install
~~~

Test
----------------------------------------
~~~
    > python setup.py test
~~~

How to use
----------------------------------------
./sample/Tutorial 1_13C-MFAEcoli.py explains how to use mfapy for 13C-MFA.

./sample/Tutorial 2_INST-MFAtoymodel.py explains how to use mfapy for INST-MFA.

./sample/Tutorial 3_mfapy_functions.py is a advanced exlanation on many mfapy function.


Detailed explanations are available from docstrings

In the ./sample folder, example python codes are also available.


Learn more <https://fumiomatsuda.github.io/mfapy-document/>



Version 054
----------------------------------------
19/06/28  class MetabolicModel: search_ci Grid search procedure is updated.

19/06/28  class MdvData: add_gaussian_noise  "normalize" option is newly added.

19/07/01  optimize: def fit_r_mdv_deep add global optimization by "GN_CRS2_LM" before iteration

Version 055
----------------------------------------
19/12/27 Support Python 3.8

19/12/27 Example_1_toymodel_INST.py is updated to explain a method to construct time course mdv data from mdv files.

Version 057
----------------------------------------
20/7/12 initializing_Rm_fitting, fit_r_mdv_scipy, fit_r_mdv_nlopt in optimizaton: Expection is newly raised to avoid error in  paralell processing.

20/7/13 joblib instead of pp is employed for paralell proceccing.

20/7/30 Format of model definition file was updated to support external with a backward compatibility.

20/7/30 load_metabolic_model_reactions in mfapyio: Support external id.

20/7/30 load_metabolic_model_metabolites in mfapyio: Support external id.

20/7/30 load_metabolic_model_reversibles in mfapyio: Support external id.

20/7/30 load_metabolic_model_fragments in mfapyio: Support external id.

20/7/30 External id data was added to Example_0_toymodel_model.txt,

20/7/30 External id data was added to Example_1_toymodel_model.txt.

20/7/30 External id data was added to Example_2_MCF7_taxol_model.txt.

20/7/30 show_results in metablicmodel: Output format was modified for more beautiful alignment.

20/7/30 show_results in metablicmodel: "checkrss" option was added to check RSS levels of each fragment and "fitting" reactions and metabolites.

20/8/19 Files for Example 6 (the Metropolis-Hastings algorithm reported in Scientific Reports volume 10, Article number: 286 (2020)) were newly included in sample folder.

20/8/19 mfapy.carbonsource.set_carbonsources method is newly developed for batch loading of labelling pattern information of carbon sources.

20/8/31 Files for Explanation 1 were newly included in sample folder to demonstrate the functions of mfapy.carbonsource.set_carbonsources

20/8/31 mfapy.metabolicmodel.show_results_in_map is newly developed to project flux data on the metabolic map (.GML) available in Vanted. Files of Example 2 were updated to generate "Example_2_cancer_map_mapped.gml" from a blank map "Example_2_cancer_map_new.gml".

Version 058
----------------------------------------
20/9/3 Files for Example 7 (the Metropolis-Hastings algorithm reported in Scientific Reports volume 10, Article number: 286 (2020)) were updated to reproduce MCF-7 example.

20/9/3 A bug in mfapy.metabolicmodel.calc_rss was fixed.

20/9/6 test_metabolicmodel.py was updated to support mfapy.metabolicmodel.generarate_flux_distribution => mfapy.metabolicmodel.generarate_state

Version 059
----------------------------------------
20/12/24 Model check functions were added to MetabolicModel.

21/1/9 Example files were totally reorganized.

21/1/9 Explanation 1 files were newly added to explain mfapy functions.

Version 060
----------------------------------------
21/2/26 Support Google docstring format.

21/2/26 'mfapy-document' repository was created. https://fumiomatsuda.github.io/mfapy-document/

Version 061
----------------------------------------
22/2/14 Limiation in the maximum EMU size was removed.

22/11/7 Bugs in MetabolicModel.show_results was fixed.

Version 061
----------------------------------------
24/3/4  DocString contents were updated.

24/3/4  "Tutorial 1_13C-MFAEcoli.py" was created to expain how to use mfapy for 13C-MFA.

24/3/4  "Tutorial 2_INST-MFAtoymodel.py" was created to expain how to use mfapy for INST-MFA.

24/3/4  "Tutorial_3_mfapy_functions.py" was created to expain many mfapy functions.


If you want to learn more about ``setup.py`` files, check out `this repository <https://github.com/fumiomatsuda/mfapy>

