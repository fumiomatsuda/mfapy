mfapy: A Python toolbox for 13C-metabolic flux analysis
================================================

mfapy is a Python toolbox for 13C-metabolic flux analysis developed by Matsuda and Shimizu lab group in Osaka university, Japan.
mfapy supports:

1. The elementary metabolite unit (EMU) framework developed by Antoniewicz et al 2007.
2. Isotopically nonstationary 13C metabolic flux analysis (INST-MFA)
3. Parallel labeling experiment without limitation of number of experiments.
4. Metabolic flux, net metabolic flux of reversible reaction, and metabolite concentration could be considered as fitting parameters.
5. G-index by introducing the "pseudo" reaction
6. Goodness-of-fit analysis
7. Parallel performing non-linear optimization jobs using parallel python.
8. Automated model construction from model description file
9. MS/MS spectra data

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
    > conda create -n mfapy python=3.6 numpy scipy matplotlib=2.1 joblib
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
Sample scripts in 'sample' explain how to used mfapy. 
Detailed explanations are available from docstrings

1. Example_1_toymodel.py  
2. Example_1_toymodel_INST.py 
3. Example_2_MCF7_taxol.py 
4. Example_3_Coryne.py 

Learn more <https://github.com/fumiomatsuda/mfapy>

Version 054
----------------------------------------
19/06/28  
class MetabolicModel: search_ci  
Grid search procedure is updated.  
19/06/28  
class MdvData: add_gaussian_noise  
"normalize" option is newly added.  
19/07/01  
optimize: def fit_r_mdv_deep  
add global optimization by "GN_CRS2_LM" before iteration  


Version 055
----------------------------------------
19/12/27  
Support Python 3.8
Example_1_toymodel_INST.py is updated to explain a method to construct time course mdv data from mdv files.  

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

20/8/19 Files for Example 4 (the Metropolis-Hastings algorithm reported in Scientific Reports volume 10, Article number: 286 (2020)) were newly included in sample folder.

20/8/19 mfapy.carbonsource.set_carbonsources method is newly developed for batch loading of labelling pattern information of carbon sources.

20/8/31 Files for Example 5 were newly included in sample folder to demonstrate the functions of mfapy.carbonsource.set_carbonsources

20/8/31 mfapy.metabolicmodel.show_results_in_map is newly developed to project flux data on the metabolic map (.GML) available in Vanted. Files of Example 2 were updated to generate "Example_2_cancer_map_mapped.gml" from a blank map "Example_2_cancer_map_new.gml".

Version 058
----------------------------------------
20/9/3 Files for Example 4 (the Metropolis-Hastings algorithm reported in Scientific Reports volume 10, Article number: 286 (2020)) were updated to reproduce MCF-7 example.

20/9/3 A bug in mfapy.metabolicmodel.calc_rss was fixed.

20/9/6 test_metabolicmodel.py was updated to support mfapy.metabolicmodel.generarate_flux_distribution => mfapy.metabolicmodel.generarate_state

If you want to learn more about ``setup.py`` files, check out `this repository <https://github.com/fumiomatsuda/mfapy>

