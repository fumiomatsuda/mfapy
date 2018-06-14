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
mfapy required the packages:
NumPy and SciPy (van der Walt, Colbert & Varoquaux, 2011)
nlopt https://nlopt.readthedocs.io/en/latest/
parallel python (pp-1.6.4.4) https://www.parallelpython.com/ developed by Vitalii Vanovschi.
mkl-service

mfapy was developed and tested using the PyScripter IDE, and was tested in 64 bit version of Anaconda3 distribution in Windows 10.


Install
----------------------------------------
1. Install 64 bit version of Anaconda3
2. Create virtual environment such as named "mfapy" and install required packages
~~~
    > conda create -n mfapy python=3.6 numpy scipy matplotlib=2.1
    > conda activate mfapy 
    > conda install -c conda-forge nlopt
    > conda install -c anaconda mkl-service
~~~
3. Download pp-1.6.4.4.zip from https://www.parallelpython.com/content/view/18/32/
Unzip and set as current folder
~~~
    > python setup.py install
~~~
3. Install mfapy
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

Example_1_toymodel.py  
Example_1_toymodel_INST.py 
Example_2_MCF7_taxol.py 

Learn more <https://github.com/fumiomatsuda/mfapy>

---------------

If you want to learn more about ``setup.py`` files, check out `this repository <https://github.com/fumiomatsuda/mfapy>



