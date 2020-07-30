#-------------------------------------------------------------------------------
# Name:        mfapyio.py
# Purpose:     Input/Output funtions in mfapy
#
# Author:      Fumio_Matsuda
#
# Created:     12/06/2018
# Copyright:   (c) Fumio_Matsuda 2018
# Licence:     MIT license
#-------------------------------------------------------------------------------

import numpy as numpy
from . import mdv
from . import metabolicmodel
import re, csv

def load_metabolic_model_reactions(filename, format = "text", output = "normal"):
    """
    Load metabolic reaction information from a text or CSV file with following format.

    //Reactions
    v1	AcCoA + OAC --> Cit	AcCoA + OAC --> Cit	AB + CDEF --> FEDBAC	(kegg:R00351)	0.0	300
    v2	Cit --> AKG + CO2ex	Cit --> AKG + CO2ex	ABCDEF --> ABCDE + F	(kegg:R00709)	0.0	300
    v3	AKG --> Glu	AKG --> Glu	ABCDE --> ABCDE 	(kegg:R00243)	0.0	300
    v4	AKG --> Suc + CO2ex	AKG --> Suc + CO2ex	ABCDE --> BCDE + A	(kegg:R01197)	0.0	300
    v5	Suc --> Fum	Suc --> Fum	ABCD --> ABCD	(kegg:R02164)	0.0	300
    v6	Fum --> OAC	Fum --> OAC	ABCD --> ABCD	(kegg:R01082)	0.0	300
    v7	OAC --> Fum	OAC --> Fum	ABCD --> ABCD	(kegg:R01082)	0.0	300
    v8	Asp --> OAC	Asp --> OAC	ABCD --> ABCD	(kegg:R00355)	0.0	300

    Parameters
    ----------
    filename : metabolic model file
    format : "text" (defalut, tab deliminated) or "csv"
    output : "normal" (defalut) or "debug" (to show loaded metabolic file data)

    Returns
    -------
    reactions : Dictionary describing metabolite reactions


    Examples
    --------
    >>> reactions = load_metabolic_model_reaction('filename.txt')
    >>> print reactions
    {'r4': {'atommap': 'AB+CDEF-->FEDBAC', 'reaction': 'ACCOA+OAA-->ICI', 'use': 'use', 'lb': 0.1, 'flux_value': 0.0, 'reversible': 'no', 'flux_var': 1.0, 'type': 'free', 'order': 3, 'ub': 300.0},...}

    """
    #
    # Distionary for reaction data
    #
    reactions = {}
    #
    # Conter to keep reaction order
    #
    counter = 0
    #
    # Ititialize mode
    #
    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            if len(row) == 0:
                continue
            if "#" == row[0][0]:
                continue
            row = [item for item in row if "#" not in item]

            if len(row) < 1:
                continue

            if "//" in row[0]:
                rid = row[0].replace(" ", "")
                if rid == "//Reactions":
                    mode = "Reactions"
                    continue
                if rid == "//End":
                    break
                mode = "other"
            # Remove '' from row
            row = list(filter(lambda s:s != '', row))


            if mode == "Reactions":
                if output == "debug":
                    print(row)
                if len(row) < 6:
                    continue
                rid = row[0].replace(" ", "")
                stoichiometry = row[1].replace(" ", "")
                reaction = row[2].replace(" ", "")
                atommap = row[3].replace(" ", "")
                exid = row[4]
                if rid == "":
                    continue


                reactions[rid] = {
                'stoichiometry': stoichiometry,
                'reaction':reaction,
                'atommap':atommap,
                'externalids' :exid,
                'order':int(counter)
                }
                if len(row) >= 6:
                    reactions[rid]['lb'] = float(row[5].replace(" ", ""))
                    reactions[rid]['ub'] = float(row[6].replace(" ", ""))
                counter = counter + 1

    return(reactions)

def load_metabolic_model_metabolites(filename, format = "text", output = "normal"):
    """
    Load Metabolite information from a file with following format.

    //Metabolites
    CO2ex	1	no	no	excreted	(kegg:C00011)	0.0	300
    AcCoA	2	no	carbonsource	no	(kegg:C00024)	0.0	300
    OAC	4	no	no	no	(kegg:C00036)	0.0	300
    Cit	6	no	no	no	(kegg:C00158)	0.0	300
    AKG	5	no	no	no	(kegg:C00026)	0.0	300
    Suc	4	symmetry	no	no	(kegg:C00042)	0.0	300
    Fum	4	symmetry	no	no	(kegg:C00122)	0.0	300
    Glu	5	no	no	no	(kegg:C00025)	0.0	300
    Asp	4	no	carbonsource	no	(kegg:C00049)	0.0	300

    Parameters
    ----------
    filename : metabolic model file
    format : "text" (defalut, tab deliminated) or "csv"
    output : "normal" (defalut) or "debug" (to show loaded metabolic file data)

    Returns
    -------
    metabolites : Dictionary including metabolite information

    Examples
    --------
    >>> metabolites = load_metabolic_model_metabolites('filename.txt')
    >>> print metabolites
    {'GLUEX': {'excreted': 'no', 'carbonsource': 'carbonsource', 'C_number': 5, 'symmetry': 'no', 'order': 4}, ...}

    """
    metabolites = {}
    counter = 0

    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            if len(row) == 0:
                continue
            if "#" == row[0][0]:
                continue
            row = [item for item in row if "#" not in item]
            if len(row) < 1:
                continue

            if "//" in row[0]:
                if row[0] == "//Metabolites":
                    mode = "Metabolites"
                    continue
                if row[0] == "//End":
                    break
                mode = "other"
            # Remove '' from row
            row = list(filter(lambda s:s != '', row))




            #
            # ?????????
            #
            if mode == "Metabolites":
                if output == "debug":
                    print(row)
                if len(row) < 5:
                    continue
                name = row[0].replace(" ", "")
                C_number = row[1]
                symmetry = row[2].replace(" ", "")
                carbonsource = row[3].replace(" ", "")
                excreted = row[4].replace(" ", "")
                if name=="": continue
                metabolites[name] = {
                'C_number' :int(C_number),
                'symmetry':symmetry,
                'carbonsource':carbonsource,
                'excreted':excreted,
                'order' :int(counter),
                'externalids' :"no external id",
                }
                if len(row) >= 6:
                    metabolites[name]['externalids'] = row[5]
                if len(row) >= 7:
                    metabolites[name]['lb'] = float(row[6].replace(" ", ""))
                    metabolites[name]['ub'] = float(row[7].replace(" ", ""))
                counter = counter + 1



                counter = counter + 1
    return(metabolites)

def load_metabolic_model_reversibles(filename, format = "text", output = "normal"):
    """
    Load definitions of reversible reactions from a metabolic model file with following format.

    //Reversible_reactions
    FUM	v6	v7	(kegg:R01082)	0.0	300

    Parameters
    ----------
    filename : metabolic model file
    format : "text" (defalut, tab deliminated) or "csv"
    output : "normal" (defalut) or "debug" (to show loaded metabolic file data)

    Returns
    -------
    reversible : Dictionary for defining reversible reactions

    Examples
    --------
    >>> reversible = load_metabolic_model_reversibles('filename.txt')
    >>> print reversible
    {'MDH': {'flux_value': 0.0, 'reverse': 'r27', 'flux_var': 1.0, 'forward': 'r26', 'type': 'free', 'order': 6}, ...}

    """
    dic = {}
    counter = 0

    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            if len(row) == 0:
                continue
            if "#" == row[0][0]:
                continue
            row = [item for item in row if "#" not in item]

            if len(row) < 1:
                continue

            if "//" in row[0]:
                if row[0] == "//Reversible_reactions":
                    mode = "Reversible_reactions"
                    continue
                if row[0] == "//End":
                    break
                mode = "other"
            # Remove '' from row
            row = list(filter(lambda s:s != '', row))



            if mode == "Reversible_reactions":
                if output == "debug":
                    print(row)
                if len(row) < 3:
                    continue
                name = row[0].replace(" ", "")
                forward = row[1].replace(" ", "")
                reverse = row[2].replace(" ", "")
                #rtype = row[3].replace(" ", "")
                dic[name] = {
                'forward':forward,
                'reverse':reverse,
                'type': "free",
                'order' :int(counter),
                'externalids' :"no external id",
                }
                if len(row) >= 4:
                    dic[name]['externalids'] = row[3]
                if len(row) >= 5:
                    dic[name]['lb'] = float(row[4].replace(" ", ""))
                    dic[name]['ub'] = float(row[5].replace(" ", ""))
                counter = counter + 1

    return(dic)

def load_metabolic_model_fragments(filename, format = "text", output = "normal"):
    """
    Load mass fragment information from a metabolic model file with following format.

    //Target_fragments
    Glue	gcms	Glu_1:2:3:4:5	use	C5H10N2O3
    Gluee	gcms	Glu_1:2:3+Glu_4:5	use	C5H10N2O3

    Parameters
    ----------
    filename : metabolic model file
    format : "text" (defalut, tab deliminated) or "csv"
    output : "normal" (defalut) or "debug" (to show loaded metabolic file data)

    Returns
    -------
    target_fragments : Dictionary of target_fragments


    Examples
    --------
    >>> target_fragments = load_metabolic_model_fragments('filename.txt')
    >>> print target_fragments
    {'Thr302': {'atommap': 'Thr_12', 'use': 'no', 'type': 'gcms', 'order': 35, 'number': 3}, ...}

    """
    dic = {}
    counter = 0

    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            if len(row) == 0:
                continue
            if "#" == row[0][0]:
                continue
            row = [item for item in row if "#" not in item]

            if len(row) < 1:
                continue

            if "//" in row[0]:
                if row[0] == "//Target_fragments":
                    mode = "Target_fragments"
                    continue
                if row[0] == "//End":
                    break
                mode = "other"
            row = list(filter(lambda s:s != '', row))


            if mode == "Target_fragments":
                if output == "debug":
                    print(row)
                if len(row) < 4:
                    continue
                formula = ""
                name = row[0].replace(" ", "")
                mtype = row[1].replace(" ", "")
                atommap = row[2].replace(" ", "")
                use = row[3].replace(" ", "")
                if len(row) >= 5:
                    formula = row[4].replace(" ", "")
                if name=="": continue
                dic[name] = {
                'type':mtype,
                'atommap':atommap,
                'use':use,
                'order' :int(counter),
                'formula':formula
                }
                counter = counter + 1

    return(dic)

def load_metabolic_model(filename, format = "text",output = "normal"):
    """
    Read metabolic model information from a text file with following format.

    CAUTION: Now this function has no error checking.

    # Example metabolite model
    //Reactions
    v1	AcCoA + OAC --> Cit	AB + CDEF --> FEDBAC	0.1	300	use
    v2	Cit --> AKG + CO2ex	ABCDEF --> ABCDE + F	0.1	300	use
    v3	AKG --> Glu	ABCDE --> ABCDE 	0.1	300	use
    v4	AKG --> Suc + CO2ex	ABCDE --> BCDE + A	0.1	300	use
    v5	Suc --> Fum	ABCD --> ABCD	0.1	300	use
    v6	Fum --> OAC	ABCD --> ABCD	0.1	300	use
    v7	OAC --> Fum	ABCD --> ABCD	0.1	300	use
    v8	Asp --> OAC	ABCD --> ABCD	0.1	300	use
    e1	Glu --> Gluex	nd 	0.1	300	use
    //Metabolites
    CO2ex	1	no	no	excreted
    AcCoA	2	no	carbonsource	no
    OAC	4	no	no	no
    Cit	6	no	no	no
    AKG	5	no	no	no
    Suc	4	symmetry	no	no
    Fum	4	symmetry	no	no
    Glu	5	no	no	no
    Gluex	5	no	no	excreted
    Asp	4	no	carbonsource	no
    //Reversible_reactions
    //Target_fragments
    Glue	gcms	Glu_12345	use
    //End

    Parameters
    ----------
    filename : metabolic model file
    format : "text" (defalut, tab deliminated) or "csv"
    output : "normal" (defalut) or "debug" (to show loaded metabolic file data)

    Returns
    -------
    reactions : Dictionary describing metabolite reactions
    reversible : Dictionary for defining reversible reactions
    metabolites : Dictionary including metabolite information
    target_fragments : Dictionary of target_fragments


    Examples
    --------
    >>> reactions, reversible, metabolites, target_fragments = load_metabolic_model("filename.txt')
    # The obtaind data (dictionaries) are directly used for generation of new Metabolic Model object
    >>> model = MetabolicModel(reactions, reversible, metabolites, target_fragments)


    """
    reactions = load_metabolic_model_reactions(filename, format=format, output = output)
    metabolites = load_metabolic_model_metabolites(filename, format=format, output = output)
    reversible_reactions = load_metabolic_model_reversibles(filename, format=format , output = output)
    target_fragments = load_metabolic_model_fragments(filename, format=format, output = output)


    return(reactions, reversible_reactions, metabolites, target_fragments)






