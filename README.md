# Analysis of two renal cancer cell lines using model pools
This repository contains the code and data for the generation and analysis of pools of logical models of two renal cancer cell lines.
The model and the results are published as 

https://www.frontiersin.org/articles/10.3389/fphys.2018.01335/full

## Software
The model definition and analysis is done using the tool `TomClass`:

```R
https://github.com/hklarner/TomClass
```

and dependencies within, including `python2` and `NuSMV`, but not GUROBI (see README of TomClass).
After installing NuSMV, its path needs to be set in the variable NUSMV_CMD in TomClass:

```R
TomClass/Engine/ModelChecking.py
```

The discretization is written for `Python 3.4` or higher.

## Discretization

The folder data contains the raw data and discretized data as CSV files.
The Python script `Discretization.py` imports, discretizes and exports the data, where the threshold can be defined as mean or median.
The output files are again CSV file, with the measured species as header and the measurement time given in the first column.

## Model building, model-checking and classification

The script for defining the model, performing model-checking and analyzing the pool by classification is `Paper.py`, which runs on
`Python 2`. There are four main functions in the script that can be used:
```
create - for building the model (or pool) and creating the database, only necessary once
```
Inside the function, the components and interactions are defined as parameters. The edge labels and logical functions are formulated as clauses (see Manual of TomClass for details).
Running the function gives the number of all possible parametrizations and the number of valid models according the the clauses, which is the initial model pool.
```
annotate_compatible - is the function that performs model checking of data on pool
```
For the model-checking process, the discretized data needs to be translated into CTL formulas, for examples and details see [Thobe 2017](https://d-nb.info/1136608877/34).
```
annotate_crosstalk - annotates presence/absence of optional edges to database, prerequisite for classification
```
This function is used to annotate whether an optional edge is present or absent in a model, thus all interactions that have the clauses NotActivating or NotInhibiting should be added here in oder to add this as a feature to a model. Later we can select for these feature (which are also datasets) to classify and analyze the model pool.
```
analyse_classes - performs classification analysis
export_csv - name ouput file
```
Finally, this function runs the classification, where the features of interest need to be entered under CLASSES. The classification can also be done for a subset of models using the option RESTRICTION, e.g. selecting all datasets from one cell line gives the cell line specific pool. By adding datasets from both cell line, the classification show an empty result.

The files `full.csv`, `1257.csv` and `1851.csv` show the output of the program for either no restriction and cell line specific restrictions for MZ1257RC and MZ1851RC, respectively.

## Further Analysis

We used the output files to further analyze the model pools. In a first step, the models are sorted by the first column CROSSTALK, which is the sum of all present optional edges, thus the smallest models are on the top and the largest on the bottom. The last column SIZE gives information on how many models are in one class, which is one row of the table. The sum of the column SIZE gives the number of models in the selected pool.
Also a subset of the table can be selected, e.g. all models that have an active edge Ct_Sora_Raf. From this table, the minimal models and necessary edges can be extracted. A selection like this can also be done directly in `Paper.py` by adding Ct_Sora_Raf to the list in RESTRICTION.
In order to find the frequency of an optional edge in a pool, the column can be multiplied by SIZE, then summed up and divided by the pool size.  
