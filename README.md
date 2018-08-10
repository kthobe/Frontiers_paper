# Analysis of two renal cancer cell lines using model pools
This repository contains the code and data for the generation and analysis of pools of logical models of two renal cancer cell lines.
The model and the result are currently under submission.

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

The discretization is written for `Python 3.4`.

## Discretization

The folder data contains the raw data and discretized data as CSV files.
The Python script `Discretization.py` imports, discretizes and exports the data, where the threshold can be defined as mean or median.
The output files are again CSV file, with the measured species as header and the measurement time given in the first column.

## Model building, model-checking and classification
