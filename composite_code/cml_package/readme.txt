ConsensusML package 
Formalizes conesus feature selection for genomics data

Main functions
* impCML - An importance function suitable for declaration with Boruta. ImpCML returns a consensus imporance metric that is a mean or median of ranks from several distinct algorithms, including (1) lasso; (2) svm; (3) random forest; and (4) xgboost. This includes an option to omit one or several functions.
* impBorutaCML - Implementation of impCML for use with Boruta wrapper function.

Utilities for Plotting and Data Summaries
Plotting and utilites functions for analyzing CML data and outputs