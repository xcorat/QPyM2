# CUORE Python M2 Analysis

This repo is arranged as follows,

data : input and output for fitting,
  |-- testname1
  |-- groupname2
         |-- input.h5 : staged data and MC -  stored as dataframes
         |-- fit.h5 : fit output and component histograms. (Everything needed to replicate a fit?)
         |-- fit.nc : arviz inference data
         
qpym2 : python package that contains most of the code and sub packages. The software documentation is also found within here

doc : Documentation for the M2 analysis, including example notebooks
  |-- Introduction
  |-- Procedure
  |--

tests : tests and other related code, playground


 
