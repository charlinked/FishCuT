# FishCuT

FishCuT is a software package to analyze mCT image stacks in zebrafish, extract quantitative information, perform statistical analysis and visualize patterns and trends between different groups (mutants, treated/untreated).

Please see our Elife paper for more information on FishCuT:
https://elifesciences.org/articles/26014

The first part of the FishCuT program runs in Matlab through a GUI, takes in a mCT image-stack in the DICOM format and produces a .txt file with quantitative information. This data can be further analyzed in R - We are using multiple regression analysis (the Global test) to identify differences between 2 groups. Posted in this folder are the R scripts needed for this (incl. sample data), as well as the scripts used for validation of the FishCut methodology and software (details below) - rationale, results and methodology discussed in detail in our paper.


