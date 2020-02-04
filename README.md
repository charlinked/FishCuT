# FishCuT
Scripts to run FishCut software 
FishCuT is a software package to analyze mCT image stacks in zebrafish, extract quantitative information, perform statistical analysis, visualize this data and reveal patterns and trends between different groups (mutants, treated/untreated).

Please see our Elife paper for more information on FishCuT:
https://elifesciences.org/articles/26014

The FishCuT program which runs in Mathlab through a GUI will take in a mCT image-stack in the DICOM format. Please make sure all files are in the same folder. When FishCuT is done it will auto-generate an excel file containing the results of the analysis, in addition to a number of images generated from the DICOM file, such as a maximal projection of the mCT image Stack.

Statistical analysis on FishCut data uses R and the Global test. The analysis will take into account 2 groups (eg. mutant and WT) and will analyze and visualize the quantative data produced by FishCuT. The Global test is used to identify differences between 2 groups: see our paper for more details and the following article on bioconducter for more information on teh global test: https://bioconductor.org/packages/release/bioc/vignettes/globaltest/inst/doc/GlobalTest.pdf

