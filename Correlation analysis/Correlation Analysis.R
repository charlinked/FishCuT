library(gplots)

# list of barcode file names
barcodeFiles = c(
"bmp1a_11.txt",
"bmp1a_12.txt",
"bmp1a_14.txt",
"bmp1a_16.txt",
"bmp1a_17.txt",
"bmp1a_19.txt",
"plod2_16.txt",
"plod2_17.txt",
"plod2_18.txt",
"plod2_13.txt",
"plod2_14.txt",
"plod2_15.txt",
"opallus_06.txt",
"opallus_07.txt",
"pissaro_05.txt",
"pissaro_06.txt",
"AB_06.txt",
"AB_07.txt",
"op1.txt",
"op2.txt",
"op3.txt",
"op4.txt",
"op5.txt",
"wt1.txt",
"wt2.txt",
"wt3.txt",
"wt4.txt",
"wt5.txt",
"zf_034.txt",
"zf_036.txt",
"zfm18.txt",
"zfm19.txt",
"zfm21.txt",
"zfm25.txt")


# this is the columns of the different parameters
# 1. Row
# 2. VertebralVolumes
# 3. CentrumVolumes
# 4. HaemalVolumes
# 5. NeuralVolumes
# 6. VertebralSAs
# 7. CentrumSAs
# 8. HaemalSAs
# 9. NeuralSAs
# 10. VertebralTMDs
# 11. CentrumTMDs
# 12. HaemalTMDs
# 13. NeuralTMDs
# 14. VertebralISs
# 15. CentrumISs
# 16. HaemalISs
# 17. NeuralISs
# 18. VertebralMTs
# 19. CentrumMTs
# 20. HaemalMTs
# 21. NeuralMTs
# 22. VertebralTSs
# 23. CentrumTSs
# 24. HaemalTSs
# 25. NeuralTSs
# 26. CentrumLength

paramsToAnalyze = c(2:26) 
vertsToAnalyze = c(1:24)

allData = array(0,dim=c(length(barcodeFiles),length(paramsToAnalyze)*length(vertsToAnalyze)))

for (k in 1:length(barcodeFiles)) {

	# retrieve values
	fishDataInit = read.table(barcodeFiles[k],header=T)
	allData[k,] = unlist(fishDataInit[vertsToAnalyze,paramsToAnalyze])
}

res = cor(allData, use="pairwise.complete.obs")	

heatmap.2(abs(res),trace="none",col=colorpanel(64, "magenta", "black", "green"),symm=T,scale="none")

median(abs(res))
