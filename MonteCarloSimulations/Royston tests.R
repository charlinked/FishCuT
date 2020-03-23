library(MVN)

#####################################################################################

barcodeFiles = c("AB-01-033017 fish01.DCM.txt",
"AB-01-033017 fish02.DCM.txt",
"AB-01-033017 fish03.DCM.txt",
"AB-01-033017 fish04.DCM.txt",
"AB-01-033017 fish05.DCM.txt",
"AB-01-033017 fish06.DCM.txt",
"AB-01-033017 fish07.DCM.txt",
"AB-01-033017 fish08.DCM.txt",
"AB-01-033017 fish09.DCM.txt",
"AB-01-033017 fish10.DCM.txt",
"AB-01-033017 fish11.DCM.txt",
"AB-01-033017 fish12.DCM.txt",
"AB-01-033017 fish13.DCM.txt",
"AB-01-033017 fish14.DCM.txt",
"AB-01-033017 fish15.DCM.txt",
"AB-01-033017 fish16.DCM.txt")

# genotype of each animal
groupIDs = t(c("AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","AB"))
standardLengths = c(21.1,21.0,20.7,19.7,20.9,19.3,20.0,21.4,20.4,20.8,20.7,18.0,19.2,21.1,21.4,20.6)


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

paramsToAnalyze = c(2:26) #params to analyze
#vertsToAnalyze = c(1,3,5,7,9,11,13,15)
vertsToAnalyze = c(2,4,6,8,10,12,14,16)

roystonResults = array(0,dim=length(paramsToAnalyze))

for (k in 1:length(paramsToAnalyze)) {
	
	#read data
	allData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze))
		
	print(k)
	p = paramsToAnalyze[k]
			
	# retrieve values
	for (i in 1:length(groupIDs)) {
		allDataInit = read.table(barcodeFiles[i],header=T)
		allData[i,] = unlist(allDataInit[vertsToAnalyze,p])
	}
	
	roystonResults[k] = attributes(roystonTest(allData))$p.value
}

 
