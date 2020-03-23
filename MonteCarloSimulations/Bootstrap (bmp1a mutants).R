library(globaltest)

# get list of low resolution barcodes to analyze
# list of barcode file names
barcodeHighFiles = c(
"bmp1a_11_highres.txt",
"bmp1a_12_highres.txt",
"bmp1a_14_highres.txt",
"bmp1a_16_highres.txt",
"bmp1a_17_highres.txt",
"bmp1a_19_highres.txt"
)

# get list of low resolution barcodes to analyze
# list of barcode file names
barcodeMedFiles = c(
"bmp1a_11_medres.txt",
"bmp1a_12_medres.txt",
"bmp1a_14_medres.txt",
"bmp1a_16_medres.txt",
"bmp1a_17_medres.txt",
"bmp1a_19_medres.txt"
)

# genotype of each animal
groupIDs = c("mutant","mutant","mutant","control","control","control")

paramsToAnalyze = c(13) #params to analyze

# simulation parameters
sampleSize = 3
numRuns = 1000
groups = c(rep(0,sampleSize),rep(1,sampleSize))


# define parameter names
paramNames = c(
"VertebralVolumes",
"CentrumVolumes",
"HaemalVolumes",
"NeuralVolumes",
"VertebralSAs",
"CentrumSAs",
"HaemalSAs",
"NeuralSAs",
"VertebralTMDs",
"CentrumTMDs",
"HaemalTMDs",
"NeuralTMDs",
"VertebralISs",
"CentrumISs",
"HaemalISs",
"NeuralISs",
"VertebralMTs",
"CentrumMTs",
"HaemalMTs",
"NeuralMTs",
"VertebralTSs",
"CentrumTSs",
"HaemalTSs",
"NeuralTSs",
"CentrumLength")


#Intialize medium resolution matrix
vertsToAnalyze = c(1:16) #vertebrae to include
highVertToAnalyze = 2

groupings1 = t(array(c(
1,2,2,
1,3,3,
2,1,1,
2,3,3,
1,2,3,
3,1,1,
2,2,3),dim=c(3,7)))

groupings2 = t(array(c(
1,2,
1,3,
1,4,
1,5,
1,6,
1,7,
2,3,
2,4,
2,5,
2,6,
2,7,
3,4,
3,5,
3,6,
3,7,
4,5,
4,6,
4,7,
5,6,
5,7,
6,7),dim=c(2,21)))

groupings3= groupings1+3

groupings4 = rbind(groupings2,
t(array(c(
1,1,
2,2,
3,3,
4,4,
5,5,
6,6,
7,7),dim=c(2,7))))


allMedData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze)*length(paramsToAnalyze))

#Add data for each fish into the matrix
for (i in 1:length(groupIDs)) {
	allDataInit = read.table(barcodeMedFiles[i],header=T)
	allMedData[i,] = unlist(allDataInit[vertsToAnalyze,paramsToAnalyze+1])
}

#Initialize high resolution matrix
vertsToAnalyze = c(1:3) #vertebrae to include
allHighData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze)*length(paramsToAnalyze))

#Add data for each fish into the matrix
for (i in 1:length(groupIDs)) {
	allDataInit = read.table(barcodeHighFiles[i],header=T)
	allHighData[i,] = unlist(allDataInit[vertsToAnalyze,paramsToAnalyze+1])
}

gtSens = array(0,length(paramsToAnalyze))
gtSpec = array(0,length(paramsToAnalyze))
tlSens = array(0,length(paramsToAnalyze))
tlSpec = array(0,length(paramsToAnalyze))
thSens = array(0,length(paramsToAnalyze))
thSpec = array(0,length(paramsToAnalyze))

for (p in 1:length(paramsToAnalyze)) {
	
	lowSpecP = array(0,dim(groupings2)[1])
	lowAveSpecP = array(0,dim(groupings2)[1])
	highSpecP = array(0,dim(groupings2)[1])
	
	for (n in 1:dim(groupings2)[1]) {
		
		group1 = groupings1[groupings2[n,1],]
		group2 = groupings1[groupings2[n,2],]
		
		lowSpecP[n] = p.value(gt(groups~1,~.,data=allMedData[c(group1,group2),]))
		lowAveSpecP[n] = t.test(rowMeans(allMedData[group1,],na.rm=T),rowMeans(allMedData[group2,],na.rm=T))$p.value
		highSpecP[n] = t.test(allHighData[group1,highVertToAnalyze],allHighData[group2,highVertToAnalyze])$p.value
		
	}
	
	gtSpec[p] = sum(lowSpecP<0.05)/dim(groupings2)[1]
	tlSpec[p] = sum(lowAveSpecP<0.05)/dim(groupings2)[1]
	thSpec[p] = sum(highSpecP<0.05)/dim(groupings2)[1]
	
	lowSensP = array(0,dim(groupings4)[1])
	lowAveSensP = array(0,dim(groupings4)[1])
	highSensP = array(0,dim(groupings4)[1])
	
	for (n in 1:dim(groupings4)[1]) {
		
		group1 = groupings1[groupings4[n,1],]
		group2 = groupings3[groupings4[n,2],]
		
		lowSensP[n] = p.value(gt(groups~1,~.,data=allMedData[c(group1,group2),]))
		lowAveSensP[n] = t.test(rowMeans(allMedData[group1,],na.rm=T),rowMeans(allMedData[group2,],na.rm=T))$p.value
		highSensP[n] = t.test(allHighData[group1,highVertToAnalyze],allHighData[group2,highVertToAnalyze])$p.value
		
	}
	
	gtSens[p] = sum(lowSensP<0.05)/dim(groupings4)[1]
	tlSens[p] = sum(lowAveSensP<0.05)/dim(groupings4)[1]
	thSens[p] = sum(highSensP<0.05)/dim(groupings4)[1]
	


}
