# get list of high resolution barcodes to analyze
# list of barcode file names
barcodeHighFiles = c(
"bmp1a_11_medres.txt",
"bmp1a_12_medres.txt",
"bmp1a_14_medres.txt",
"bmp1a_16_medres.txt",
"bmp1a_17_medres.txt",
"bmp1a_19_medres.txt"
)

# get list of low resolution barcodes to analyze
# list of barcode file names
barcodeMedFiles = c(
"bmp1a_11_highres.txt",
"bmp1a_12_highres.txt",
"bmp1a_14_highres.txt",
"bmp1a_16_highres.txt",
"bmp1a_17_highres.txt",
"bmp1a_19_highres.txt"
)

groupIDs = c("Mutant","Mutant","Mutant","Control","Control","Control")
paramsToAnalyze = c(2:4,10:12,18:20,25) #params to analyze
vertsToAnalyze = c(2) #vertebrae to include

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

#Initialize high resolution matrix
allHighData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze)*length(paramsToAnalyze))

#Add data for each fish into the matrix
for (i in 1:length(groupIDs)) {
	allDataInit = read.table(barcodeHighFiles[i],header=T)
	allHighData[i,] = unlist(allDataInit[vertsToAnalyze,paramsToAnalyze+1])
}

#Intialize medium resolution matrix
allMedData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze)*length(paramsToAnalyze))

#Add data for each fish into the matrix
for (i in 1:length(groupIDs)) {
	allDataInit = read.table(barcodeMedFiles[i],header=T)
	allMedData[i,] = unlist(allDataInit[vertsToAnalyze,paramsToAnalyze+1])
}

#Initialize data vectors
rsquaredvec = array(0,dim(allHighData)[2])
slopevec = array(0,dim(allHighData)[2])
sevec = array(0,dim(allHighData)[2])
pValVecMed = array(0,dim(allHighData)[2])
pValVecHigh = array(0,dim(allHighData)[2])


for (k in 1:dim(allHighData)[2]) {
	y = allHighData[,k]
	x = allMedData[,k]
	fit = lm(y~x-1)
	rsquaredvec[k] = summary(fit)$r.squared
	slopevec[k] = summary(fit)$coef[1]
	sevec[k] = summary(fit)$coef[2]
	pValVecMed[k] = t.test(allMedData[1:3,k],allMedData[4:6,k])$p.val
	pValVecHigh[k] = t.test(allHighData[1:3,k],allHighData[4:6,k])$p.val
}
	