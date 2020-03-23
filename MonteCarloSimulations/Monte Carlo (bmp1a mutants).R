library(globaltest)
library(psych)
library(VineCopula)
library(copula)
library(MVN)
library(MASS)
library(matrixStats)
set.seed(100)

#####################################################################################

bmp1aHighBarcodeFiles = c(
"bmp1a_11_highres.txt",
"bmp1a_12_highres.txt",
"bmp1a_14_highres.txt",
"bmp1a_16_highres.txt",
"bmp1a_17_highres.txt",
"bmp1a_19_highres.txt"
)

# get list of low resolution barcodes to analyze

# list of barcode file names
bmp1aMedBarcodeFiles = c(
"bmp1a_11_medres.txt",
"bmp1a_12_medres.txt",
"bmp1a_14_medres.txt",
"bmp1a_16_medres.txt",
"bmp1a_17_medres.txt",
"bmp1a_19_medres.txt"
)

bmp1aGroupIDs = c("mutant","mutant","mutant","control","control","control")

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

paramsToAnalyze = c(10:13) #params to analyze
vertsToAnalyze = c(1:16)
numSimulations = 10000
sampleSizeVector = c(3)
alphaVector = c(0.05)
singleVertToAnalyze = 2


gtSensData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))
gtSpecData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))
tlSensData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))
tlSpecData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))
thSensData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))
thSpecData = array(0,dim=c(length(paramsToAnalyze),length(sampleSizeVector),length(alphaVector)))

for (k in 1:length(paramsToAnalyze)) {
	
	print(k)
	p = paramsToAnalyze[k]
			
	#read data
	bmp1aMedAllData = matrix(0,nrow=length(bmp1aGroupIDs),ncol=length(vertsToAnalyze))
	bmp1aHighAllData = matrix(0,nrow=length(bmp1aGroupIDs),ncol=length(vertsToAnalyze))
		
	# retrieve values
	for (i in 1:length(bmp1aGroupIDs)) {
		allDataInit = read.table(bmp1aMedBarcodeFiles[i],header=T)
		bmp1aMedAllData[i,] = unlist(allDataInit[vertsToAnalyze,p])
		allDataInit = read.table(bmp1aHighBarcodeFiles[i],header=T)
		bmp1aHighAllData[i,] = unlist(allDataInit[vertsToAnalyze,p])
	}
	
	
	bmp1aMedWTData = bmp1aMedAllData[bmp1aGroupIDs=="control",]
	bmp1aMedMutantData = bmp1aMedAllData[bmp1aGroupIDs=="mutant",]
	
	bmp1aMedWTCor = cor(bmp1aMedWTData)
	bmp1aMedMutantCor = cor(bmp1aMedMutantData)
	
	bmp1aMedWTStartVec = bmp1aMedWTCor[lower.tri(bmp1aMedWTCor)]
	bmp1aMedMutantStartVec = bmp1aMedMutantCor[lower.tri(bmp1aMedMutantCor)]
	
	bmp1aHighWTData = bmp1aHighAllData[bmp1aGroupIDs=="control",1:3]
	bmp1aHighMutantData = bmp1aHighAllData[bmp1aGroupIDs=="mutant",1:3]
	
	bmp1aHighWTCor = cor(bmp1aHighWTData)
	bmp1aHighMutantCor = cor(bmp1aHighMutantData)
	
	bmp1aHighWTStartVec = bmp1aHighWTCor[lower.tri(bmp1aHighWTCor)]
	bmp1aHighMutantStartVec = bmp1aHighMutantCor[lower.tri(bmp1aHighMutantCor)]
	
	
	wtMedCopula <- mvdc(copula=normalCopula(bmp1aMedWTStartVec,dim=length(vertsToAnalyze),dispstr="un"), margins=rep("norm",length(vertsToAnalyze)),paramMargins=list(
	list(mean=mean(bmp1aMedWTData[,1]), sd=sd(bmp1aMedWTData[,1])),
	list(mean=mean(bmp1aMedWTData[,2]), sd=sd(bmp1aMedWTData[,2])),
	list(mean=mean(bmp1aMedWTData[,3]), sd=sd(bmp1aMedWTData[,3])),
	list(mean=mean(bmp1aMedWTData[,4]), sd=sd(bmp1aMedWTData[,4])),
	list(mean=mean(bmp1aMedWTData[,5]), sd=sd(bmp1aMedWTData[,5])),
	list(mean=mean(bmp1aMedWTData[,6]), sd=sd(bmp1aMedWTData[,6])),
	list(mean=mean(bmp1aMedWTData[,7]), sd=sd(bmp1aMedWTData[,7])),
	list(mean=mean(bmp1aMedWTData[,8]), sd=sd(bmp1aMedWTData[,8])),
	list(mean=mean(bmp1aMedWTData[,9]), sd=sd(bmp1aMedWTData[,9])),
	list(mean=mean(bmp1aMedWTData[,10]), sd=sd(bmp1aMedWTData[,10])),
	list(mean=mean(bmp1aMedWTData[,11]), sd=sd(bmp1aMedWTData[,11])),
	list(mean=mean(bmp1aMedWTData[,12]), sd=sd(bmp1aMedWTData[,12])),
	list(mean=mean(bmp1aMedWTData[,13]), sd=sd(bmp1aMedWTData[,13])),
	list(mean=mean(bmp1aMedWTData[,14]), sd=sd(bmp1aMedWTData[,14])),
	list(mean=mean(bmp1aMedWTData[,15]), sd=sd(bmp1aMedWTData[,15])),
	list(mean=mean(bmp1aMedWTData[,16]), sd=sd(bmp1aMedWTData[,16]))
	))
	
	wtHighCopula <- mvdc(copula=normalCopula(bmp1aHighWTStartVec,dim=3,dispstr="un"), margins=rep("norm",3),paramMargins=list(
	list(mean=mean(bmp1aHighWTData[,1]), sd=sd(bmp1aHighWTData[,1])),
	list(mean=mean(bmp1aHighWTData[,2]), sd=sd(bmp1aHighWTData[,2])),
	list(mean=mean(bmp1aHighWTData[,3]), sd=sd(bmp1aHighWTData[,3]))
	))
			
	mutantMedCopula <- mvdc(copula=normalCopula(bmp1aMedMutantStartVec,dim=length(vertsToAnalyze),dispstr="un"), margins=rep("norm",length(vertsToAnalyze)),paramMargins=list(
	list(mean=mean(bmp1aMedMutantData[,1]), sd=sd(bmp1aMedMutantData[,1])),
	list(mean=mean(bmp1aMedMutantData[,2]), sd=sd(bmp1aMedMutantData[,2])),
	list(mean=mean(bmp1aMedMutantData[,3]), sd=sd(bmp1aMedMutantData[,3])),
	list(mean=mean(bmp1aMedMutantData[,4]), sd=sd(bmp1aMedMutantData[,4])),
	list(mean=mean(bmp1aMedMutantData[,5]), sd=sd(bmp1aMedMutantData[,5])),
	list(mean=mean(bmp1aMedMutantData[,6]), sd=sd(bmp1aMedMutantData[,6])),
	list(mean=mean(bmp1aMedMutantData[,7]), sd=sd(bmp1aMedMutantData[,7])),
	list(mean=mean(bmp1aMedMutantData[,8]), sd=sd(bmp1aMedMutantData[,8])),
	list(mean=mean(bmp1aMedMutantData[,9]), sd=sd(bmp1aMedMutantData[,9])),
	list(mean=mean(bmp1aMedMutantData[,10]), sd=sd(bmp1aMedMutantData[,10])),
	list(mean=mean(bmp1aMedMutantData[,11]), sd=sd(bmp1aMedMutantData[,11])),
	list(mean=mean(bmp1aMedMutantData[,12]), sd=sd(bmp1aMedMutantData[,12])),
	list(mean=mean(bmp1aMedMutantData[,13]), sd=sd(bmp1aMedMutantData[,13])),
	list(mean=mean(bmp1aMedMutantData[,14]), sd=sd(bmp1aMedMutantData[,14])),
	list(mean=mean(bmp1aMedMutantData[,15]), sd=sd(bmp1aMedMutantData[,15])),
	list(mean=mean(bmp1aMedMutantData[,16]), sd=sd(bmp1aMedMutantData[,16]))
	))	
	
	mutantHighCopula <- mvdc(copula=normalCopula(bmp1aHighMutantStartVec,dim=3,dispstr="un"), margins=rep("norm",3),paramMargins=list(
	list(mean=mean(bmp1aHighMutantData[,1]), sd=sd(bmp1aHighMutantData[,1])),
	list(mean=mean(bmp1aHighMutantData[,2]), sd=sd(bmp1aHighMutantData[,2])),
	list(mean=mean(bmp1aHighMutantData[,3]), sd=sd(bmp1aHighMutantData[,3]))
	))	

	
				
	for (n in 1:length(sampleSizeVector)) {

		numSamples = sampleSizeVector[n]			
		groups = c(rep(0,numSamples),rep(1,numSamples))
		gtPVector_WTvsMutant = array(0,length(numSimulations))
		gtPVector_WTvsWT = array(0,length(numSimulations))
		thPVector_WTvsMutant = array(0,length(numSimulations))
		thPVector_WTvsWT = array(0,length(numSimulations))
		tlPVector_WTvsMutant = array(0,length(numSimulations))
		tlPVector_WTvsWT = array(0,length(numSimulations))
			
		for (r in 1:numSimulations) {
		
			wtSimData = rMvdc(numSamples,wtMedCopula)
			mutantSimData = rMvdc(numSamples,mutantMedCopula)		
			wtSimDataNew = rMvdc(numSamples,wtMedCopula)
			
			gtPVector_WTvsMutant[r] = p.value(gt(groups~1,~.,data=rbind(mutantSimData,wtSimData)))
			gtPVector_WTvsWT[r] = p.value(gt(groups~1,~.,data=rbind(wtSimDataNew,wtSimData)))
			
			tlPVector_WTvsMutant[r] = t.test(rowMeans(mutantSimData,na.rm=T),rowMeans(wtSimData,na.rm=T))$p.value
			tlPVector_WTvsWT[r] = t.test(rowMeans(wtSimDataNew,na.rm=T),rowMeans(wtSimData,na.rm=T))$p.value
			
			wtSimData = rMvdc(numSamples,wtHighCopula)
			mutantSimData = rMvdc(numSamples,mutantHighCopula)		
			wtSimDataNew = rMvdc(numSamples,wtHighCopula)
			
			thPVector_WTvsMutant[r] = t.test(mutantSimData[,singleVertToAnalyze],wtSimData[,singleVertToAnalyze],na.rm=T)$p.value
			thPVector_WTvsWT[r] = t.test(wtSimDataNew[,singleVertToAnalyze],wtSimData[,singleVertToAnalyze],na.rm=T)$p.value
			
			
						
			for (a in 1:length(alphaVector)) {
				
				alpha = alphaVector[a]
				gtSensData[k,n,a] = sum(gtPVector_WTvsMutant<alpha)/numSimulations	
				gtSpecData[k,n,a] = sum(gtPVector_WTvsWT<alpha)/numSimulations	
				thSensData[k,n,a] = sum(thPVector_WTvsMutant<alpha)/numSimulations	
				thSpecData[k,n,a] = sum(thPVector_WTvsWT<alpha)/numSimulations	
				tlSensData[k,n,a] = sum(tlPVector_WTvsMutant<alpha)/numSimulations	
				tlSpecData[k,n,a] = sum(tlPVector_WTvsWT<alpha)/numSimulations	
					
			}
			
		}
			
	}	
	
}			
	
save(gtSensData,file="gtSensData.RData")	
save(gtSpecData,file="gtSpecData.RData")	
save(thSensData,file="thSensData.RData")	
save(thSpecData,file="thSpecData.RData")	
save(tlSensData,file="tlSensData.RData")	
save(tlSpecData,file="tlSpecData.RData")	

# load(file="gtSensData (p2_26v1_16n10000s3d4a05p123s2).RData")	
# load(file="gtSpecData (p2_26v1_16n10000s3d4a05p123s2).RData")	
# load(file="thSensData (p2_26v1_16n10000s3d4a05p123s2).RData")	
# load(file="thSpecData (p2_26v1_16n10000s3d4a05p123s2).RData")	
# load(file="tlSensData (p2_26v1_16n10000s3d4a05p123s2).RData")	
# load(file="tlSpecData (p2_26v1_16n10000s3d4a05p123s2).RData")	

# tableData = array(0,dim=c(length(paramsToAnalyze),18))
# p = 1
# tableData[,1] = gtSensData[,1,p,1,1]
# tableData[,2] = thSensData[,1,p,1,1]
# tableData[,3] = tlSensData[,1,p,1,1]
# tableData[,4] = 1-gtSpecData[,1,p,1,1]
# tableData[,5] = 1-thSpecData[,1,p,1,1]
# tableData[,6] = 1-tlSpecData[,1,p,1,1]
# p = 2
# tableData[,7] = gtSensData[,1,p,1,1]
# tableData[,8] = thSensData[,1,p,1,1]
# tableData[,9] = tlSensData[,1,p,1,1]
# tableData[,10] = 1-gtSpecData[,1,p,1,1]
# tableData[,11] = 1-thSpecData[,1,p,1,1]
# tableData[,12] = 1-tlSpecData[,1,p,1,1]
# p = 3
# tableData[,13] = gtSensData[,1,p,1,1]
# tableData[,14] = thSensData[,1,p,1,1]
# tableData[,15] = tlSensData[,1,p,1,1]
# tableData[,16] = 1-gtSpecData[,1,p,1,1]
# tableData[,17] = 1-thSpecData[,1,p,1,1]
# tableData[,18] = 1-tlSpecData[,1,p,1,1]

# write(t(tableData),file="test.txt",ncolumns=18)




