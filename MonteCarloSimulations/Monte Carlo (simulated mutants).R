library(globaltest)
library(psych)
library(VineCopula)
library(copula)
library(MVN)
library(MASS)
library(matrixStats)
set.seed(100)

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

#paramsToAnalyze = c(2:26) #params to analyze
#vertsToAnalyze = c(1:16)
#numSimulations = 10000
#sampleSizeVector = c(3)
#effectSizeVector = c(4)
#alphaVector = c(0.05)
#patternTypeVector = c(1,2,3)
#singleVertToAnalyze = 2

# paramsToAnalyze = c(2,10:13,18) #params to analyze
# vertsToAnalyze = c(1:16)
# numSimulations = 10000
# sampleSizeVector = c(3)
# effectSizeVector = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4)
# alphaVector = c(0.05)
# patternTypeVector = c(1,2,3)
# singleVertToAnalyze = 2

paramsToAnalyze = c(10) #params to analyze
vertsToAnalyze = c(1:16)
numSimulations = 10000
sampleSizeVector = c(2,3,4,5,6,7,8,9,10)
effectSizeVector = c(4)
alphaVector = c(0.05,0.01,0.001)
patternTypeVector = c(1)
singleVertToAnalyze = 2

# paramsToAnalyze = c(10) #params to analyze
# vertsToAnalyze = c(1:16)
# numSimulations = 10000
# sampleSizeVector = c(3)
# effectSizeVector = c(4)
# alphaVector = c(0.05,0.01,0.001)
# patternTypeVector = c(1)
# singleVertToAnalyze = 2


gtSensData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))
gtSpecData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))
tlSensData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))
tlSpecData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))
thSensData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))
thSpecData = array(0,dim=c(length(paramsToAnalyze),length(effectSizeVector),length(patternTypeVector),length(sampleSizeVector),length(alphaVector)))

for (k in 1:length(paramsToAnalyze)) {
	
	#read data
	allData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze))
	bmp1aMedAllData = matrix(0,nrow=length(bmp1aGroupIDs),ncol=length(vertsToAnalyze))
	bmp1aHighAllData = matrix(0,nrow=length(bmp1aGroupIDs),ncol=length(vertsToAnalyze))
		
	print(k)
	p = paramsToAnalyze[k]
			
	# retrieve values
	for (i in 1:length(groupIDs)) {
		allDataInit = read.table(barcodeFiles[i],header=T)
		allData[i,] = unlist(allDataInit[vertsToAnalyze,p])
	}
	
	for (i in 1:length(bmp1aGroupIDs)) {
		allDataInit = read.table(bmp1aMedBarcodeFiles[i],header=T)
		bmp1aMedAllData[i,] = unlist(allDataInit[vertsToAnalyze,p])
		allDataInit = read.table(bmp1aHighBarcodeFiles[i],header=T)
		bmp1aHighAllData[i,] = unlist(allDataInit[vertsToAnalyze,p])
	}
	
	
	bmp1aMutantMeans = colMeans(bmp1aMedAllData[bmp1aGroupIDs=="mutant",])
	bmp1aWTMeans = colMeans(bmp1aMedAllData[bmp1aGroupIDs=="control",])
	bmp1aWTSDs = colSds(bmp1aMedAllData[bmp1aGroupIDs=="control",])
	bmp1aEffectSize = abs(bmp1aMutantMeans-bmp1aWTMeans)/bmp1aWTSDs
			
	allDataCor = cor(allData)
	startVec = allDataCor[lower.tri(allDataCor)]
	wtCopula <- mvdc(copula=normalCopula(startVec,dim=length(vertsToAnalyze),dispstr="un"), margins=rep("norm",length(vertsToAnalyze)),paramMargins=list(
	list(mean=mean(allData[,1]), sd=sd(allData[,1])),
	list(mean=mean(allData[,2]), sd=sd(allData[,2])),
	list(mean=mean(allData[,3]), sd=sd(allData[,3])),
	list(mean=mean(allData[,4]), sd=sd(allData[,4])),
	list(mean=mean(allData[,5]), sd=sd(allData[,5])),
	list(mean=mean(allData[,6]), sd=sd(allData[,6])),
	list(mean=mean(allData[,7]), sd=sd(allData[,7])),
	list(mean=mean(allData[,8]), sd=sd(allData[,8])),
	list(mean=mean(allData[,9]), sd=sd(allData[,9])),
	list(mean=mean(allData[,10]), sd=sd(allData[,10])),
	list(mean=mean(allData[,11]), sd=sd(allData[,11])),
	list(mean=mean(allData[,12]), sd=sd(allData[,12])),
	list(mean=mean(allData[,13]), sd=sd(allData[,13])),
	list(mean=mean(allData[,14]), sd=sd(allData[,14])),
	list(mean=mean(allData[,15]), sd=sd(allData[,15])),
	list(mean=mean(allData[,16]), sd=sd(allData[,16]))
	))
	
	for (d in 1:length(effectSizeVector)) {
		
		print(d)
	
		effectSize = effectSizeVector[d]
		
		for (t in 1:length(patternTypeVector)) {
			
			if (patternTypeVector[t] == 1) {patternVector = array(1,length(vertsToAnalyze))}
			if (patternTypeVector[t] == 2) {patternVector = seq(from=0,to=1,length=16)}
			if (patternTypeVector[t] == 3) {patternVector = seq(from=1,to=0,length=16)}
					
			mutantCopula <- mvdc(copula=normalCopula(startVec,dim=length(vertsToAnalyze),dispstr="un"), margins=rep("norm",length(vertsToAnalyze)),paramMargins=list(
			list(mean=mean(allData[,1])+patternVector[1]*effectSize*sd(allData[,1]), sd=sd(allData[,1])),
			list(mean=mean(allData[,2])+patternVector[2]*effectSize*sd(allData[,2]), sd=sd(allData[,2])),
			list(mean=mean(allData[,3])+patternVector[3]*effectSize*sd(allData[,3]), sd=sd(allData[,3])),
			list(mean=mean(allData[,4])+patternVector[4]*effectSize*sd(allData[,4]), sd=sd(allData[,4])),
			list(mean=mean(allData[,5])+patternVector[5]*effectSize*sd(allData[,5]), sd=sd(allData[,5])),
			list(mean=mean(allData[,6])+patternVector[6]*effectSize*sd(allData[,6]), sd=sd(allData[,6])),
			list(mean=mean(allData[,7])+patternVector[7]*effectSize*sd(allData[,7]), sd=sd(allData[,7])),
			list(mean=mean(allData[,8])+patternVector[8]*effectSize*sd(allData[,8]), sd=sd(allData[,8])),
			list(mean=mean(allData[,9])+patternVector[9]*effectSize*sd(allData[,9]), sd=sd(allData[,9])),
			list(mean=mean(allData[,10])+patternVector[10]*effectSize*sd(allData[,10]), sd=sd(allData[,10])),
			list(mean=mean(allData[,11])+patternVector[11]*effectSize*sd(allData[,11]), sd=sd(allData[,11])),
			list(mean=mean(allData[,12])+patternVector[12]*effectSize*sd(allData[,12]), sd=sd(allData[,12])),
			list(mean=mean(allData[,13])+patternVector[13]*effectSize*sd(allData[,13]), sd=sd(allData[,13])),
			list(mean=mean(allData[,14])+patternVector[14]*effectSize*sd(allData[,14]), sd=sd(allData[,14])),
			list(mean=mean(allData[,15])+patternVector[15]*effectSize*sd(allData[,15]), sd=sd(allData[,15])),
			list(mean=mean(allData[,16])+patternVector[16]*effectSize*sd(allData[,16]), sd=sd(allData[,16]))
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
				
					wtSimData = rMvdc(numSamples,wtCopula)
					mutantSimData = rMvdc(numSamples,mutantCopula)		
					wtSimDataNew = rMvdc(numSamples,wtCopula)
					
					gtPVector_WTvsMutant[r] = p.value(gt(groups~1,~.,data=rbind(mutantSimData,wtSimData)))
					gtPVector_WTvsWT[r] = p.value(gt(groups~1,~.,data=rbind(wtSimDataNew,wtSimData)))
					
					thPVector_WTvsMutant[r] = t.test(mutantSimData[,singleVertToAnalyze],wtSimData[,singleVertToAnalyze],na.rm=T)$p.value
					thPVector_WTvsWT[r] = t.test(wtSimDataNew[,singleVertToAnalyze],wtSimData[,singleVertToAnalyze],na.rm=T)$p.value
					
					tlPVector_WTvsMutant[r] = t.test(rowMeans(mutantSimData,na.rm=T),rowMeans(wtSimData,na.rm=T))$p.value
					tlPVector_WTvsWT[r] = t.test(rowMeans(wtSimDataNew,na.rm=T),rowMeans(wtSimData,na.rm=T))$p.value
					
								
					for (a in 1:length(alphaVector)) {
						
						alpha = alphaVector[a]
						gtSensData[k,d,t,n,a] = sum(gtPVector_WTvsMutant<alpha)/numSimulations	
						gtSpecData[k,d,t,n,a] = sum(gtPVector_WTvsWT<alpha)/numSimulations	
						thSensData[k,d,t,n,a] = sum(thPVector_WTvsMutant<alpha)/numSimulations	
						thSpecData[k,d,t,n,a] = sum(thPVector_WTvsWT<alpha)/numSimulations	
						tlSensData[k,d,t,n,a] = sum(tlPVector_WTvsMutant<alpha)/numSimulations	
						tlSpecData[k,d,t,n,a] = sum(tlPVector_WTvsWT<alpha)/numSimulations	
										
								
					}
					
				}
					
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




