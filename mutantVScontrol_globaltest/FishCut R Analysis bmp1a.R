# load libraries
library(ggplot2)
library(globaltest) #global test requires installation from Bioconductor (bioconductor.org)

##########################################################################################
# in this section user input is required


# list data files output by FishCut that you wish to be analyzed
dataFiles = c(
"bmp1a_11.txt",
"bmp1a_12.txt",
"bmp1a_14.txt",
"bmp1a_16.txt",
"bmp1a_17.txt",
"bmp1a_19.txt"
)

# enter "Mutant" or "Control" for each data file
groupIDs = t(c("Mutant","Mutant","Mutant","Control","Control","Control"))

# enter standard lengths for each data file
standardLengths = c(23.9,22.0,23.6,24.8,27.1,25.4)

# enter labels for figure legends, first for Mutant, then for Control
legendLabels = c("WT","bmp1a-/-")

# this is the numbering of the different parameters
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

# vector of parameters to analyze from the list above
paramsToAnalyze = c(2:26)

# vertebrae to analyze
vertsToAnalyze = c(1:16)

# p-value threshold
pValueThresh = 0.05

# option to set y axis to zero
setYAxisToZero = TRUE

# user input not required hereafter

##########################################################################################
# function definitions

#this function computes standard errors for each group within a data frame
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# this function collates multiple plots into a single plot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########################################################################################
# plot each parameter as a function of vertebra number

# initiate variables
plotVector = vector('list', length(paramsToAnalyze))
pValVector = array(0,length(paramsToAnalyze))
controlMeanVector = array(0,length(paramsToAnalyze))
controlSDVector = array(0,length(paramsToAnalyze))
mutantMeanVector = array(0,length(paramsToAnalyze))
mutantSDVector = array(0,length(paramsToAnalyze))

# specify parameter names, to be used as plot titles
paramNames = c( 
"Row",
expression(paste("Vert.Vol (",mu,"m3)")),
expression(paste("Cent.Vol (",mu,"m3)")),
expression(paste("Haem.Vol (",mu,"m3)")),
expression(paste("Neur.Vol (",mu,"m3)")),
expression(paste("Vert.SA (",mu,"m2)")),
expression(paste("Cent.SA (",mu,"m2)")),
expression(paste("Haem.SA (",mu,"m2)")),
expression(paste("Neur.SA (",mu,"m2)")),
expression(paste("Vert.TMD (mgHA/cm3)")),
expression(paste("Cent.TMD (mgHA/cm3)")),
expression(paste("Haem.TMD (mgHA/cm3)")),
expression(paste("Neur.TMD (mgHA/cm3)")),
expression(paste("Vert.TMD.sd (mgHA/cm3)")),
expression(paste("Cent.TMD.sd (mgHA/cm3)")),
expression(paste("Haem.TMD.sd (mgHA/cm3)")),
expression(paste("Neur.TMD.sd (mgHA/cm3)")),
expression(paste("Vert.Th (",mu,"m)")),
expression(paste("Cent.Th (",mu,"m)")),
expression(paste("Haem.Th (",mu,"m)")),
expression(paste("Neur.Th (",mu,"m)")),
expression(paste("Vert.Th.sd (",mu,"m)")),
expression(paste("Cent.Th.sd (",mu,"m)")),
expression(paste("Haem.Th.sd (",mu,"m)")),
expression(paste("Neur.Th.sd (",mu,"m)")),
expression(paste("Cent.Le (",mu,"m)"))
)

# for each parameter, generate summary values, compute p-values, and create plot object
for (k in 1:length(paramsToAnalyze)) {
	
	print(k)
	
	p = paramsToAnalyze[k]
		
	# read data files and enter values into a matrix (rows: different animals, columns: different vertebrae)
	allData = matrix(0,nrow=length(groupIDs),ncol=length(vertsToAnalyze))
	for (i in 1:length(groupIDs)) {
		allDataInit = read.table(dataFiles[i],header=T)
		allData[i,] = unlist(allDataInit[vertsToAnalyze,p])
	}
	
	# format values into data table
	idVector = array(0,length(groupIDs)*length(vertsToAnalyze))
	vertNumVector = idVector
	valVector = idVector
	subVector = idVector
	
	index = 1
	for (i in 1:length(groupIDs)) {
			for (j in 1:length(vertsToAnalyze)) {
				if (groupIDs[i] == "Mutant") idVector[index] = 1
				vertNumVector[index] = j
				subVector[index] = i
				valVector[index] = allData[i,j]
				index = index + 1
			}	
	}
	
	allDataTable = data.frame(Genotype=idVector,Vertebra=vertNumVector,Subject=subVector,Value=valVector)
	allDataTable$Subject = factor(allDataTable$Subject)
	allDataTable$Genotype = factor(allDataTable$Genotype,labels=legendLabels)
	allDataTable$Value = as.vector(allDataTable$Value)
	allDataTable$Vertebra = as.vector(allDataTable$Vertebra)
	
	# perform global test and store p values
	dataMatrix = matrix(allDataTable$Value,nrow=length(vertsToAnalyze))
	gtData = t(dataMatrix)
	colnames(gtData) = as.character(vertsToAnalyze)
	y = t((groupIDs=="Mutant")*1)
	fit = gt(y~1,~.,data=gtData)
	pValVector[k] = p.value(fit)
	
	
	# create plot of parameter as a function of vertebra number
	dataToPlot = summarySE(allDataTable,measurevar="Value",groupvars=c("Genotype","Vertebra"))
	plotObj = ggplot(dataToPlot,aes(x=Vertebra,y=Value,group=Genotype,color=Genotype))
	plotObj = plotObj + geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=.1)
	plotObj = plotObj + ggtitle(paramNames[p])
	plotObj = plotObj + geom_line() + geom_point()
	
	if (setYAxisToZero==TRUE) {
		plotObj = plotObj + scale_y_continuous(limits=c(0,1.3*max(dataToPlot$Value)))	
		
	}
		
	if (pValVector[k] > pValueThresh) {
		plotObj = plotObj + scale_color_hue(l=40,c=35)
	}
	
	plotVector[[k]] = plotObj
	
	# create summary values
	controlMeanVector[k] = mean(allData[groupIDs=="Control",],na.rm=T)
	controlSDVector[k] = sd(allData[groupIDs=="Control",],na.rm=T)
	mutantMeanVector[k] = mean(allData[groupIDs=="Mutant",],na.rm=T)
	mutantSDVector[k] = sd(allData[groupIDs=="Mutant",],na.rm=T)
	
}

# create multiplot
dev.new(width=14,height=6)
multiplot(plotlist=plotVector,cols=ceiling(sqrt(length(paramsToAnalyze))))

##########################################################################################
# for each data file, create a new txt file with _bc appended containing standard scores. these "barcode files" can be used to create heatmaps. 

for (m in 1:length(groupIDs)) {
	fishDataInit = read.table(dataFiles[m],header=T)
	fishData = fishDataInit[vertsToAnalyze,paramsToAnalyze]
		
	barcodeData = array(0,dim(fishData))
	for (n in 1:length(paramsToAnalyze)) {
			barcodeData[,n] = (fishData[,n]-controlMeanVector[n])/controlSDVector[n]
	}
	
	write(barcodeData,file=gsub(".txt","_bc.txt",dataFiles[m]),ncolumns=length(vertsToAnalyze))
}

