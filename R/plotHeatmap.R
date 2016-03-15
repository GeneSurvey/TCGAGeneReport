#TCGAGeneReport Copyright 2014, 2015, 2016 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#################################################################
#################################################################
#################################################################
#################################################################

plotHeatmapOutput <- function(theGene, theOutputDir, theProbeData, theBarcodeDiseases, theBarcodeSampleType,
															theDataType, theDataTypeLabel, theZipFile, theVerboseFlag, theReadProbeFunction, 
															theTag="")
{
	stopifnot(length(theProbeData)>0)	
	dir.create(theOutputDir, recursive = TRUE, showWarnings=FALSE)
	returnFilename <- c()
	filesub <- paste(theGene,"_",theDataType,sep="")
	if (""!=theTag)
	{
		filesub <- paste(filesub, "_", theTag,sep="")
	}
	filename <- file.path(theOutputDir, compressIntoFilename(paste(filesub, "_Heatmap.PNG", sep="")))
	returnFilename <- filename
	namesToUse <- rownames(theProbeData)[rownames(theProbeData) %in% names(theBarcodeDiseases)]
	namesToUse <- sort(namesToUse)
	plotData <- theProbeData[namesToUse, , drop=FALSE]
	plotBarcodeDiseases <- theBarcodeDiseases[namesToUse, drop=FALSE]
	diseaseLabels <- c()
	diseaseRows <- c()
	tmpLabelList <- c()
	tempDis <- as.vector(unlist(plotBarcodeDiseases))
	for (diseaseIndex in 1:length(tempDis))
	{
		disease <- tempDis[diseaseIndex]
		diseasePlace <- floor(median(which(tempDis==disease)))
		if (FALSE==disease %in% tmpLabelList)
		{
			diseaseRows <- c(diseaseRows, diseaseIndex)
			tmpLabelList <- c(tmpLabelList, disease)
		}
		if (diseaseIndex==diseasePlace)
		{
			diseaseLabels <- c(diseaseLabels, paste(disease, " (", length(which(tempDis==disease)), ")", sep=""))
		}
		else
		{
			diseaseLabels <- c(diseaseLabels, "")
		}
	}
	mainText <- paste(theGene, " : ", theDataType, "\n", theDataTypeLabel, " (N=", dim(plotData)[1], ")", sep="")
	if (""!=theTag)
	{
		mainText <- paste(theGene, " : ", theDataType, " : ", theTag, "\n", theDataTypeLabel, " (N=", dim(plotData)[1], ")", sep="")
	}
	myrowsep <- diseaseRows
	rownames(plotData) <- diseaseLabels
	CairoPNG(filename=filename, width = 2400, height = 2400, pointsize=36)
	on.exit(dev.off(), add = TRUE)
	locations <- unlist(lapply(colnames(plotData), function(theProbe, theZipFile, theVerboseFlag)
	{
		probeData <- theReadProbeFunction(theProbe, theZipFile=theZipFile)
		# stuck as a list :-(
		paste(theProbe, " @ ", probeData[[1]]@mProbeLocation, sep="")
	},theZipFile=theZipFile,theVerboseFlag=theVerboseFlag))
	locOnly <- unlist(lapply(colnames(plotData), function(theProbe, theZipFile, theVerboseFlag)
	{
		probeData <- theReadProbeFunction(theProbe, theZipFile=theZipFile)
		# stuck as a list :-(
		probeData[[1]]@mProbeLocation
	},theZipFile=theZipFile,theVerboseFlag=theVerboseFlag))
	colnames(plotData) <- locations
	plotData <- plotData[,order(locOnly), drop=FALSE]
	if ((length(locOnly)>1) && (length(as.vector(plotData))!=(sum(is.nan(as.vector(plotData))))))
	{
		heatmap.2(plotData, Rowv=NULL, Colv=NULL, dendrogram="none", main=mainText, 
							sepcolor="black", rowsep=myrowsep, trace="none",
							cexRow = (0.5 + 1/log10((length(diseaseLabels)/2.0))),
							cexCol = (0.5 + 1/log10((length(diseaseLabels)/2.0))),
							margins=c(10,5))
	}
	else
	{
		image(plotData, main=mainText)
	}
	returnFilename
}
