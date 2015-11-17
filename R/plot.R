#TCGAGeneReport Copyright 2014, 2015 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#################################################################
#################################################################
# internal
#################################################################
#################################################################

plot_fromGene_internal <- function(theGene, theOutputDir, theGeneReportDir, 
																	 theMapGeneToProbesFunction, theReadProbeFunction, theReadGeneFunction,
																	 theDataType, theDataTypeLabel, theVerboseFlag, theUseDeltaFlag=FALSE,
																	 theReplicateFlag=FALSE)
{
	dir.create(theOutputDir, recursive = TRUE, showWarnings=FALSE)
	returnFilename <- c()
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myPathCombined <- NULL
	myGeneEqPre <- NULL
	# read gene data
	if (FALSE==theUseDeltaFlag)
	{
		verboseMessage("plot_fromGene_internal FALSE==theUseDeltaFlag", theVerboseFlag=theVerboseFlag)
		myTag <- ""
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	else
	{
		verboseMessage("plot_fromGene_internal TRUE==theUseDeltaFlag", theVerboseFlag=theVerboseFlag)
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	probeList <- as.vector(unlist(theMapGeneToProbesFunction(theGene, file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)))
	locations <- as.vector(unlist(lapply(probeList, function(theProbe, theDataDir, theVerboseFlag=theVerboseFlag)
	{
		probeData <- as.vector(unlist(theReadProbeFunction(theProbe, theDataDir=theDataDir, theVerboseFlag=theVerboseFlag)))
		# this is only one value, but you cannot unlist a list/vector of classes, I think
		probeData[[1]]@mProbeLocation
	},theDataDir=file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)))
	probeList <- probeList[order(locations)]
	maxLocationSize <- max(as.vector(unlist(lapply(locations, nchar))))
	probeDataList <- lapply(probeList, function(theProbe, thePath, theVerboseFlag)
	{
		theReadGeneFunction(theProbe, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	}, thePath=myPathCombined, theVerboseFlag=theVerboseFlag)
	# first write individual probe plots
	for (index in 1:length(probeDataList))
	{
		myGeneData <- probeDataList[[index]]
		probe <- probeList[index]
		myProbeData <- as.vector(unlist(theReadProbeFunction(probe, theDataDir=file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)))
		# this is only one value, but you cannot unlist a list/vector of classes, I think
		myLocation <- myProbeData[[1]]@mProbeLocation
		if(!is.na(as.numeric(myLocation)))
		{
			myLocation <- sprintf(paste("%0", maxLocationSize, "d", sep=""), as.numeric(myLocation))
		}
		# this is only one value, but you cannot unlist a list/vector of classes, I think
		geneEqPre <- paste(padChromosomeName(myProbeData[[1]]@mChromosome), "-", myLocation, "-", sep="")
		returnFilename <- c(returnFilename, 
												plotGenericOutput(paste(theGene, probe, sep="-"), 
																					theOutputDir, myGeneData, myBarcodeDiseases, 
																					myBarcodeSampleType, theDataType, theDataTypeLabel, 
																					theVerboseFlag=theVerboseFlag, 
																					theTag=myTag, theGeneEqPre=geneEqPre,
																					theReplicateFlag=theReplicateFlag))
	}
	# now do combined probe plots
	# first write individual probe plots
	#for (index in 1:length(probeDataList))
	#{
	#	myGeneData <- probeDataList[[index]]
	#	probe <- probeList[index]
	#	message(paste("probe=",probe," ",paste(myGeneData[1:3],collapse=", ", sep=""), sep=""))
	#}
	numberOfRows <- length(probeDataList[[1]])
	sampleList <- colnames(probeDataList[[1]])
	probeDataList <- as.vector(unlist(probeDataList))
	probeData <- matrixWithIssues(probeDataList, nrow=numberOfRows)
	colnames(probeData) <- probeList
	rownames(probeData) <- sampleList
	# probeData[sample,probe]
	returnFilename <- c(returnFilename,
											plotHeatmapOutput(theGene, theOutputDir, probeData, myBarcodeDiseases, myBarcodeSampleType, 
																				theDataType, theDataTypeLabel, file.path(theGeneReportDir,"data"), 
																				theVerboseFlag, theReadProbeFunction, myTag))
	returnFilename
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

####
#### uses gene equivalent from data file
####

plot_GeneSymbol_Mutations <- function(theGeneEq, theOutputDir, 
																		theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																		theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_GeneSymbol_Mutations(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "Mutations", "Count", theVerboseFlag=theVerboseFlag, 
										theTag=myTag, theGeneEqPre="", theReplicateFlag=theReplicateFlag)
}

plot_GeneSymbol_RnaSeq2 <- function(theGeneEq, theOutputDir, 
																		theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
												 theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_GeneSymbol_RnaSeq2(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "RNASeqV2", "Log Normalized Count", theVerboseFlag=theVerboseFlag, 
										theTag=myTag, theGeneEqPre="", theReplicateFlag=theReplicateFlag)
}

plot_GeneSymbol_RnaSeq <- function(theGeneEq, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
												 theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_GeneSymbol_RnaSeq(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "RNASeq", "Log Normalized Count", theVerboseFlag=theVerboseFlag, 
										theTag=myTag, theGeneEqPre="", theReplicateFlag=theReplicateFlag)
}

plot_GeneSymbol_SNP6 <- function(theGeneEq, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
											theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_GeneSymbol_SNP6(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "SNP6", "Log Normalized Copy Number", theVerboseFlag=theVerboseFlag, 
										theTag=myTag, theGeneEqPre="", theReplicateFlag=theReplicateFlag)
}

plot_Probe_Meth450 <- function(theGeneEq, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
												 theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myProbe <- getMetadata_Probe_Meth450(theGeneEq, file.path(theGeneReportDir,"data"))
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_Probe_Meth450(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	geneEqPre <- paste(padChromosomeName(myProbe[[1]]@mChromosome), "-", myProbe[[1]]@mProbeLocation, "-", sep="")
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "Meth450", "Beta Value", theVerboseFlag=theVerboseFlag, 
										theTag=myTag, theGeneEqPre=geneEqPre, theReplicateFlag=theReplicateFlag)
}

plot_Probe_Meth27 <- function(theGeneEq, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
												theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myProbe <- getMetadata_Probe_Meth27(theGeneEq, file.path(theGeneReportDir,"data"))
	stopifnot(FALSE==is.null(myProbe[[1]]))
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_Probe_Meth27(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	geneEqPre <- paste(padChromosomeName(myProbe[[1]]@mChromosome), "-", myProbe[[1]]@mProbeLocation, "-", sep="")
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "Meth27", "Beta Value", theVerboseFlag=theVerboseFlag,
										theTag=myTag, theGeneEqPre=geneEqPre, theReplicateFlag=theReplicateFlag)
}

plot_CombinedHsaMimat_miRNASeq <- function(theGeneEq, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
													theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(1==length(theGeneEq))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	# read gene data
	myGeneData <- NULL
	myBarcodeDiseases <- NULL
	myBarcodeSampleType <- NULL
	myTag <- ""
	myGeneData <- getData_CombinedHsaMimat_miRNASeq(theGeneEq, theGeneReportDir, theVerboseFlag=theVerboseFlag, theUseDeltaFlag=theUseDeltaFlag)
	if(TRUE==theUseDeltaFlag)
	{
		myTag <- "Delta"
		myBarcodeDiseases <- getMetadataPop_PatientDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- myBarcodeDiseases
		myBarcodeSampleType[myBarcodeSampleType] <- "Delta"
	}
	else
	{
		myBarcodeDiseases <- getMetadataPop_BarcodeDisease(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
		myBarcodeSampleType <- getMetadataPop_BarcodeSamplecode(file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	}
	plotGenericOutput(theGeneEq, theOutputDir, myGeneData, myBarcodeDiseases, 
										myBarcodeSampleType, "miRNASeq", "isoform", theVerboseFlag=theVerboseFlag,
										theTag=myTag, theGeneEqPre="", theReplicateFlag=theReplicateFlag)
}

####
#### uses mapping functions
####


plot_GeneSymbol_Meth450 <- function(theGene, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																		theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(1==length(theGene))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	plot_fromGene_internal(theGene, theOutputDir, theGeneReportDir, 
												 theMapGeneToProbesFunction=getNames_ProbeFromGeneSymbol_Meth450, 
												 theReadProbeFunction=getMetadata_Probe_Meth450, 
												 theReadGeneFunction=getData_Probe_Meth450,
												 theDataType="Meth450", 
												 theDataTypeLabel="Beta Value", 
												 theVerboseFlag=theVerboseFlag,
												 theUseDeltaFlag=theUseDeltaFlag,
												 theReplicateFlag=theReplicateFlag)
}

plot_GeneSymbol_Meth27 <- function(theGene, theOutputDir, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																 theVerboseFlag=FALSE, theUseDeltaFlag=FALSE, theReplicateFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(1==length(theGene))
	stopifnot(isValidDirectoryPath(theOutputDir))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	stopifnot((TRUE==theReplicateFlag)||(FALSE==theReplicateFlag))
	plot_fromGene_internal(theGene, theOutputDir, theGeneReportDir, 
												 theMapGeneToProbesFunction=getNames_ProbeFromGeneSymbol_Meth27, 
												 theReadProbeFunction=getMetadata_Probe_Meth27, 
												 theReadGeneFunction=getData_Probe_Meth27,
												 theDataType="Meth27", 
												 theDataTypeLabel="Beta Value", 
												 theVerboseFlag=theVerboseFlag,
												 theUseDeltaFlag=theUseDeltaFlag,
												 theReplicateFlag=theReplicateFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
