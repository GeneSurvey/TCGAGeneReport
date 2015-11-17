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

convertFileToList <- function(theRDataFile)
{
	dataObject <- NULL
	load(theRDataFile)
	myObj <- dataObject
	dataList <- list(myObj@mData, myObj@mSampleAnnotations, myObj@mGenomeAnnotations, myObj@mTissueAnnotations)
	names(dataList) <- c("data", "sample.annotations", "genome.annotations", "tissue.annotations")
	save(dataList, file=theRDataFile)
	dataList
}

readMatrix_internal <- function(theGeneEqList, theGeneReportDir, theRemoveDupFlag, theMethodString, theVerboseFlag, theDeltaFunction)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/GetDataMatrix", file.path(theGeneReportDir, "combined"))
	result <- .jcall(jObj, returnSig = "Z", method=theMethodString, 
									 .jarray(as.vector(as.character(theGeneEqList))))
	if(TRUE==result)
	{
		results <- matrixWithIssues(jObj$mGenesBySamplesValues, nrow=length(jObj$mGenes))
		colnames(results) <- jObj$mSamples
		rownames(results) <- jObj$mGenes
	}
	if (theRemoveDupFlag==TRUE)
	{
		rNames <- rownames(results)
		cNames <- colnames(results)[!duplicated(colnames(results))]
		results <- results[,!duplicated(colnames(results))]
		# have to do this as above line removes "matrixness" from matrix with single row
		results <- matrixWithIssues(as.vector(unlist(results)), nrow=length(rNames))
		colnames(results) <- cNames
		rownames(results) <- rNames
	}
	if(!is.null(theDeltaFunction))
	{
		results <- theDeltaFunction(results, theGeneReportDir=theGeneReportDir, theVerboseFlag=theVerboseFlag)
	}
	results
}

readSampleAnnotations_internal <- function (theSamples, theDataDir, theVerboseFlag, theUsePatientsFlag)
{
	# barcode (sample)
	### theSamples
	# patient id
	# shorten to patient only
	patientList <- substring(theSamples, 1, 12)
	# disease code
	diseaseList <- NULL
	if (FALSE==theUsePatientsFlag)
	{
		diseaseList <- getMetadataPop_BarcodeDisease_forList(theSamples, theDataDir, theVerboseFlag)
	}
	else
	{
		diseaseList <- getMetadataPop_PatientDisease_forList(theSamples, theDataDir, theVerboseFlag)
	}
	# disease string
	diseaseStringList <- as.vector(unlist(lapply(diseaseList, function(theId)
	{
		getMetadataTcga_Name_Disease(theId, theDataDir, theVerboseFlag)
	})))
	# sample type code
	# sample type string
	typeCodeList <- NULL
	typeStringList <- NULL
	if (FALSE==theUsePatientsFlag)
	{
		typeCodeList <- getMetadataPop_BarcodeSamplecode_forList(theSamples, theDataDir, theVerboseFlag)
		typeStringList <- as.vector(unlist(lapply(typeCodeList, function(theId)
		{
			getMetadataTcga_Name_SampleType(theId, theDataDir, theVerboseFlag)
		})))
	}
	else
	{
		typeCodeList <- as.vector(unlist(lapply(theSamples, function(theId)
		{
			"Delta"
		})))
		typeStringList <- as.vector(unlist(lapply(theSamples, function(theId)
		{
			"Delta(Tumor-Normal)"
		})))
	}
	# data frame
	dataFrame <- data.frame(sample=theSamples,
													patient=patientList,
													diseaseCode=diseaseList,
													diseaseString=diseaseStringList,
													typeCode=typeCodeList,
													typeString=typeStringList,
													check.names=FALSE,
													stringsAsFactors=FALSE)
}

readProbeAnnotations_internal <- function (theGeneEqList, theMetadataFunction, theDataDir, theVerboseFlag, 
																					 theGeneEqToGeneMap)
{
	# get metadata
	metadataList <- as.vector(unlist(lapply(theGeneEqList, function(theGeneEq)
	{
		theMetadataFunction(theGeneEq, theDataDir=theDataDir, theVerboseFlag=theVerboseFlag)
	})))
	names(metadataList) <- theGeneEqList
	#probe id
	probeId <- as.vector(unlist(theGeneEqList))
	#fromGene
	fromGene <- NULL
	if (is.null(theGeneEqToGeneMap))
	{
		fromGene <- as.vector(unlist(lapply(theGeneEqList, function(theGeneEq)
		{
			NA
		})))
	}
	else
	{
		fromGene <- as.vector(unlist(lapply(theGeneEqList, function(theGeneEq)
		{
			genes <- unique(sort(as.vector(unlist(lapply(names(theGeneEqToGeneMap), function(theParentGene, theGeneEqToGeneMap)
			{
				parentGene <- NULL
				if (theGeneEq %in% as.vector(unlist(theGeneEqToGeneMap[theParentGene])))
				{
					parentGene <- theParentGene
				}
				parentGene
			}, theGeneEqToGeneMap)))))
			if (length(genes)>1)
			{
				genes <- paste(genes, collapse=" | ")
			}
			genes
		})))
	}
	#chromosome
	chromosome <- as.vector(unlist(lapply(metadataList, function(theProbe) { theProbe@mChromosome })))
	#strand
	strand <- as.vector(unlist(lapply(metadataList, function(theProbe) { theProbe@mStrand })))
	#bp start - bp stop
	bpLocation <- as.vector(unlist(lapply(metadataList, function(theProbe) { theProbe@mProbeLocation })))
	#geneStructure
	geneStructure <- as.vector(unlist(lapply(metadataList, function(theProbe) 
	{
		stru <- theProbe@mGeneStructure
		result <- NULL
		for(myName in names(stru))
		{
			if (is.null(result))
			{
				result <- paste(myName, stru[myName], sep=" | ")
			}
			else
			{
				result <- paste(result, paste(myName, stru[myName], sep=" | "), sep = " : ")
			}
		}
		result
	})))
	# build dataframe
	dataFrame <- data.frame(probeId=probeId,
													fromGene=fromGene,
													chromosome=chromosome,
													strand=strand,
													bpLocation=bpLocation,
													geneStructure=geneStructure,
													check.names=FALSE,
													stringsAsFactors=FALSE)
}

readGeneAnnotations_internal <- function (theGeneEqList, theMetadataFunction, theDataDir, theVerboseFlag)
{
	# get metadata
	metadataList <- as.vector(unlist(lapply(theGeneEqList, function(theGeneEq)
	{
		theMetadataFunction(theGeneEq, theDataDir=theDataDir, theVerboseFlag=theVerboseFlag)
	})))
	#gene symbol
	geneSymbol <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mGeneSymbol })))
	geneId <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mGeneId })))
	versionIndex <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mVersionIndex })))
	chromosome <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mChromosome })))
	locationStart <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mLocationStart })))
	locationEnd <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mLocationEnd })))
	strand <- as.vector(unlist(lapply(metadataList, function(theGene) { theGene@mStrand })))
	# build dataframe
	dataFrame <- data.frame(geneSymbol=geneSymbol,
													geneId=geneId,
													versionIndex=versionIndex,
													chromosome=chromosome,
													locationStart=locationStart,
													locationEnd=locationEnd,
													strand=strand,
													check.names=FALSE,
													stringsAsFactors=FALSE)
}

readHsaMimatAnnotations_internal <- function (theGeneEqList, theMetadataFunction, theDataDir, theVerboseFlag, theGeneEqToGeneMap)
{
	# TODO: redo using theMetadataFunction when available
	#shortHsaId
	shortHsaId <- lapply(theGeneEqList, function(theId) 
	{
		name <- NULL
		split <- strsplit(theId, ".", fixed=TRUE)[[1]]
		if (2==length(split))
		{
			name <- split[1]
		}
		else if (TRUE==grepl("has", theId, fixed = TRUE))
		{
			name <- theId
		}
		else
		{
			# TODO look this up from metadata function when available
			name <- "" 
		}
		name
	})
	# mimatId
	mimatId <- lapply(theGeneEqList, function(theId) 
	{
		name <- NULL
		split <- strsplit(theId, ".", fixed=TRUE)[[1]]
		if (2==length(split))
		{
			name <- split[2]
		}
		else if (TRUE==grepl("MIMAT", theId, fixed = TRUE))
		{
			name <- theId
		}
		else
		{
			# TODO look this up from metadata function when available
			name <- "" 
		}
		name
	})
	# combinedId
	combinedId <- lapply(theGeneEqList, function(theId) 
	{
		name <- NULL
		split <- strsplit(theId, ".", fixed=TRUE)[[1]]
		if (2==length(split))
		{
			name <- theId
		}
		else
		{
			# TODO look this up from metadata function when available
			name <- "" 
		}
		name
	})
	#chromosome
	#strand
	#bp start - bp stop
	#geneStructure
	# build dataframe
	dataFrame <- data.frame(shortHsaId=shortHsaId,
													mimatId=mimatId,
													combinedId=combinedId,
													check.names=FALSE,
													stringsAsFactors=FALSE)
}

readTissueAnnotations_internal <- function(theSampleIds, theDataDir, theVerboseFlag, theUsePatientsFlag)
{
	# diseaseType
	diseaseType <- NULL
	if(FALSE==theUsePatientsFlag)
	{
		diseaseType <- getMetadataPop_BarcodeDisease_forList(theSampleIds, theDataDir=theDataDir, theVerboseFlag=theVerboseFlag)
	}
	else
	{
		diseaseType <- getMetadataPop_PatientDisease_forList(theSampleIds, theDataDir=theDataDir, theVerboseFlag=theVerboseFlag)
	}
	diseaseType <- unique(sort(as.vector(unlist(diseaseType[names(diseaseType) %in% theSampleIds]))))
	# diseaseString
	diseaseString <- as.vector(unlist(lapply(diseaseType, function(theId)
	{
		getMetadataTcga_Name_Disease(theId, theDataDir, theVerboseFlag)
	})))
	# plotColor
	plotColor <- getColorsForDiseases(diseaseType)
	dataFrame <- data.frame(diseaseType=diseaseType,
													diseaseString=diseaseString,
													plotColor=plotColor,
													check.names=FALSE,
													stringsAsFactors=FALSE)
}

buildGeneReportDataObject_internal <- function(theGeneEqList, 
																							 theGeneReportDir, theDataDir,
																							 theReadMatrixFunction, theRemoveDupFlag, 
																							 theProbeMetadataFunction, theGeneMetadataFunction, theHsaMimatMetadataFunction, 
																							 theVerboseFlag, theDeltaFunction)
{
	geneEqVector <- theGeneEqList
	eqToGeneMap <- NULL
	if (is.list(theGeneEqList))
	{
		# means this is a list of vectors each named by a gene
		geneEqVector <- unique(sort(as.vector(unlist(theGeneEqList))))
		eqToGeneMap <- theGeneEqList
	}
	#	mData="matrix",
	dataMatrix <- readMatrix_internal(geneEqVector, theGeneReportDir, theRemoveDupFlag, theReadMatrixFunction, theVerboseFlag=theVerboseFlag, theDeltaFunction=theDeltaFunction)
	#	mSampleAnnotations="data.frame",
	sampleAnnotations <- readSampleAnnotations_internal(colnames(dataMatrix), theDataDir, theVerboseFlag, !is.null(theDeltaFunction))
	#	mTissueAnnotations="data.frame"))
	tissueAnnotations <- readTissueAnnotations_internal(colnames(dataMatrix), theDataDir, theVerboseFlag, !is.null(theDeltaFunction))
	#	mGenomeAnnotations="data.frame",
	genomeAnnotations <- NULL
	if (!is.null(theProbeMetadataFunction))
	{
		genomeAnnotations <- readProbeAnnotations_internal(rownames(dataMatrix), theProbeMetadataFunction, theDataDir, theVerboseFlag, eqToGeneMap)
	}
	else if (!is.null(theGeneMetadataFunction))
	{
		genomeAnnotations <- readGeneAnnotations_internal(rownames(dataMatrix), theGeneMetadataFunction, theDataDir, theVerboseFlag)
	}
	else if (!is.null(theHsaMimatMetadataFunction))
	{
		genomeAnnotations <- readHsaMimatAnnotations_internal(rownames(dataMatrix), theHsaMimatMetadataFunction, theDataDir, theVerboseFlag, eqToGeneMap)
	}
	#setClass("ReadGeneList", representation(
	#	mData="matrix",
	#	mSampleAnnotations="data.frame",
	#	mGenomeAnnotations="data.frame",
	#	mTissueAnnotations="data.frame"))
	dataObject <- new ("GeneReportDataObject", dataMatrix, sampleAnnotations, genomeAnnotations, tissueAnnotations)
}

getDataObj_internal <- function(theGeneEqList, 
																theGeneReportDir, theDataDir,
																theReadMatrixFunction, theRemoveDupFlag, 
																theProbeMetadataFunction, theGeneMetadataFunction, theHsaMimatMetadataFunction,
																theVerboseFlag, theRDataFile, theListOnlyFlag, theDeltaFunction)
{
	dataObject <- buildGeneReportDataObject_internal(theGeneEqList, 
																									 theGeneReportDir, theDataDir,
																									 theReadMatrixFunction, theRemoveDupFlag,
																									 theProbeMetadataFunction, theGeneMetadataFunction, theHsaMimatMetadataFunction,
																									 theVerboseFlag, theDeltaFunction)
	if(FALSE==is.null(theRDataFile))
	{
		save(dataObject, file=theRDataFile)
		if (TRUE==theListOnlyFlag)
		{
			dataObject <- convertFileToList(theRDataFile)
		}
	}
	dataObject	
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

loadDataObject <- function(theRDataFile)
{
	stopifnot(isValidDirectoryPath(theRDataFile))
	dataObject <- NULL
	load(theRDataFile)
	dataObject
}

####
#### use gene eq
####

getDataObject_GeneSymbol_RnaSeq2 <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																						 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																						 theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_RnaSeq2
	}
	getDataObj_internal(theGeneEqList=theGeneEqList, 
											theGeneReportDir=theGeneReportDir,theDataDir=file.path(theGeneReportDir, "data"), 
											theReadMatrixFunction='getDataMatrix_RnaSeq2', theRemoveDupFlag=theRemoveDupFlag, 
											theProbeMetadataFunction=NULL, theGeneMetadataFunction=getMetadata_Gene_RnaSeq2, theHsaMimatMetadataFunction=NULL,
											theVerboseFlag=theVerboseFlag, theRDataFile=theRDataFile, 
											theListOnlyFlag=theListOnlyFlag, theDeltaFunction=deltaFunction)
}

getDataObject_GeneSymbol_Mutations <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																							 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																							 theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_Mutations
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_Mutations', theRemoveDupFlag, 
											NULL, getMetadata_Gene_Mutations, NULL,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}


getDataObject_GeneSymbol_RnaSeq <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																						theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																						theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_RnaSeq
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_RnaSeq', theRemoveDupFlag, 
											NULL, getMetadata_Gene_RnaSeq, NULL,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}

getDataObject_GeneSymbol_SNP6 <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																					theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																					theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_SNP6
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_SNP6', theRemoveDupFlag, 
											NULL, getMetadata_Gene_HG19, NULL,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}

getDataObject_Probe_Meth450 <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																				theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																				theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_Meth450
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_Meth450', theRemoveDupFlag, 
											getMetadata_Probe_Meth450, NULL, NULL,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}

getDataObject_Probe_Meth27 <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																			 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																			 theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_Meth27
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_Meth27', theRemoveDupFlag,
											getMetadata_Probe_Meth27, NULL, NULL,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}

getDataObject_CombinedHsaMimat_miRNASeq <- function(theGeneEqList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																										theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																										theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneEqList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	deltaFunction <- NULL
	if(TRUE==theUseDeltaFlag)
	{
		deltaFunction <- getDelta_miRNASeq
	}
	getDataObj_internal(theGeneEqList, 
											theGeneReportDir, file.path(theGeneReportDir, "data"), 
											'getDataMatrix_miRNASeq', theRemoveDupFlag, 
											NULL, NULL, TRUE,
											theVerboseFlag, theRDataFile, theListOnlyFlag, deltaFunction)
}

####
#### use gene symbol mapping
####

getDataObject_GeneSymbol_Meth450 <- function(theGeneList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																						 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																						 theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	probeList <- as.vector(unlist(lapply(theGeneList, function(theGene)
	{
		getNames_ProbeFromGeneSymbol_Meth450(theGene, file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	})))
	names(probeList) <- theGeneList
	getDataObject_Probe_Meth450(probeList, theGeneReportDir, theRemoveDupFlag, theVerboseFlag, theRDataFile, theListOnlyFlag, theUseDeltaFlag)
}

getDataObject_GeneSymbol_Meth27 <- function(theGeneList, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																						theRemoveDupFlag=TRUE, theVerboseFlag=FALSE, theRDataFile=NULL,
																						theListOnlyFlag=FALSE, theUseDeltaFlag=FALSE)
{
	stopifnot(is.character(theGeneList))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	stopifnot((is.null(theRDataFile))||(is.character(theRDataFile)))
	stopifnot((TRUE==theListOnlyFlag)||(FALSE==theListOnlyFlag))
	stopifnot((TRUE==theUseDeltaFlag)||(FALSE==theUseDeltaFlag))
	probeList <- as.vector(unlist(lapply(theGeneList, function(theGene)
	{
		getNames_ProbeFromGeneSymbol_Meth27(theGene, file.path(theGeneReportDir,"data"), theVerboseFlag=theVerboseFlag)
	})))
	names(probeList) <- theGeneList
	getDataObject_Probe_Meth27(probeList, theGeneReportDir, theRemoveDupFlag, theVerboseFlag, theRDataFile, theListOnlyFlag, theUseDeltaFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
