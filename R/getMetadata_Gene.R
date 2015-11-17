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

getMetadata_Gene_JavaToR <- function(theJavaObj)
{
	strandVal <- "unknown"
	if(FALSE==is.jnull(theJavaObj$mStrand))
	{
		strandVal <- theJavaObj$mStrand
	}
	return(new("Metadata_Gene", 
						 theGeneSymbol=theJavaObj$mGeneSymbol,
						 theGeneId=theJavaObj$mGeneId,
						 theVersionIndex=theJavaObj$mVersionIndex,
						 theChromosome=theJavaObj$mChromosome,
						 theLocationStart=as.numeric(theJavaObj$mLocationStart),
						 theLocationEnd=as.numeric(theJavaObj$mLocationEnd),
						 theStrand=strandVal))
}

getMetadata_Gene_internal_OLD_SINGLERESULT <- function(theGene, theDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	# class objects in a vector are always a list, not a vector
	listResults <- lapply(theGene, function(myGene)
	{
		results <- NULL
		jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/MetadataGene", theDataDir)
		foundFlag <- .jcall(jReadGeneObj, returnSig = "Z", method=theMethodString, 
												.jnew("java/lang/String",myGene))
		if(TRUE==foundFlag)
		{
			results <- getMetadata_Gene_JavaToR(jReadGeneObj)
		}
		results
	})
	names(listResults) <- theGene
	listResults
}

getMetadata_Gene_internal <- function(theGene, theDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	# class objects in a vector are always a list, not a vector
	listResults <- lapply(theGene, function(myGene)
	{
		results <- NULL
		jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/MetadataGene", theDataDir)
		resultArray <- .jcall(jReadGeneObj, returnSig = "[Lorg/mda/bcb/tcgagsdata/retrieve/MetadataGene;",
													method=theMethodString, 
													.jnew("java/lang/String",myGene))
		for(jReadGeneObj in resultArray)
		{
			results <- c(results, getMetadata_Gene_JavaToR(jReadGeneObj))
		}
		results
	})
	names(listResults) <- theGene
	listResults
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

####
#### specific probes
####

getMetadata_Gene_Mutations <- function(theGene, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadata_Gene_internal(theGene, theDataDir, 'getMetadataList_Mutations', theVerboseFlag=theVerboseFlag)
}


getMetadata_Gene_RnaSeq <- function(theGene, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadata_Gene_internal(theGene, theDataDir, 'getMetadataList_RNASeq', theVerboseFlag=theVerboseFlag)
}

getMetadata_Gene_RnaSeq2 <- function(theGene, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadata_Gene_internal(theGene, theDataDir, 'getMetadataList_RNASeqV2', theVerboseFlag=theVerboseFlag)
}

####
#### probes by location
####

getMetadata_GeneByNeighbor_RnaSeq <- function(theStartPosition, theStopPosition, 
																							theChromosome, theStrand,
																							theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", 
																							theVerboseFlag=FALSE)
{
	stopifnot(is.numeric(theStartPosition))
	stopifnot(is.numeric(theStopPosition))
	stopifnot(is.character(theChromosome))
	stopifnot(is.character(theStrand))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq started", theVerboseFlag=theVerboseFlag)
	jFindNeighborObj <- .jnew("org/mda/bcb/tcgagsdata/neighbors/FN_RNASeq", theDataDir)
	resultArray <- .jcall(jFindNeighborObj, 
												returnSig = "[Lorg/mda/bcb/tcgagsdata/retrieve/MetadataGene;", 
												method="findNeighbors", 
												.jlong(theStartPosition), 
												.jlong(theStopPosition),
												.jnew("java/lang/String",theChromosome),
												.jnew("java/lang/String",theStrand)
	)
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq after call", theVerboseFlag=theVerboseFlag)
	if(length(resultArray)>0)
	{
		verboseMessage("getMetadata_GeneByNeighbor_RnaSeq length(resultArray)=", length(resultArray), theVerboseFlag=theVerboseFlag)
		for(jReadGeneObj in resultArray)
		{
			results <- c(results, getMetadata_Gene_JavaToR(jReadGeneObj))
		}
	}
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq finished", theVerboseFlag=theVerboseFlag)
	results
}

getMetadata_GeneByNeighbor_RnaSeq2 <- function(theStartPosition, theStopPosition, 
																							theChromosome, theStrand,
																							theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", 
																							theVerboseFlag=FALSE)
{
	stopifnot(is.numeric(theStartPosition))
	stopifnot(is.numeric(theStopPosition))
	stopifnot(is.character(theChromosome))
	stopifnot(is.character(theStrand))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq2 started", theVerboseFlag=theVerboseFlag)
	jFindNeighborObj <- .jnew("org/mda/bcb/tcgagsdata/neighbors/FN_RNASeqV2", theDataDir)
	resultArray <- .jcall(jFindNeighborObj, 
												returnSig = "[Lorg/mda/bcb/tcgagsdata/retrieve/MetadataGene;", 
												method="findNeighbors", 
												.jlong(theStartPosition), 
												.jlong(theStopPosition),
												.jnew("java/lang/String",theChromosome),
												.jnew("java/lang/String",theStrand)
	)
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq2 after call", theVerboseFlag=theVerboseFlag)
	if(length(resultArray)>0)
	{
		verboseMessage("getMetadata_GeneByNeighbor_RnaSeq2 length(resultArray)=", length(resultArray), theVerboseFlag=theVerboseFlag)
		for(jReadGeneObj in resultArray)
		{
			results <- c(results, getMetadata_Gene_JavaToR(jReadGeneObj))
		}
	}
	verboseMessage("getMetadata_GeneByNeighbor_RnaSeq2 finished", theVerboseFlag=theVerboseFlag)
	results
}

#################################################################
#################################################################
#################################################################
#################################################################
