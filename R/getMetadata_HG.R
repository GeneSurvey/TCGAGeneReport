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

getMetadata_HG_JavaToR <- function(theJavaObj)
{
	strandVal <- "unknown"
	if(FALSE==is.jnull(theJavaObj$mStrand))
	{
		strandVal <- theJavaObj$mStrand
	}
	message("theJavaObj$mGeneSymbol = ",theJavaObj$mGeneSymbol )
	message("theJavaObj$mGeneId = ",theJavaObj$mGeneId )
	message("theJavaObj$mVersionIndex = ", theJavaObj$mVersionIndex)
	return(new("Metadata_Gene", 
						 theGeneSymbol=theJavaObj$mGeneSymbol,
						 theGeneId=theJavaObj$mGeneId,
						 theVersionIndex=theJavaObj$mVersionIndex,
						 theChromosome=theJavaObj$mChromosome,
						 theLocationStart=as.numeric(theJavaObj$mLocationStart),
						 theLocationEnd=as.numeric(theJavaObj$mLocationEnd),
						 theStrand=strandVal))
}

getMetadata_HG_internal <- function(theGene, theDataDir, theMethodString, theVerboseFlag)
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
		for(myObj in resultArray)
		{
			results <- c(results, getMetadata_HG_JavaToR(myObj))
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

getMetadata_Gene_HG19 <- function(theGene, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadata_HG_internal(theGene, theDataDir, 'getMetadataList_HG19', theVerboseFlag=theVerboseFlag)
}

getMetadata_Gene_HG18 <- function(theGene, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGene))
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadata_HG_internal(theGene, theDataDir, 'getMetadataList_HG18', theVerboseFlag=theVerboseFlag)
}

####
#### probes by location
####

getMetadata_GeneByNeighbor_HG18 <- function(theStartPosition, theStopPosition, 
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
	verboseMessage("getMetadata_GeneByNeighbor_HG18 started", theVerboseFlag=theVerboseFlag)
	jFindNeighborObj <- .jnew("org/mda/bcb/tcgagsdata/neighbors/GS_HG18", theDataDir)
	resultArray <- .jcall(jFindNeighborObj, 
												returnSig = "[Lorg/mda/bcb/tcgagsdata/retrieve/MetadataGene;", 
												method="findNeighbors", 
												.jlong(theStartPosition), 
												.jlong(theStopPosition),
												.jnew("java/lang/String",theChromosome),
												.jnew("java/lang/String",theStrand)
	)
	verboseMessage("getMetadata_GeneByNeighbor_HG18 after call", theVerboseFlag=theVerboseFlag)
	if(length(resultArray)>0)
	{
		verboseMessage("getMetadata_GeneByNeighbor_HG18 length(resultArray)=", length(resultArray), theVerboseFlag=theVerboseFlag)
		for(jReadGeneObj in resultArray)
		{
			results <- c(results, getMetadata_HG_JavaToR(jReadGeneObj))
		}
	}
	verboseMessage("getMetadata_GeneByNeighbor_HG18 finished", theVerboseFlag=theVerboseFlag)
	results
}

getMetadata_GeneByNeighbor_HG19 <- function(theStartPosition, theStopPosition, 
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
	verboseMessage("getMetadata_GeneByNeighbor_HG19 started", theVerboseFlag=theVerboseFlag)
	jFindNeighborObj <- .jnew("org/mda/bcb/tcgagsdata/neighbors/GS_HG19", theDataDir)
	resultArray <- .jcall(jFindNeighborObj, 
												returnSig = "[Lorg/mda/bcb/tcgagsdata/retrieve/MetadataGene;", 
												method="findNeighbors", 
												.jlong(theStartPosition), 
												.jlong(theStopPosition),
												.jnew("java/lang/String",theChromosome),
												.jnew("java/lang/String",theStrand)
	)
	verboseMessage("getMetadata_GeneByNeighbor_HG19 after call", theVerboseFlag=theVerboseFlag)
	if(length(resultArray)>0)
	{
		verboseMessage("getMetadata_GeneByNeighbor_HG19 length(resultArray)=", length(resultArray), theVerboseFlag=theVerboseFlag)
		for(jReadGeneObj in resultArray)
		{
			results <- c(results, getMetadata_HG_JavaToR(jReadGeneObj))
		}
	}
	verboseMessage("getMetadata_GeneByNeighbor_HG19 finished", theVerboseFlag=theVerboseFlag)
	results
}

#################################################################
#################################################################
#################################################################
#################################################################
