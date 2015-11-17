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

getImputedNAs_GeneEq_internal <- function(theGeneEqList, theGeneReportDir, theRemoveDupFlag, theMethodString, 
																		theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/GetImputedNAsMatrix", file.path(theGeneReportDir, "combined"))
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
	results
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

####
#### uses gene equivalent from data file
####

getImputedNAs_GeneSymbol_RnaSeq2 <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																			 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_RnaSeq2 -> getImputedNAsMatrix_RnaSeq2
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_RnaSeq2', 
													theVerboseFlag=theVerboseFlag)
}

getImputedNAs_GeneSymbol_RnaSeq <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																			theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_RnaSeq -> getImputedNAsMatrix_RnaSeq
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_RnaSeq', 
													theVerboseFlag=theVerboseFlag)
}

getImputedNAs_GeneSymbol_SNP6 <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																		theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_SNP6 -> getImputedNAsMatrix_SNP6
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_SNP6', 
													theVerboseFlag=theVerboseFlag)
}

getImputedNAs_Probe_Meth450 <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																	theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_Meth450 -> getImputedNAsMatrix_Meth450
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_Meth450', 
													theVerboseFlag=theVerboseFlag)
}

getImputedNAs_Probe_Meth27 <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																 theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_Meth27 -> getImputedNAsMatrix_Meth27
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_Meth27', 
													theVerboseFlag=theVerboseFlag)
}

getImputedNAs_CombinedHsaMimat_miRNASeq <- function(theGeneEq, theGeneReportDir="/rsrch1/bcb/batcheffects/GENE_REPORT", 
																							theRemoveDupFlag=TRUE, theVerboseFlag=FALSE)
{
	stopifnot(is.character(theGeneEq))
	stopifnot(isValidDirectoryPath(theGeneReportDir))
	stopifnot((TRUE==theRemoveDupFlag)||(FALSE==theRemoveDupFlag))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	# getImputedNAs_miRNASeq -> getImputedNAsMatrix_miRNASeq
	getImputedNAs_GeneEq_internal(theGeneEq, theGeneReportDir, theRemoveDupFlag, 'getImputedNAsMatrix_miRNASeq', 
													theVerboseFlag=theVerboseFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
