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

readAsGenericMatrix <- function(theZipFile, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jObj <- .jnew("org/mda/bcb/tcgagsdata/CallFromR", theZipFile)
	result <- .jcall(jObj, returnSig = "Lorg/mda/bcb/tcgagsdata/retrieve/GetMatrixPlatform;", method=theMethodString)
	if(FALSE==is.jnull(result))
	{
		results <- matrixWithIssues(result$mGenesBySamplesValues, nrow=length(result$mGenes))
		colnames(results) <- result$mSamples
		rownames(results) <- result$mGenes
	}
	results
}

getDataPlatform_internal <- function(theZipFile, theMethodString, theVerboseFlag)
{
	platformData <- readAsGenericMatrix(theZipFile, theMethodString, theVerboseFlag)
	platformData
}

getDataPlatform_GeneSymbol_Mutations <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_MutationsPlatform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq2 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_RnaSeq2Platform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_RnaSeqPlatform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_SNP6 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_SNP6Platform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth450 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_Meth450Platform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth27 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_Meth27Platform', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_CombinedHsaMimat_miRNASeq <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'getDataMatrix_miRNASeqPlatform', theVerboseFlag=theVerboseFlag)
}
