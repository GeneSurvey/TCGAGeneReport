#TCGAGeneReport Copyright 2014, 2015 University of Texas MD Anderson Cancer Center
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

getDataPlatform_internal <- function(theZipFile, thePlatformName, theVerboseFlag)
{
	platformData <- NULL
	zippedDataStream <- unz(theZipFile, paste("combined", thePlatformName, "platform.RData", sep="/"))
	load(filezippedDataStream, verbose=theVerboseFlag)
	platformData
}

getDataPlatform_GeneSymbol_Mutations <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'mutations', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq2 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'illuminahiseq_rnaseqv2_gene', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'illuminahiseq_rnaseq_uncGeneRPKM', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_SNP6 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'genome_wide_snp_6_hg19nocnvWxy', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth450 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'humanmethylation450_level3', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth27 <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'humanmethylation27_hg19Wxy', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_CombinedHsaMimat_miRNASeq <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	stopifnot(file.exists(theZipFile))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theZipFile, 'illuminahiseq_mirnaseq_isoform', theVerboseFlag=theVerboseFlag)
}
