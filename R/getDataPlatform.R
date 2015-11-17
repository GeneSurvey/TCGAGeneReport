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

getDataPlatform_internal <- function(theCombinedDir, thePlatformName, theVerboseFlag)
{
	platformData <- NULL
	load(file.path(theCombinedDir, thePlatformName, "platform.RData"), verbose=theVerboseFlag)
	platformData
}

getDataPlatform_GeneSymbol_Mutations <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'mutations', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq2 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'illuminahiseq_rnaseqv2_gene', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_RnaSeq <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'illuminahiseq_rnaseq_uncGeneRPKM', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_GeneSymbol_SNP6 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'genome_wide_snp_6_hg19nocnvWxy', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth450 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'humanmethylation450_level3', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_Probe_Meth27 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'humanmethylation27_hg19Wxy', theVerboseFlag=theVerboseFlag)
}

getDataPlatform_CombinedHsaMimat_miRNASeq <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataPlatform_internal(theCombinedDir, 'illuminahiseq_mirnaseq_isoform', theVerboseFlag=theVerboseFlag)
}
