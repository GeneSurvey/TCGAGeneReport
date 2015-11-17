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

getNames_internal <- function(theCombinedDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jListGenesObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/GetNamesGeneEq", theCombinedDir)
	result <- .jcall(jListGenesObj, returnSig = "Z", method=theMethodString)
	if(TRUE==result)
	{
		results <- jListGenesObj$mGenes
	}
	results
}

getNamesFromMapping_internal <- function(theDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/GetMapGeneEq", theDataDir)
	result <- .jcall(jReadGeneObj, returnSig = "[S", method=theMethodString)
	if(TRUE==is.jnull(result))
	{
		result <- NULL
	}
	else
	{
		result <- unique(sort(result))
	}
	result
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

####
#### uses gene equivalent from data file
####

getNames_GeneSymbol_Mutations <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_Mutations', theVerboseFlag=theVerboseFlag)
}

getNames_GeneSymbol_RnaSeq2 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_RnaSeq2', theVerboseFlag=theVerboseFlag)
}

getNames_GeneSymbol_RnaSeq <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_RnaSeq', theVerboseFlag=theVerboseFlag)
}

getNames_GeneSymbol_SNP6 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_SNP6', theVerboseFlag=theVerboseFlag)
}

getNames_Probe_Meth450 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_Meth450', theVerboseFlag=theVerboseFlag)
}

getNames_Probe_Meth27 <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_Meth27', theVerboseFlag=theVerboseFlag)
}

getNames_CombinedHsaMimat_miRNASeq <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	getNames_internal(theCombinedDir, 'getNames_miRNASeq', theVerboseFlag=theVerboseFlag)
}

####
#### uses mapping functions
####

getNames_GeneSymbol_Meth450 <- function(theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	getNamesFromMapping_internal(theDataDir, 'getMappingGeneSymbols_Meth450', theVerboseFlag=theVerboseFlag)
}

getNames_GeneSymbol_Meth27 <- function(theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	getNamesFromMapping_internal(theDataDir, 'getMappingGeneSymbols_Meth27', theVerboseFlag=theVerboseFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
