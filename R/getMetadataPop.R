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

getMetadataPop_internal <- function(theDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/MetadataPop", theDataDir)
	result <- .jcall(jReadGeneObj, returnSig = "Z", method=theMethodString)
	if(TRUE==result)
	{
		results <- jReadGeneObj$mValues
		names(results) <- jReadGeneObj$mIds
	}
	results
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

getMetadataPop_BarcodeDisease <- function(theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadataPop_internal(theDataDir, 'getMetadataPop_BarcodeDisease', theVerboseFlag=theVerboseFlag)
}

getMetadataPop_BarcodeSamplecode <- function(theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadataPop_internal(theDataDir, 'getMetadataPop_BarcodeSamplecode', theVerboseFlag=theVerboseFlag)
}

getMetadataPop_PatientDisease <- function(theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theDataDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getMetadataPop_internal(theDataDir, 'getMetadataPop_PatientDisease', theVerboseFlag=theVerboseFlag)
}

### used outside this file, but not exported

getMetadataPop_BarcodeDisease_forList <- function(theList, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	results <- getMetadataPop_BarcodeDisease(theDataDir, theVerboseFlag=theVerboseFlag)
	results <- results[theList]
	names(results) <- theList
	results <- as.vector(unlist(lapply(results, function(theVal)
	{
		if (is.na(theVal))
		{
			theVal <- "UNK"
		}
		theVal
	})))
	names(results) <- theList
	results
}

getMetadataPop_BarcodeSamplecode_forList <- function(theList, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	results <- getMetadataPop_BarcodeSamplecode(theDataDir, theVerboseFlag=theVerboseFlag)
	results <- results[theList]
	names(results) <- theList
	results <- as.vector(unlist(lapply(results, function(theVal)
	{
		if (is.na(theVal))
		{
			theVal <- "UNK"
		}
		theVal
	})))
	names(results) <- theList
	results
}

getMetadataPop_PatientDisease_forList <- function(theList, theDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/data", theVerboseFlag=FALSE)
{
	results <- getMetadataPop_PatientDisease(theDataDir, theVerboseFlag=theVerboseFlag)
	results <- results[theList]
	names(results) <- theList
	results <- as.vector(unlist(lapply(results, function(theVal)
	{
		if (is.na(theVal))
		{
			theVal <- "UNK"
		}
		theVal
	})))
	names(results) <- theList
	results
}

#################################################################
#################################################################
#################################################################
#################################################################
