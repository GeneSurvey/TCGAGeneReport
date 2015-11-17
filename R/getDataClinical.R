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

getDataClinical_internal <- function(theCombinedDir, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/GetDataClinical", theCombinedDir)
	result <- .jcall(jObj, returnSig = "Z", method="getDataClinical")
	if(TRUE==result)
	{
		results <- matrixWithIssues(jObj$mGenesBySamplesValues, ncol=length(jObj$mColumnLabels))
		rownames(results) <- jObj$mPatientIds
		colnames(results) <- jObj$mColumnLabels
	}
	results
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

getDataClinical <- function(theCombinedDir="/rsrch1/bcb/batcheffects/GENE_REPORT/combined", theVerboseFlag=FALSE)
{
	stopifnot(isValidDirectoryPath(theCombinedDir))
	stopifnot((TRUE==theVerboseFlag)||(FALSE==theVerboseFlag))
	getDataClinical_internal(theCombinedDir, theVerboseFlag=theVerboseFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
