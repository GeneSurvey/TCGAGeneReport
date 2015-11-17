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

getOneToOne_List_internal <- function(theIdDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	results <- NULL
	jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/OneToOneUcscHgnc", theIdDataDir)
	results <- .jcall(jReadGeneObj, returnSig = "[S", method=theMethodString)
	results
}

getOneToOne_Name_internal <- function(theId, theIdDataDir, theMethodString, theVerboseFlag)
{
	setJavaVerboseFlag(theVerboseFlag)
	listResults <- as.vector(unlist(lapply(theId, function(myId)
	{
		jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/retrieve/OneToOneUcscHgnc", theIdDataDir)
		return(.jcall(jReadGeneObj, returnSig = "S", method=theMethodString, 
										.jnew("java/lang/String",myId)))
	})))
	names(listResults) <- theId
	listResults
}

#################################################################
#################################################################
# exported
#################################################################
#################################################################

getOneToOne_UCSC_List <- function(theIdDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", theVerboseFlag=FALSE)
{
	getOneToOne_List_internal(theIdDataDir, 'getOneToOne_UCSC_List', theVerboseFlag=theVerboseFlag)
}

getOneToOne_GeneSymbol_List <- function(theIdDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", theVerboseFlag=FALSE)
{
	getOneToOne_List_internal(theIdDataDir, 'getOneToOne_GeneSymbol_List', theVerboseFlag=theVerboseFlag)
}

getOneToOne_GeneSymbol_UCID <- function(theUCIDId, theIdDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", theVerboseFlag=FALSE)
{
	getOneToOne_Name_internal(theUCIDId, theIdDataDir, 'getOneToOne_GeneSymbol_UCID', theVerboseFlag=theVerboseFlag)
}

getOneToOne_UCID_GeneSymbol <- function(theGeneSymbol, theIdDataDir="/rsrch1/bcb/batcheffects/GENE_REPORT/iddata/downloads", theVerboseFlag=FALSE)
{
	getOneToOne_Name_internal(theGeneSymbol, theIdDataDir, 'getOneToOne_UCID_GeneSymbol', theVerboseFlag=theVerboseFlag)
}

#################################################################
#################################################################
#################################################################
#################################################################
