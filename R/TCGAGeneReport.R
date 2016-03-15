#TCGAGeneReport Copyright 2014, 2015, 2016 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(rJava)
library(Cairo)
library(gplots)
library(methods)

compressIntoFilename<-function(theString)
{
	### listing whole list of characters out looks wrong, but is locale independent
	theString <- gsub("[^ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789|_/\\/.-]", "", theString)
	theString <- gsub("\\", "_", theString, fixed=TRUE)
	theString <- gsub("/", "_", theString, fixed=TRUE)
	return(theString)
}

padChromosomeName <- function(theChr)
{
	if(1==nchar(theChr))
	{
		if (("x"!=theChr)&&("X"!=theChr)&&("y"!=theChr)&&("Y"!=theChr))
		{
			theChr <- paste("0", theChr, sep="")
		}
	}
	theChr
}

initGeneReport <- function(theParameters="-Xms2400m")
{
	stopifnot(is.character(theParameters))
	message("TCGAGeneReport 2016-03-11-1300")
	message("note: parameters must be in terms of 'm' or .jinit will fail on Linux")
	message("In testing, -Xmx42000m worked, while larger values caused core exceptions")
	myJavaJars <- file.path(
		system.file("TcgaGSData", "commons-codec-1.9.jar", package="TCGAGeneReport"),
		system.file("TcgaGSData", "commons-compress-1.7.jar", package="TCGAGeneReport"),
		system.file("TcgaGSData", "commons-io-2.4.jar", package="TCGAGeneReport"),
		system.file("TcgaGSData", "commons-net-3.3.jar", package="TCGAGeneReport"),
		system.file("TcgaGSData", "TcgaGSData.jar", package="TCGAGeneReport"),
		system.file("TcgaGSData", "TcgaIdConverter.jar", package="TCGAGeneReport"),
		fsep=.Platform$path.sep)
	message("Calling .jinit theParameters=", theParameters)
	message(myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters=theParameters)
	message("After .jinit")
	cfr <- .jnew("org/mda/bcb/tcgagsdata/CallFromR", "")
	.jcall(cfr, returnSig = "S", method="printVersion")
}

verboseMessage <- function(..., theVerboseFlag=FALSE)
{
	if(isTRUE(theVerboseFlag))
	{
		message(...)
	}
}

setJavaVerboseFlag <- function(theVerboseFlag)
{
	objFlag <- .jnew("java/lang/Boolean", theVerboseFlag)
	.jcall("org/mda/bcb/tcgagsdata/TcgaGSData", returnSig = "V", method="setVerboseFlag", objFlag)
}

matrixWithIssues<-function(...)
{
	warnLevel<-getOption("warn")
	on.exit(options(warn=warnLevel))
	options(warn=3) # warnings are errors
	return(matrix(...))
}

isValidDirectoryPath <- function(thePath)
{
	marker <- file.info(thePath)[1,"isdir"]
	if (is.na(marker))
	{
		dir.create(thePath, recursive = TRUE, showWarnings=FALSE)
		marker <- file.info(thePath)[1,"isdir"]
		unlink(thePath)
	}
	marker
}

#################################################################
#################################################################
#################################################################
#################################################################

getDataVersion <- function(theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip")
{
	stopifnot(file.exists(theZipFile))
	zippedDataStream <- unz(theZipFile, "time.txt")
	read.csv(zippedDataStream, header=FALSE, stringsAsFactors=FALSE)[[1]]
}

#################################################################
#################################################################
#################################################################
#################################################################

geneSymbolsForDataset <- function(theRequestVector, theGeneSymbolVector)
{
	matches <- NULL
	unmatched <- NULL
	for(theGene in theRequestVector)
	{
		##message("check ", theGene)
		if (theGene %in% theGeneSymbolVector)
		{
			message("matches 1-1 ", theGene)
			matches <- c(matches, theGene)
		}
		else
		{
			matchCheck <- grepl(paste("^", theGene, "\\|.*$", sep=""), theGeneSymbolVector, perl=TRUE)
			if (sum(matchCheck)>0)
			{
				##message("matches 1-1 ", theGene, " and ", packageList[matchCheck])
				matches <- c(matches, theGeneSymbolVector[matchCheck])
			}
			else
			{
				##message("unmatched ", theGene)
				unmatched <- c(unmatched, theGene)
			}
		}
	}
	results <- list(matches, unmatched)
	names(results) <- c("matched", "unmatched")
	results
}

geneSymbolsMatchDataset <- function(theRequestVector, theGeneSymbolVector)
{
	matches <- (0==sum((theRequestVector %in% theGeneSymbolVector)==FALSE))
	matches
}

getSynonyms <- function(theId, theZipFile="/rsrch1/bcb/batcheffects/GENE_REPORT/GeneSurvey.zip", theVerboseFlag=FALSE)
{
	setJavaVerboseFlag(theVerboseFlag)
	verboseMessage("getSynonyms theId=", theId, theVerboseFlag=theVerboseFlag)
	results <- NULL
	jReadGeneObj <- .jnew("org/mda/bcb/tcgagsdata/CallFromR", theZipFile)
	results <- .jcall(jReadGeneObj, returnSig = "[S", method="getList_GeneSymbol_Synonym", 
										.jnew("java/lang/String",theId))
	results
}


#################################################################
#################################################################
#################################################################
#################################################################
