# AllClasses
# 
# 
###############################################################################

##' @name SeqData
##' @aliases SeqData
##' @title SeqData-constructor and supported functions
##' @rdname SeqData-class
##' @docType methods
##' @description The SeqData function is a constructor for constructing SeqData instances. Four supporting functions are available to conduct basic sequence analysis upon SeqData.
##' @usage
##' SeqData(bamFile=character(0),annotationFile=character(0))
##' @param bamFie Full path to bam file (optional).
##' @param annotationFile Full path to annotation file (optional).
##' @return a SeqData object contains 8 slots with four supporting functions. 
##' @details 
##' This package provides a data structure and four supporting functions for sequencing data (eg.ChIPseq and RNAseq) analysis.
##' \itemize{
##' \item{SeqData} an S4 data structure, contains 8 slots:
##'  
##' bamFile="character",
##' 
##' annotationFile="character",
##' 
##' chrSize="integer",
##' 
##' libSize="numeric",  
##' 
##' readAlignment="GAlignments",
##' 
##' readCoverage="RleList",
##' 
##' featureAnnotation="GRanges",  
##'       
##' coverageView="SimpleRleViewsList"  
##' 
##'  \item{Supporting fucntions} Four supporting functions includes:
##'   
##' viewCoverage
##' 
##' countReads
##' 
##' findPeaks
##' 
##' annotatePeaks 
##' }
##' 
##' @section Usage:{SeqData(bamFile=character(0),annotationFile=character(0))} 
##' @examples
##' ## seq=SeqData(bamFile="bamFileLocation",annotationFile="annotationFileLocation")
##' @seealso
##' See corresponding function documentation for details.


##' @import IRanges
##' @import GenomicRanges
##' 
## dplyr has masked  intersect, setdiff, setequal, union from base and other packages, try to use importFrom instead of import package
##' @importFrom dplyr summarise group_by select
##' 
##' 
# library needed:
library(methods)
#library(IRanges)
#library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
# library(dplyr)


# .onLoad <- function(SeqData,Seqata) {
#     cat("message from .onLoad via cat\n")
#     message("message from .onLoad via message")
#     packageStartupMessage("message from .onLoad via
# packageStartupMessage\n")
# }

## .onLoad, .onUnload, .onAttach and .onDetach are looked for as internal objects in the namespace and should not be exported

.onLoad <- function(SeqData,Seqata) {
    packageStartupMessage("Loading package SeqData")
}

.onAttach <- function(SeqData,Seqata) {
    packageStartupMessage("Done")
}

# Class
##' @exportClass SeqData
setClass(
    Class="SeqData",
    representation=representation(
        
        bamFile="character",
        annotationFile="character",
        chrSize="integer",
        libSize="numeric",    
        
        readAlignment="GAlignments",
        readCoverage="RleList",
        featureAnnotation="GRanges",        
        coverageView="SimpleRleViewsList",
        baseAlignment="GRanges"
    ),
    prototype=prototype(
        chrSize=integer(0),
        bamFile=character(0),
        annotationFile=character(0),
        libSize=integer(0),
        
        readAlignment=GAlignments(),
        readCoverage=RleList(),
        featureAnnotation=GRanges(),
        coverageView=RleViewsList(),
        baseAlignment=GRanges()
        
    )
)

# the SeqData constructor 
##' @export SeqData
SeqData=function(bamFile=character(0),annotationFile=character(0)){
    new("SeqData",bamFile=bamFile,annotationFile=annotationFile)
}








