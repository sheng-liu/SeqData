# AllClasses
# 
# 
###############################################################################
##' @exportClass SeqData
##' @import IRanges
##' @import GenomicRanges

# library needed:
library(methods)
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)


# Class
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
        coverageViews="SimpleRleViewsList"  
    ),
    prototype=prototype(
        chrSize=integer(0),
        bamFile=character(0),
        annotationFile=character(0),
        libSize=integer(0),
        
        readAlignment=GAlignments(),
        readCoverage=RleList(),
        featureAnnotation=GRanges(),
        coverageViews=RleViewsList()
        
    )
)
