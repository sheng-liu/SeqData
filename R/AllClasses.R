# Classes
# 
# 
###############################################################################

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
        libSize="integer",    

        readAlignment="GappedAlignments",
        readCoverage="RleList",
        featureAnnotation="GRanges",        
        peakRegions="SimpleRleViewsList"  
        
        #(rleViews,reads more than 1 and ranges longer than 100bp,which are informational, not noise, baseline analysis)
        # later may make it peakRegions, which is more useful
        # peaks are called within these regions
        # every enrichedRegions has a peak basically, it is sorted based on fdr in the end
        
        
        
        # remove countTable, peakTable slots
        # slots are great for passing variable between functions, countTable and PeakTable can be output and input as file, do not need to taken up a slot.
        
        
        
        
        
    ),
    prototype=prototype(
        chrSize=integer(0),
        bamFile=character(0),
        annotationFile="character",
        libSize=integer(0),
      
        readAlignment=GappedAlignments(),
        readCoverage=RleList(),
        featureAnnotation=GRanges(),
        peakRegions=RleViewsList()
        
    )
)
