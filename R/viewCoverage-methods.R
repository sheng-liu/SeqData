## viewCoverage-methods
## 
## 
###############################################################################

##' @name viewCoverage 
##' @aliases viewCoverage
##' @title viewCoverage
##' @rdname viewCoverage-method
##' @docType methods
##' @description Generate a bigWig file for viewing read coverage in a genome browser, or a coverageTable.csv for genomic ranges specified. 
##' @usage #not working
##' viewCoverage(obj,bamFile)
##' viewCoverage(bamFile=bamFile) 
##' viewCoverage(bamFile=bamFile,ranges=GRanges) 
##' viewCoverage(obj)          
##' @param obj An SeqData object 
##' @param bamFile The full path to the sequencing data. The data needs to be in bam file format.
##' @details Generate bigWig file from a bam file for viewing read coverage in a genome browser; it can also generate coverage views on any ranges specified, output as a coverageTable.csv file. 
##' coverageTable containing read maxPosition, maxHeight, readCountSums at specific genomic ranges. 
##' It also set readCoverage slot, and coverageViews slot if ranges are provided.
##' @return 
##' \itemize{
##' \item{bamFileName.bigWig} a bigWig file for viewing the coverage of the bamFile provided
##' \item{bamFileName.csv} a csv file contains the coverage information about the GRanges provided.
##' \item{coverageViews slot} if ranges are passed in, a RleViews object are returned and coverage corresponding to each of those ranges are stored in the coverageViews slot.
##' }
##' @section Usage:{viewCoverage(obj=NULL,bamFile=character(0),ranges=GRanges())}
##' @examples
##' # viewCoverage(obj)  
##' # viewCoverage(obj,bamFile) 
##' # viewCoverage(obj,ranges=GRanges)
##' # viewCoverage(bamFile=bamFile) 
##' # viewCoverage(bamFile=bamFile,ranges=GRanges)
##' 
##' # viewCoverage(obj,annotationFile=annotationFile) 
##' # viewCoverage(bamFile=bamFile,annotationFile=annotationFile) 


##-----------------------------------------------------------------------------
## Methods:

setMethod(
    f="viewCoverage",
    signature=c(obj="SeqData",bamFile="character",ranges="missing"),
    definition=function(obj,bamFile=character(1)){
        
        # check readCoverage slot
        if (length(readCoverage(obj))==0){                   
            obj=getReadCoverage(obj,bamFile)    
        }
        
        # smoothing with runing window average
        cat("Smoothing readCoverage...","\n")
        smoothedCoverage=runmean(readCoverage(obj),k=51,endrule="constant") 
        
        # export file for viewing
        cat("Output bigWig file...","\n")
        file.name=unlist(strsplit(bamFile,split="/"))
        file.name.bigWig=paste(file.name[length(file.name)],".bigWig",sep="")    
        export(smoothedCoverage,file.name.bigWig)
        
        cat("Done.","\n\n")
        return(obj)
        
    }
)

# a dispatcher when input is SeqData
setMethod(
    f="viewCoverage",
    signature=c(obj="SeqData",bamFile="missing",ranges="missing"),
    definition=function(obj){
        
        # check whether bamFile slot is empty
        if (length(bamFile(obj))==0) {
            stop ("Slot bamFile is empty, please fill in bamFile slot first","\n")    
        }
        
        # get coverage
        obj=viewCoverage(obj,bamFile(obj))        
        return(obj)
    })

# a dispatcher when input is bamFile
setMethod(
    f="viewCoverage",
    signature=c(obj="missing",bamFile="character",ranges="missing"),
    #signature=c(obj="missing",bamFile="character"),
    definition=function(bamFile){
        
        obj=new("SeqData")        
        bamFile(obj)=bamFile        
        obj=viewCoverage(obj)         
        return(obj)
        
    })


# a dispatcher when input is SeqData and GRanges
setMethod(
    f="viewCoverage",
    
    signature=c(obj="SeqData",bamFile="missing",ranges="GRanges"),
    definition=function(obj,ranges){
        
        
        # check readCoverage slot is non-empty
        
        
        if (length(readCoverage(obj))==0&&length(bamFile(obj))==0) {
            stop("readCoverage slot and bamFile slot are both empty, please fill in either slot first\n")
        }else{
            obj=getReadCoverage(obj)
        }
        
        # matching chr names and chr number in the readAlignment and featureAnnotation "GRanges"
        
        intersect=intersectChr(readCoverage(obj),ranges)        
        readCoverage=intersect$reads        
        ranges=intersect$features
        
        #create coverageView
        chrl=seqlevels(ranges)        
        rangesList=split(ranges,seqnames(ranges))        
        coverageView=RleViewsList(sapply(chrl,function(chr){
            Views(readCoverage[[chr]],
                  start=start(rangesList)[[chr]],
                  end=end(rangesList)[[chr]])
        }))
        
        #set coverageViews slot
        coverageViews(obj)=coverageView
        
        readCountSums=unlist(viewSums(coverageView))
        maxHeight=unlist(viewMaxs(coverageView))
        maxPosition=unlist(viewWhichMaxs(coverageView))
        
        start=start(ranges)
        end=end(ranges)
        
        #rpkm normalization of readCountSums
        numBasesCovered=sum(sum(width(ranges(coverageView))))
        readCountSums.rpkm=.rpkm.libSize(readCountSums,libSize(obj),numBasesCovered)
        coverageTable=as.matrix(cbind(start,end,maxPosition,maxHeight,readCountSums,readCountSums.rpkm))
        
        #output file
        cat("Output coverageView ...","\n")
        fileName=paste("coverageTable-",.timeStamp(bamFile(obj)),sep="")
        write.csv(file=fileName,coverageTable)
        
        return(obj)

    }) 

# a dispatcher when input is bamFile and GRanges
setMethod(
    f="viewCoverage",    
    signature=c(obj="missing",bamFile="character",ranges="GRanges"),
    definition=function(bamFile,ranges){
        
        obj=new("SeqData")        
        bamFile(obj)=bamFile        
        obj=getReadCoverage(obj)         
        obj=viewCoverage(obj,ranges=ranges)        
        return(obj)
        
    })     


# a dispatcher for 
#  viewCoverage(obj,annotationFile=annotationFile) 
setMethod(
    f="viewCoverage",
    signature=c(obj="SeqData",bamFile="missing",ranges="missing",annotationFile="character"),
    definition=function(obj,annotationFile){
        
        obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
        obj=viewCoverage(obj,ranges=featureAnnotation(obj))
        return(obj)
    })

# a dispatcher for 
#  viewCoverage(bamFile=bamFile,annotationFile=annotationFile) 
setMethod(
    f="viewCoverage",    
    signature=c(obj="missing",bamFile="character",ranges="missing",annotationFile="character"),
    definition=function(bamFile,annotationFile){
        obj=new("SeqData")        
        bamFile(obj)=bamFile        
        obj=getReadCoverage(obj)
        obj=viewCoverage(obj,annotationFile=annotationFile)
        return(obj)
    })
        






##-----------------------------------------------------------------------------
## TODO:

# add a summmary of percentage of coverage of the genome
# slice /sum(seqlengths())



