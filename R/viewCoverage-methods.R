## viewCoverage-methods
## 
## 
###############################################################################

##' @name viewCoverage
##' @aliases viewCoverage
##' @title viewCoverage
##' @rdname viewCoverage-methods
##' @docType methods
##' @description Generate a bigWig file for viewing read coverage in a genome browser, or a coverageTable.csv for genomic ranges specified. 
##' @usage #not working
##' viewCoverage(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,smooth=T,description=character(0))

##' @param obj A SeqData object.
##' @param obj.control  Optional. A SeqData object used as a control to see differential coverage with obj.
##' @param bamFile The full path to the sequencing data. The data needs to be in bam file format.
##' @param bamFile.control Optional. A bam file used as a control to see differential coverage with bam file. 
##' @param weight.control A parameter used to weight control data (obj.control or bamFile.control). 
##' @param feature interal feature annotation data set, includes: feature=c("gene","erv","repeats","transposon")
##' @param ranges GRanges object specifying the ranges the coverageTable will produce.
##' @param annotationFile csv annotation file for specifying ranges of interest. It must includes the miminum five field: 'chrom','txSart','txEnd','strand' and 'name'. 
##' @param split.metaData information in GRanges,or annotation file other than 'chrom','txSart','txEnd','strand' are treated as metaData. if split.metaData=T, the ranges will be splited based on the first column of metaData, and the coverage will be calculated as the total coverage on all ranges that have the same metaData. It can be used in conjunction with parameter "feature","ranges" and "annotationFile". 
##' @param bigWig a switch for output bigWig file. The data is from readCoverage slot of SeqData.     
##' @param description Optional. User can input short desciption for the output file, it will show up in the file name of the output file. 
##' @details Generate bigWig file from a bam file for viewing read coverage in a genome browser; it can also generate coverage views on any ranges specified, output as a coverageTable.csv file. 
##' coverageTable containing read maxPosition, maxHeight, coverageSums at specific genomic ranges. 
##' It also set readCoverage slot, and coverageView slot if ranges are provided.
##' 
##' Note, parameter 'feature','ranges','annotationFile' or 'slit.metaData' are used only with single SeqData/bamFile.
##' @note Parameter 'feature','ranges','annotationFile' or 'slit.metaData' are used only with single SeqData/bamFile.
##' @return 
##' \itemize{
##' \item{bamFileName.bigWig} a bigWig file for viewing the coverage of the bamFile provided
##' \item{bamFileName.csv} a csv file contains the coverage information about the GRanges provided.
##' \item{coverageView slot} if ranges are passed in, a RleViews object are returned and coverage corresponding to each of those ranges are stored in the coverageView slot.
##' }
##' @section Usage:{
##' viewCoverage(obj=NULL,obj.control=NULL,
##' 
##' bamFile=character(0),bamFile.control=character(0),weight.control=1,
##' 
##' feature=character(0),ranges=GRanges(),annotationFile=character(0),
##' 
##' split.metaData=F,bigWig=F,smooth=T,description=character(0))}
##' @examples
##' # viewCoverage(obj)  
##' # viewCoverage(obj,ranges=GRanges)
##' # viewCoverage(obj,annotationFile=annotationFile) 
##' 
##' # viewCoverage(bamFile=bamFile) 
##' # viewCoverage(bamFile=bamFile,ranges=GRanges)
##' # viewCoverage(bamFile=bamFile,annotationFile=annotationFile)
##'  
##' # viewCoverage(obj,obj.control) 
##' # viewCoverage(bamFile,bamFile.control)
##' # to view coverageView slot 
##' # coverageView(obj)[[1]]
##' # viewCoverage(obj,feature="transposon",split.metaData=T,description="obj_transposon")

# the advantage of using "ranges" over "feature" is GRanges are more frequently used and are easy to assemble. while features are easy to use, but if change it to features, the program lacks a big input gate--GRanges. featureAnnotation(obj) is GRanges as well. It is for the ease of use I didn't put "feature" here, as it blocks a more general use --GRanges.


## ranges can be GRanges or GRangesList
## it must contain start, end, name (in mcols) etc.
## split.metaData=T, means it will split on metaData

## ranges need to always be GRanges, can't pass in GRangesList, but you can choose split.metaData=T, it will split on first column of meta data to make an internal GRangesList



##-----------------------------------------------------------------------------
## Methods:

setMethod(
    f="viewCoverage",
    signature=c(obj="SeqData",obj.control="missing",bamFile="missing",bamFile.control="missing"),
    definition=function(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,smooth=T,description=character(0)){
        
        print("Compute read coverage.")
        
        # check whether readCoverage slot is non-empty
        if (length(readCoverage(obj))==0){                   
            obj=getReadCoverage(obj) 
            # Note bamFile slot non-empty check is done in getReadCoverage
        }
        
        ## ---------------------------------------------------------------------
        ## using coverageView to generate the views on coverage
        
        ## if annotationFile are supplied       
        if (length(annotationFile)!=0){            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
            ranges=featureAnnotation(obj) 
            if (split.metaData==T){   
                rangesList=split(ranges,mcols(ranges)[[1]])
            }       
        }else {
            ## if annotationFile are not supplied
            ## use internal featureAnnotation
            if(length(featureAnnotation(obj))==0){
                obj=getFeatureAnnotation(obj)
                ranges=featureAnnotation(obj)
                
            }else ranges=featureAnnotation(obj)
             }
        
        
        
        
        ## if ranges are supplied
        if (length(ranges)!=0 && split.metaData==T){            
            rangesList=split(ranges,mcols(ranges)[[1]])
            # this step can be saved as there is no use of rangesList but rangesList.sum
        }
        
        
        ## if feature are supplied
        if (length(feature)!=0){
            obj=getFeatureAnnotation(obj,feature=feature)
            ranges=featureAnnotation(obj)           
            if (split.metaData==T){
                rangesList=split(ranges,mcols(ranges)[[1]])
            }
        }
        
        ## calculate and dispatch viewCoverage on Ganges and GRangesList 
        ## TODO: another easier way of doing this maybe just add a step to make GRanges a single element GRangeslist, so all of them are GRangesList
        if (length(ranges)!=0){
            
            covView=.covView(readCoverage(obj),ranges)
            #set coverageView slot
            coverageView(obj)=covView$coverageView
            #reset ranges after intersect by .viewCov
            ranges.ints=covView$ranges
            
            # approximate readCounts from coverage, by coverageSums/qwidth
            readLength=if(length(readAlignment(obj))!=0)
                mean(qwidth(readAlignment(obj))) else 100
            
            
            
            #if (exists("rangesList")){   
            if (split.metaData==T){
                cat("split.metaData=T, summarize coverageView based on metaData...\n")
                rangesList.ints=split(ranges.ints,mcols(ranges.ints)[[1]])
                
                rangesList.sums=split(unlist(viewSums(coverageView(obj))),mcols(ranges.ints)[[1]])
                
                coverageSums=lapply(rangesList.sums,sum)            
                coverageSums=as.data.frame(do.call(rbind,coverageSums)) 
                names(coverageSums)="coverageSums"
                
                coverageMeans=lapply(rangesList.sums,mean)            
                coverageMeans=as.data.frame(do.call(rbind,coverageMeans)) 
                names(coverageMeans)="coverageMeans"
                
                
                # number of bases covered for each element in the rangesList
                rangesList.width=split(unlist(width(coverageView(obj))),mcols(ranges.ints)[[1]])
                
                numBasesCovered=lapply(rangesList.width,sum)
                numBasesCovered=as.data.frame(do.call(rbind,numBasesCovered))
                
                readCountSums=coverageSums/readLength
                readCountMeans=coverageMeans/readLength
                
                #normalize to libSize, for in library comparison
                #rpkm=.rpkm.libSize(readCountSums,libSize(obj),numBasesCovered)
                
                #normalize to total counts, given it is total exon reads counts or total repeats read counts
                rpkm.cov.sum=if(length(libSize(obj))!=0) 
                    .rpkm.libSize(readCountSums,libSize(obj),numBasesCovered) else
                        .rpkm(readCountSums,numBasesCovered)
                
                rpkm.cov.mean=if(length(libSize(obj))!=0) 
                    .rpkm.libSize(readCountMeans,libSize(obj),numBasesCovered) else
                        .rpkm(readCountMeans,numBasesCovered)
                
                #rpkm=.rpkm(readCountSums,numBasesCovered)
                #names(rpkm)="rpkm"
                
                totalRangeNums=elementLengths(rangesList.ints)
                coverageTable=as.matrix(cbind(coverageSums,rpkm.cov.sum,
                                              coverageMeans,rpkm.cov.mean,
                                              totalRangeNums))
                
            }else{
                coverageSums=unlist(viewSums(coverageView(obj)))
                coverageMeans=unlist(viewMeans(coverageView(obj)))
                maxHeight=unlist(viewMaxs(coverageView(obj)))
                maxPosition=unlist(viewWhichMaxs(coverageView(obj)))
                
                chromosome=as.vector(seqnames(ranges.ints))
                start=start(ranges.ints)
                end=end(ranges.ints)
                metaData=as.data.frame(mcols(ranges.ints))
                
                # rpkm normalization of coverageSums
                # number of bases covered for each range in ranges
                numBasesCovered=unlist(width(ranges(coverageView(obj))))
                
                readCountSums=coverageSums/readLength
                readCountMeans=coverageMeans/readLength
                
                rpkm.cov.sum=if(length(libSize(obj))!=0) 
                    .rpkm.libSize(readCountSums,libSize(obj),numBasesCovered) else
                        .rpkm(readCountSums,numBasesCovered)
                
                rpkm.cov.mean=if(length(libSize(obj))!=0) 
                    .rpkm.libSize(readCountMeans,libSize(obj),numBasesCovered) else
                        .rpkm(readCountMeans,numBasesCovered)
                
                coverageTable=as.matrix(cbind(chromosome,start,end,maxPosition,maxHeight,coverageSums,rpkm.cov.sum,coverageMeans,rpkm.cov.mean,metaData))           
            }
            
            # }
            
            ## output mandatory csv file
            #if (exists("coverageTable")) {
                cat("Output coverageView ...","\n")
                fileName=paste("coverageTable-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
                write.csv(file=fileName,coverageTable)
            #}
            
        }
        
        ## output bigWig file
        if (bigWig==T) {
            if (smooth==T){
                # smoothing with runing window average
                cat("Smoothing readCoverage...","\n")
                smoothedCoverage=runmean(readCoverage(obj),k=51,endrule="constant") 
                
                cat("Output bigWig file...","\n")
                fileName=paste("bigWig-",description,"-",.timeStamp(bamFile(obj)),".bigWig", sep="")
                export(smoothedCoverage,fileName)
                
            }else{
                cat("Output bigWig file...","\n")
                fileName=paste("bigWig-",description,"-",.timeStamp(bamFile(obj)),".bigWig", sep="")
                export(readCoverage(obj),fileName)
            }
            
            
        }
        
        cat("Done.","\n\n")
        return(obj)
        
    })


# a dispatcher when input is bamFile
setMethod(
    f="viewCoverage",
    signature=c(obj="missing",obj.control="missing",bamFile="character",bamFile.control="missing"),
    definition=function(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,smooth=T,description=character(0)) {
        
        obj=new("SeqData",bamFile=bamFile)               
        obj=viewCoverage(obj,bigWig=bigWig,description=description)         
        return(obj)
    })

##-----------------------------------------------------------------------------
## Compare two coverage 

## use one coverage as background control to generate differential coverage removed of background

setMethod(
    f="viewCoverage",
    signature=c(obj="SeqData",obj.control="SeqData", bamFile="missing",bamFile.control="missing"),
    definition=function(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,smooth=T,description=character(0)){
        
        # exclude input of feature,ranges,annotationFile or split.metaData
        if (length(feature)!=0 || length(ranges)!=0 || length(annotationFile)!=0||split.metaData==T){
            warning("Parameter 'feature','ranges','annotationFile' or 'slit.metaData' are used only with single SeqData/bamFile.\n")
        }
        
        # check whether either of obj's bamFile slot is empty
        if (length(bamFile(obj))==0||length(bamFile(obj))==0){
            stop("bamFile slot is empty, please fill in bamFile slot first\n")
        }
        
        # check whether readAlignment slot is empty
        if (length(readAlignment(obj))==0) obj=getReadAlignment(obj)
        if (length(readAlignment(obj.control))==0) obj.control=getReadAlignment(obj.control)
        
        # normalize libSize/sequencing depth difference
        # by sampling the bigger library to match the samll library
        
        chipSet=.samplingNorm(obj,obj.control)
        obj=chipSet$chip
        obj.control=chipSet$control
        
        # compute differential coverage
        if (weight.control!=1) {
            cat("Computing differential coverage, weight.control=",weight.control,"\n")
        }else{
            cat("Computing differential coverage\n")
        }
        
        
        diffCov=readCoverage(obj)-readCoverage(obj.control)*weight.control
        # identical to
        # diffCov=readCoverage(obj)-readCoverage(obj.control,weight=weight.control)
        
        # set obj coverage slot to diffCov
        readCoverage(obj)=diffCov
        
        # output bigWig file
        obj=viewCoverage(obj,bigWig=bigWig,description=description)
        return(obj)
    })

# a dispatcher for viewCoverage(bamFile,bamFile.control)
setMethod(
    f="viewCoverage",
    signature=c(obj="missing",obj.control="missing", bamFile="character",bamFile.control="character"),
    definition=function(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,smooth=T,description=character(0)){
        
        obj=new("SeqData",bamFile=bamFile)
        obj.control=new("SeqData",bamFile=bamFile.control)
        
        obj=viewCoverage(obj,obj.control,weight.control=weight.control,bigWig=bigWig,description=description)
        return(obj)
        
    })


##-----------------------------------------------------------------------------
## 
## @export .covView  # getCoverageView or coverageView should be a function output to user as it creates a slot

## if ranges are supplied or annotationFile are supplied
.covView=function(readCoverage,ranges){
    
    cat("Computing coverageView on ranges\n")
    # matching chr names and chr number in the readAlignment
    # and featureAnnotation "GRanges"
    intersect=intersectChr(readCoverage,ranges)        
    readCoverage=intersect$reads        
    ranges=intersect$features
    
    #create coverageView, "view on coverage at ranges"
    chrl=seqlevels(ranges)        
    rangesList=split(ranges,seqnames(ranges))        
    coverageView=RleViewsList(sapply(chrl,function(chr){
        Views(readCoverage[[chr]],
              start=start(rangesList)[[chr]],
              end=end(rangesList)[[chr]])
    }))
    
    covView=list(coverageView=coverageView,ranges=ranges)
    return(covView)
}


##-----------------------------------------------------------------------------
## 
## TODO:

# add a summmary of percentage of coverage of the genome
# slice /sum(seqlengths())

## TODO: another easier way of doing this maybe just add a step to make GRanges a single element GRangeslist, so all of them are GRangesList

