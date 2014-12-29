# countReads-methods
# 
# 
###############################################################################
##' @name countReads
##' @aliases countReads
##' @title countReads
##' @rdname countReads-methods
##' @docType methods
##' @description Count reads overlap different genomic features, such as gene,
##'   transciripts, exons, promoter regions (1kb upstream and downstream of
##'   tss), or any genomic ranges of interest. It will count how many reads is
##'   mapped to each genomic feature.
##' @usage countReads(obj,annotationFile,countingMode,interFeature,feature)
##'   countReads(bamFile,annotationFile,countingMode,feature)
##' @param obj A SeqData object.
##' @param annotationFile Optional full path to annotation file.
##' @param countingMode Including
##'   "Union","IntersectionStrict","IntersectionNotEmpty". See summarizeOverlaps
##'   from GenomicRanges package for details.
##' @param interFeature This parameter is accompanying countingMode. It is to
##'   decide, when multiple features overlap, and reads are mapped to that
##'   overlap regions, whether these reads should be removed using the
##'   countingMode. Default is TRUE, meaning using the countingMode to aovoid
##'   count a single reads more than once. If it is set to FALSE, read are
##'   counted to each feature they mapped to, which means they are allowed to
##'   count multiple times.
##' @param feature Genomic features including
##'   "exon","transcript","gene","erv","tss","5LTR" (#
##'   feature=c("exon","transcript","gene","erv","tss","5LTR"))
##' @param description Optional. User can input short desciption for the output
##'   file, it will show up in the file name of the output file.
##' @details This function answers how to assign reads to features that are
##'   overlapped. The advantage of this function is it does not count a read
##'   twice, and it take duplicated gene into count.The back end is the function
##'   summarizeOverlaps from GenomicRanges package.
##' @return A countTable includes feature name, read count and its rpkm
##'   normalization value.
##' @section Usage:
##'   {countReads(obj=NULL,bamFile=character(0),annotationFile=character(0), 
##'   countingMode="Union",interFeature=TRUE,feature=c("exon","transcript","gene","erv","repeats","transposon","tss","5LTR"))}
##'   
##' @examples 
##' #countReads(obj)  #readAlignment and featureAnnotation slots non-empty
##'                   #user can pass in GAlignments(readAlignment slot) and 
##'                   #GRanges(featureAnnotation slot) to countReads() via obj
##'                   
##' #countReads(obj,feature="gene")    #readAlignment slot non-empty
##' #countReads(obj,annotationFile=annotationFile)   #readAlignment slot non-empty
##' #countReads(bamFile=bamFile,annotationFile=annotationFile)
##' #countReads(bamFile=bamFile,feature="gene")
##' 

##-----------------------------------------------------------------------------
## Methods:

setMethod(
    f="countReads",
    signature=c(obj="SeqData",bamFile="missing"),
    definition=function(
        obj,annotationFile=character(0),countingMode="Union",interFeature=TRUE,
        feature=character(0),description=character(0)){
        
        annotationFile=annotationFile(obj)
        bamFile=bamFile(obj)
        
        print("Compute read counts")
        # set featureAnnotation slot, based on annotationFile
        if (length(annotationFile)!=0){
            cat ("Setting up featureAnnotation slot using provided annotationFile...\n")            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
            feature="custom"
        }else{
            if (length(feature)==0) feature="featureAnnotation"
            # length(feature)==0 designed for only pass in obj so user can define GRanges and GAlignments themselves
            
            cat("feature=",feature,"\n")
            switch(feature,
                   "exon"={obj=getFeatureAnnotation(obj,feature="exon")},
                   "gene"={obj=getFeatureAnnotation(obj,feature="exon")},
                   "transcript"={obj=getFeatureAnnotation(obj,feature="exon")},
                   "erv"={obj=getFeatureAnnotation(obj,feature="erv")},
                   "repeats"={obj=getFeatureAnnotation(obj,feature="repeats")},
                   "transposon"={obj=getFeatureAnnotation(obj,feature="transposon")},
                   "tss"={
                       obj=getFeatureAnnotation(obj,feature="tss")
                       tss.1kb=resize(
                           featureAnnotation(obj),fix="center",width=2000)
                       featureAnnotation(obj)=tss.1kb
                   },
                   "5LTR"={
                       obj=getFeatureAnnotation(obj,feature="5LTR")
                       LTR.1kb=resize(
                           featureAnnotation(obj),fix="center",width=2000)
                       featureAnnotation(obj)=LTR.1kb
                   },
                   "featureAnnotation"={
                       obj
                       annotationFile(obj)="featureAnnotation"
                   },
{stop(cat("No match was found for feature requested, please check spelling","\n"))}
            )
        }

        ##check readAlignment
        if(length(readAlignment(obj))==0){
            if(length(bamFile(obj))!=0){
                obj=getReadAlignment(obj)
            }else{
                stop("readAlignment slot is empty")
            }
        }
        
        # matching chr names and chr number in the readAlignment and featureAnnotation
        cat("Intersecting chromosomes between alignment and annotation...\n")
        intersect=intersectChr(readAlignment(obj),featureAnnotation(obj))
        intersectAlignment=intersect$reads
        intersectFeatureAnnot=intersect$features
        
        # process annotation base on feature specified then pass in summarizeOverlaps 
        cat("Processing annotation base on feature...\n")
        Annot=switch(feature,
                     "exon"={
                         exonAnnot=split(
                             intersectFeatureAnnot,
                             mcols(intersectFeatureAnnot)$exon)
                         
                     },
                     # exonic coverage of transcript
                     "transcript"={
                         transcriptAnnot=split(
                             intersectFeatureAnnot,
                             mcols(intersectFeatureAnnot)$transcript)
                     },
                     # exonic coverage of gene
                     "gene"={
                         geneAnnot=split(
                             intersectFeatureAnnot,
                             mcols(intersectFeatureAnnot)$gene) 
                     },
                     # custom provided annotation file
                     "custom"={
                         customAnnot=split(
                             intersectFeatureAnnot,
                             mcols(intersectFeatureAnnot)$symbol)},
                     # all the rest features
                     {Annot=split(
                         intersectFeatureAnnot,
                         mcols(intersectFeatureAnnot)$symbol)})
        
        #count reads over features
        cat("Counting reads overlapped with ",feature,"s"," ...\n",sep="")
        featureHits=suppressMessages(summarizeOverlaps(Annot,intersectAlignment,mode=countingMode,inter.feature=interFeature))
        counts=assays(featureHits)$counts
        
        #copyNum of each gene/exon/erv list
        
        #totalRangeNums=sapply(Annot,length) # takes much longer time
        totalRangeNums=elementLengths(Annot)
        #note copies in those chromosomes that are not specifically mapped to
        #the finished genome (chrN_random and chrUn_random) are not counted in
        #this way. (if set scanBamParam to include unmapped portion, those
        #chromosome will automatically included)
        
        #if want the total copy number in the genome
        #Annot.total=split(featureAnnotation(obj),
        #mcols(featureAnnotation(obj))$symbol)
        #totalRangeNums=sapply(Annot.total,length)

        #rpkm normalization
        numBasesCovered=sum(width(Annot))
        # this is a vecter contains each gene's numBasesCovered.
        

#         #equivalent to:
#         numBasesCovered=lapply(width(Annot),sum)
#         numBasesCovered=as.data.frame(do.call(rbind,numBasesCovered))
        
        
        #check libSize slot for rpkm normalization
        if (length(libSize(obj))==0){
            cat("counting libSize...\n")
            libSize(obj)=.libSize(bamFile(obj))
            cat("libSize=",libSize(obj),"\n")

        }
        
        #normalize to libSize, for in library comparison
        #rpkm=.rpkm.libSize(counts,libSize(obj),numBasesCovered)
        
        #normalize to total counts, given it is total exon reads counts or total repeats read counts
        rpkm=.rpkm(counts,numBasesCovered)
        countTable=cbind(counts,rpkm,totalRangeNums) # matrix
        colnames(countTable)=c("counts","rpkm","totalRangeNums")
        
        #output
        cat("output countTable...\n")
        
        fileName=paste("countTable-",feature,"-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,countTable)
        
        cat("\nDone!\n")
        
        return(obj)
        
    })


# a dispatcher when input is bamFile
setMethod(
    f="countReads",
    signature=c(obj="missing",bamFile="character"),
    definition=function(
        bamFile=character(0),annotationFile=character(0),
        countingMode="Union",interFeature=TRUE,feature=character(0),description=character(0)){
        
        obj=new("SeqData")
        bamFile(obj)=bamFile
        obj=getReadAlignment(obj,bamFile)
        
        obj=countReads(
            obj,annotationFile=annotationFile,
            countingMode=countingMode,interFeature=interFeature,feature=feature,description=description)

        return(obj)
        
    })


##-----------------------------------------------------------------------------
## 





