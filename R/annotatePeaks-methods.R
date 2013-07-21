## annotatePeaks


###############################################################################
##' @name annotatePeaks 
##' @aliases annotatePeaks
##' @title annotatePeaks
##' @rdname annotatePeaks-method
##' @docType methods
##' @description Annotate peaks with internal or user provided annotation file. 
##' @details It answers what gene are enriched with this protein/histone mark. Back end it is assaying which peak position overlaps with tss position by using "nearest" function from GenomicRanges package. 
##' It can annotate not just peaks find by annotatePeaks, but whatever RleViews is stored in coverageViews.
##' @param obj An SeqData Object.
##' @param feature Including gene, erv.
##' @param peakTable Full path to results "peakTable.csv" from findPeaks.
##' @return a targetTable. 
##' @section Usage:{
##' annotatePeaks(obj=NULL,peakTable=character(0),annotationFile=character(0),feature=c("gene","erv"))
##' }
##' @examples
##' #annotatePeaks(obj, feature)
##' #annotatePeaks(obj,annotationFile)
##' #annotatePeaks(peakTable, feature)
##' #annotatePeaks(peakTable,annotationFile)

###############################################################################

# dispatch on the presence/absense of obj ("SeqData") and peakTable("character")
# not supposed to present at the same time, it's a either or case.
setMethod(
    f="annotatePeaks",
    signature=c(obj="SeqData",peakTable="missing"),
    definition=function(
        obj,annotationFile=character(0),feature=c("gene","erv")){ 
        
        # set featureAnnotation slot based on presence of annotationFile
        if (length(annotationFile)!=0){
            cat ("Set featureAnnotation slot using provided annotationFile...\n")            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
        }else{
            
            #reset annotationFile slot based on feature inputed
            if (feature=="gene"||feature=="tss"||feature=="exon") feature="gene"
            if (feature=="erv"||feature=="5LTR") feature="erv"
                       
            cat ("Getting featureAnnotation...\n")            
            obj=getFeatureAnnotation(obj,feature=feature)
        }
        
        # get featureGRanges        
        featureGRanges=resize(featureAnnotation(obj),fix="start",width=1)
        # by not using feature="tss", maintains featureAnnotation information intact
        
        # check whether coverageViews are non-empty
        if (length(coverageViews(obj))==0){
            stop ("annotatePeaks requires coverageViews slot to be non-empty, please run findPeaks first.\n") 
        }
        
        # get peakGRanges
        peakIRangesList=ranges(coverageViews(obj))  #SimpleRangesList class
        peakGRanges=as(peakIRangesList,"GRanges")
        
        # get targetTable
        targetTable=.getTargetTable(obj,peakGRanges=peakGRanges,featureGRanges=featureGRanges)
        
        
        cat("Output file...","\n")
        
        file.name=basename(bamFile(obj))
        file.name.csv=paste("targetTable.",format(Sys.time(),"%Y%m%d.%H%M%S"),file.name,".csv",sep="")
        
        write.csv(file=file.name.csv,targetTable)

        return(obj)
        
    })


## add a dispatch so the function can be used without obj but with the output of findPeaks(),peakTable
setMethod(
    f="annotatePeaks",
    signature=c(obj="missing",peakTable="character"),
    definition=function(peakTable=character(0),annotationFile=character(0),feature=c("gene","erv")){
        
        obj=new("SeqData")
        bamFile(obj)="subsetPeakTable"
        
        # set featureAnnotation slot based on presence of annotationFile
        if (length(annotationFile)!=0){
            cat ("Setting up featureAnnotation slot using provided annotationFile...\n")            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
        }else{
            # reset annotationFile slot based on feature inputed
            if (feature=="gene"||feature=="tss"||feature=="exon") feature="gene"
            if (feature=="erv"||feature=="5LTR") feature="erv"
            
            cat ("Getting featureAnnotation...\n")            
            obj=getFeatureAnnotation(obj,feature=feature)
        }
        
        # get featureGRanges
        featureGRanges=resize(featureAnnotation(obj),fix="start",width=1)
        #featureGRanges=featureAnnotation(obj)

        # get peakGRanges
        peakTable=read.csv(file=peakTable,as.is=T,header=T)

        peakIRanges=RangedData(
            IRanges(
                start=peakTable$start,
                end=peakTable$end),
            space=peakTable$chromosome)
        
        peakGRanges=as(peakIRanges,"GRanges")
        ## peakTable deals with coverage and views, it doesn't have strand information
        
        # get targetTable
        targetTable=.getTargetTable(obj,peakTable,peakGRanges,featureGRanges)
        
        cat("Output file...","\n")
        fileName=paste("targetTable-",.timeStamp(bamFile(obj)),sep="")
        write.csv(file=fileName,targetTable)
        
        return(obj)
        
    })




##-----------------------------------------------------------------------------
## helperFunctions

## a S3 method dispatching on the presence of peakTable
.getTargetTable=function(obj=NULL,peakTable=data.frame(),peakGRanges=GRanges(),featureGRanges=GRanges(), distance.max=4000){
        
        # calculate distance between peak and feature
        distanceToNearest=suppressWarnings(distanceToNearest(peakGRanges,featureGRanges))
        
        distanceToFeature=distanceToNearest$distance
        featureIndex=distanceToNearest$subjectHits
        #peakIndex=distanceToNearest$queryHits
        
        # subset featureIndex to generate target index 
        targetIndex=abs(distanceToFeature)<=distance.max
        targetGRanges=peakGRanges[targetIndex]
        
        ##-----peak information-----##
        # using targetIndex to get peak information, either subset on peakTable, or subset on coverageViews. Both method gives same results. 
        
        if(length(peakTable)!=0){     
        #if (exists("peakTable")) doesn't work as good.
            
            #subsetting peakTable to get peak information
            targetPeakPos=peakTable$peakPositions[targetIndex] 
            targetPeakReadCount=peakTable$peakReadCountSums[targetIndex]
            #width=peakTable$end-peakTable$start+1
            #numBasesCovered=sum(width)
            rpkm.targetPeakReadCount=peakTable$rpkm.peakReadCountSums[targetIndex]
            
        }else{
            

            #subseting peakRegions to get targetRegions

            chrl=seqlevels(peakGRanges)
            
            targetGRanges=peakGRanges[targetIndex]
            targetGRangesList=split(targetGRanges,seqnames(targetGRanges))
            
            targetRegions=RleViewsList(sapply(chrl,function(chr) {
                Views(
                    readCoverage(obj)[[chr]],
                    start=start(targetGRangesList)[[chr]],
                    end=end(targetGRangesList[[chr]]))
            }))

            #get peak information from targetRegions 
            targetReadCount=viewSums(targetRegions)
            targetPeakPos=viewWhichMaxs(targetRegions)
            
            targetPeakReadCount=unlist(targetReadCount)
            targetPeakPos=unlist(targetPeakPos)
            
            numBasesCovered=sum(width(coverageViews(obj)))
            # numBasesCovered=sum(width(targetRegions))
            # annotatePeak shouldn't change the rpkm value of the peak table
            # targetPeak are among all the peaks, the rpkm of each peak should be normalized to all peaks not just targetPeak group.

            rpkm.targetPeakReadCount=unlist(
                .rpkm.libSize(targetReadCount,libSize(obj),numBasesCovered))    
        }
        
        ##-----feature information-----##
        featureInfo=mcols(featureGRanges[featureIndex])
        targetFeatureNames=featureInfo$symbol[targetIndex]
        
        if (length(featureInfo$EntrezID)!=0){
            targetFeatureIDs=featureInfo$EntrezID[targetIndex]
        }
        targetFeaturePos=start(featureGRanges[featureIndex])[targetIndex] 
        
        ##-----peak-feature distance-----##
        peakDistanceToFeature=distanceToFeature[targetIndex]
        
        ##-----construct targetTable-----##
        # didn't put chr in, but it output from one of the variables.
        if (length(featureInfo$EntrezID)!=0){
            targetTable=as.matrix(cbind(targetFeatureNames,targetFeatureIDs,targetFeaturePos,targetPeakPos,peakDistanceToFeature,targetPeakReadCount,rpkm.targetPeakReadCount))
        }else{
            targetTable=as.matrix(cbind(targetFeatureNames,targetFeaturePos,targetPeakPos,peakDistanceToFeature,targetPeakReadCount,targetPeakReadCount,rpkm.targetPeakReadCount))
        }

        return(targetTable)
        
    }
    

##-----------------------------------------------------------------------------
## TODO:
## TODO: get the distance to tss and plot it.

        
        
        