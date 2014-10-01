## annotatePeaks-methods


###############################################################################
##' @name annotatePeaks 
##' @aliases annotatePeaks
##' @title annotatePeaks
##' @rdname annotatePeaks-methods
##' @docType methods
##' @description Annotate peaks with internal or user provided annotation file. 
##' @details It answers what gene are enriched with this protein/histone mark. Back end it is assaying which peak position overlaps with tss position by using "nearest" function from GenomicRanges package. 
##' It can annotate not just peaks find by annotatePeaks, but whatever RleViews is stored in coverageView.
##' @param obj An SeqData Object.
##' @param feature=c("gene","erv","repeats","transposon")
##' including all features that are available in getFeatureAnnotation function.
##' Note gene, exon, transcripts and tss are processed the same as gene.
##' @param peakTable Full path to results "peakTable.csv" from findPeaks.
##' @param distance.max maximum peak-to-feature distance used to locate features.
##' @param description Optional. User can input short desciption for the output file, it will show up in the file name of the output file. 
##' @return a targetTable. 
##' @section Usage:{
##' annotatePeaks( obj,annotationFile=character(0),feature=character(0),distance.max=4000,description=character(0))
##' }
##' @examples
##' #if you have setup bamFile and annotationFile slot
##' #annotatePeaks(obj)
##' 
##' #if you want to use internal annotation file (feature)
##' #annotatePeaks(obj, feature="gene")
##' 
##' # specify another annotation file
##' #annotatePeaks(obj,annotationFile)
##' 
##' # use peakTable output from findPeaks
##' #annotatePeaks(peakTable, feature)
##' 
##' # use peakTable and annotation file
##' #annotatePeaks(peakTable,annotationFile)

###############################################################################

# dispatch on the presence/absense of obj ("SeqData") and peakTable("character")
# not supposed to present at the same time, it's a either or case.
setMethod(
    f="annotatePeaks",
    signature=c(obj="SeqData",peakTable="missing"),
    definition=function(
        obj,annotationFile=character(0),feature=character(0),distance.max=4000,description=character(0)){ 
        
        
        # set featureAnnotation slot based on presence of annotationFile
        if (length(annotationFile)!=0 ){
            cat ("Set featureAnnotation slot using provided annotationFile...\n")            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
        
        }else if (length(feature)!=0){
            
            #reset annotationFile slot based on feature inputed
            if (feature=="gene"||feature=="tss"||feature=="exon"||feature=="transcript") feature="gene"
            if (feature=="erv"||feature=="5LTR") feature="erv"
            
            cat ("feature=",feature,"\n")
            cat ("Getting featureAnnotation...\n")            
            obj=getFeatureAnnotation(obj,feature=feature)
        
        }else if (length(annotationFile(obj))!=0) {
                cat("featureAnnotation slot is not empty, using the current featureAnnotation...\n")
                obj=obj  
        }else{
            stop("please specify feature of interest (gene,erv,repeats,transposon) or specify custom annotationFile.\n")
        }
        # get featureGRanges        
        featureGRanges=featureAnnotation(obj)
        # featureGRanges=resize(featureAnnotation(obj),fix="start",width=1)
        # by resize featureGRanges in the helper function, maintains featureGRanges information.
        # by not using feature="tss", maintains featureAnnotation information intact
        
        # check whether coverageView are non-empty
        if (length(coverageView(obj))==0){
            stop ("annotatePeaks requires coverageView slot to be non-empty, please run findPeaks first.\n") 
        }
        
        # get peakGRanges
        peakIRangesList=ranges(coverageView(obj))  #SimpleRangesList class
        peakGRanges=as(peakIRangesList,"GRanges")
        
        # get targetTable
        targetTable=.getTargetTable(obj,peakGRanges=peakGRanges,featureGRanges=featureGRanges,distance.max=distance.max)
        
        # output targetTable
        cat("Output file...","\n")
        fileName=paste("targetTable-",feature,"-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,targetTable)

        return(obj)
        
    })


## add a dispatch so the function can be used without obj but with the output of findPeaks(),peakTable
setMethod(
    f="annotatePeaks",
    signature=c(obj="missing",peakTable="character"),
    definition=function(peakTable=character(0),annotationFile=character(0),feature=character(0),distance.max=4000,description=character(0)){
        
        obj=new("SeqData")
        bamFile(obj)="subsetPeakTable"
        
        # set featureAnnotation slot based on presence of annotationFile
        if (length(annotationFile)!=0 ){
            cat ("Set featureAnnotation slot using provided annotationFile...\n")            
            obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
            
        }else if (length(feature)!=0){
            
            #reset annotationFile slot based on feature inputed
            if (feature=="gene"||feature=="tss"||feature=="exon") feature="gene"
            if (feature=="erv"||feature=="5LTR") feature="erv"
            
            cat ("feature=",feature,"\n")
            cat ("Getting featureAnnotation...\n")            
            obj=getFeatureAnnotation(obj,feature=feature)
            
        }else if (length(annotationFile(obj))!=0) {
            cat("featureAnnotation slot is not empty, using the current featureAnnotation...\n")
            obj=obj  
        }else{
            cat("please specify feature of interest (gene,erv,repeats,transposon), or specify custom annotationFile.\n")
        }
        
        # get featureGRanges
        featureGRanges=featureAnnotation(obj)
        #featureGRanges=resize(featureAnnotation(obj),fix="start",width=1)

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
        targetTable=.getTargetTable(obj,peakTable,peakGRanges=peakGRanges,featureGRanges=featureGRanges,distance.max=distance.max)
        
        cat("Output file...","\n")
        fileName=paste("targetTable-",feature,"-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,targetTable)
        
        return(obj)
        
    })




##-----------------------------------------------------------------------------
## helperFunctions

## a S3 method dispatching on the presence of peakTable
.getTargetTable=function(obj=NULL,peakTable=data.frame(),peakGRanges=GRanges(),featureGRanges=GRanges(),distance.max=distance.max){
    
    featureGRangesStart=resize(featureGRanges,fix="start",width=1)
        
    
    
        # calculate distance between peak and feature
        distanceToNearest=suppressWarnings(distanceToNearest(peakGRanges,featureGRangesStart))
        
        distanceToFeature=distanceToNearest@elementMetadata@listData$distance
        featureIndex=distanceToNearest@subjectHits
        featureIndex=featureIndex[is.na(featureIndex)==F]
        
        #peakIndex=distanceToNearest$queryHits
        
        # subset featureIndex to generate target index 
        targetIndex=abs(distanceToFeature)<=distance.max
        targetIndex=targetIndex[is.na(targetIndex)==F]
        
        #targetGRanges=peakGRanges[targetIndex]
        
        ##-----peak information-----##
        # using targetIndex to get peak information, either subset on peakTable, or subset on coverageView. Both method gives same results. 
        
        if(length(peakTable)!=0){     
        #if (exists("peakTable")) doesn't work as good.
            
            #subsetting peakTable to get peak information
            targetPeakPos=peakTable$peakPositions[targetIndex] 
            targetPeakCoverageSums=peakTable$peakCoverageSums[targetIndex]
            #width=peakTable$end-peakTable$start+1
            #numBasesCovered=sum(width)
            rpkm=peakTable$rpkm[targetIndex]
            
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
            targetCoverageSums=viewSums(targetRegions)
            targetPeakPos=viewWhichMaxs(targetRegions)
            
            targetPeakCoverageSums=unlist(targetCoverageSums)
            targetPeakPos=unlist(targetPeakPos)
            
            targetPeakWidth=width(targetGRanges)
            
            numBasesCovered=sum(width(coverageView(obj)))
            # numBasesCovered=sum(width(targetRegions))
            # annotatePeak shouldn't change the rpkm value of the peak table
            # targetPeak are among all the peaks, the rpkm of each peak should be normalized to all peaks not just targetPeak group, as they are the same peak as in peak table, just been annotated. 

            targetReadCountSums=targetCoverageSums/mean(qwidth(readAlignment(obj)))
            rpkm=unlist(
                .rpkm.libSize(targetReadCountSums,libSize(obj),numBasesCovered))    
        }
        
        ##-----feature information-----##
        featureInfo=mcols(featureGRangesStart[featureIndex])
        targetFeatureNames=featureInfo$symbol[targetIndex]
        
        if (length(featureInfo$EntrezID)!=0){
            targetFeatureIDs=featureInfo$EntrezID[targetIndex]
        }
        targetFeaturePos=start(featureGRangesStart[featureIndex])[targetIndex] 
        targetFeatureWidth=width(featureGRanges[featureIndex])[targetIndex]
        
        ##-----peak-feature distance-----##
        peakDistanceToFeature=distanceToFeature[targetIndex]
        
        ##-----construct targetTable-----##
        # didn't put chr in, but it output from one of the variables.
        if (length(featureInfo$EntrezID)!=0){
            targetTable=as.matrix(cbind(targetFeatureNames,targetFeatureIDs,targetFeaturePos,targetFeatureWidth,targetPeakPos,targetPeakWidth,peakDistanceToFeature,targetPeakCoverageSums,rpkm))
        }else{
            targetTable=as.matrix(cbind(targetFeatureNames,targetFeaturePos,targetFeatureWidth,targetPeakPos,targetPeakWidth,peakDistanceToFeature,targetPeakCoverageSums,rpkm))
        }

        return(targetTable)
        
    }
    

##-----------------------------------------------------------------------------
## TODO:
## TODO: get the distance to tss and plot it.

        
        
        