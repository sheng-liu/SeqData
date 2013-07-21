
## findPeaks-methods
## 
## 
###############################################################################
##' @name findPeaks 
##' @aliases findPeaks
##' @title findPeaks
##' @rdname findPeaks-method
##' @docType methods 
##' @description Find reads enriched regions using poission distribution model.
##' @usage  findPeaks(bamFile.chip)
##' @param bamFie.chip Full path to chip bam file.
##' @param bamFile.control Full path to control bam file.
##' @param obj.chip SeqData object containing information on chip.
##' @param obj.control SeqData object containing information on control.
##' @return
##' \itemize{
##' \item{peakTable.csv} A csv file contaning information about the significantly enriched regions, including chromosome,start,end,peakPositions,peakReadCountSums,rpkm.peakReadCountSums,chipPeakHeights,controlPeakHeights,p.values,fdr,foldChange.
##' \item{coverageViews slot} A RleViews object containing Views of the peak regions. 
##' }
##' @section Usage : {
##' findPeaks(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(0),bamFile.control=character(0),fdr.max=1e-5,foldChange.min=2)
##' }
##' @examples
##' #findPeaks(bamFile.chip)
##' #findPeaks(bamFile.chip,bamFile.control)
##' #findPeaks(obj.chip)
##' #findPeaks(obj.chip,obj.control)
##' #findPeaks(obj.chip,obj.control,bamFile.chip,bamFile.control)
###############################################################################


setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="SeqData",obj.control="SeqData",
        bamFile.chip="character",bamFile.control="character"),
    definition=function(
        obj.chip=NULL,obj.control=NULL,
        bamFile.chip=character(0),bamFile.control=character(0),
        fdr.max=1e-5,foldChange.min=2){
        
        ## ---------------------------------------------------------------------
        ## Assign enriched Regions   

        # set bamFile slot of the obj
        bamFile(obj.chip)=bamFile.chip
        bamFile(obj.control)=bamFile.control
        
        # check readCoverage slot
        if (length(readCoverage(obj.chip))==0||length(readCoverage(obj.control))==0){                   
            obj.chip=getReadCoverage(obj.chip,bamFile=bamFile.chip)
            obj.control=getReadCoverage(obj.control,bamFile=bamFile.control)
        }
        
        cat("Defining read enriched regions ...","\n")
        islands.chip=slice(readCoverage(obj.chip),lower=5)
        chrl=names(chrSize(obj.chip))
        islands.control=RleViewsList(sapply(chrl,function(chr) {
            Views(
                readCoverage(obj.control)[[chr]],
                start=start(islands.chip[[chr]]),
                end=end(islands.chip[[chr]]))
        }))
        # (this is what Views are, start, end and readCoverage)

        peaks.chip=viewMaxs(islands.chip)
        peaks.control=viewMaxs(islands.control)
        
        enriched=peaks.chip>=peaks.control
        # (added "=" to the equation so that when obj.chip and obj.control are the same, they still can be processed)
        
        coverageViews(obj.chip)=islands.chip[enriched]
        coverageViews(obj.control)=islands.control[enriched]

        ## ---------------------------------------------------------------------
        ## findPeaks 
        
        cat("Finding peak summits ... ","\n")
        
        # get the postion of the maximum peak height (summit)
        peakPositions=viewWhichMaxs(coverageViews(obj.chip))
        
        ## center the flat peaks
        # viewWhichMaxs returns the leftmost positon for flat peaks, to center the positon of the flat peaks
        flatPeaks = width(viewRangeMaxs(coverageViews(obj.chip)))>1
        flatPeaks.increment=floor(
            width(viewRangeMaxs(coverageViews(obj.chip)))[flatPeaks])/2
        peakPositions[flatPeaks]= peakPositions[flatPeaks]+flatPeaks.increment
        
        ##-----normalize libSize difference-----##
        
        # normalize libary size difference when control has more reads
        libSizeDiff=libSize(obj.chip)/libSize(obj.control)
        cat("library size ratio (chip/control) is",libSizeDiff,"\n")
        if (libSizeDiff<1) {
            cat("Sampling control data to normalize library size difference...","\n")
            libSize.normalize=min(libSize(obj.chip),libSize(obj.control))
            readAlignment(obj.chip)=sample(
                readAlignment(obj.chip),libSize.normalize)
            readAlignment(obj.control)=sample(
                readAlignment(obj.control),libSize.normalize)
            
        }

        # normalize libary size difference when chip has more reads
        peakHeights=viewMaxs(coverageViews(obj.chip))
        libSizeDiff=libSize(obj.chip)/libSize(obj.control)
        if (libSizeDiff>=1){
            cat("libSizeDiff=",libSizeDiff,"normalizeing peakHeight to libSize","\n")
            peakHeights=peakHeights/libSizeDiff
        }
        
        ##-----estimate significance of peaks from the control-----##
        
        # calculate lambda, the average
        lambda=NumericList(sapply(chrl,function(chr) {
            
            # setting up different size of windows for calcualting lambda
            views.control.1k=Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-500,
                end=peakPositions[[chr]]+500)
            
            views.control.5k=Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-2500,
                end=peakPositions[[chr]]+2500)
            
            views.control.10k=Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-5000,
                end=peakPositions[[chr]]+5000)
            
            # calculate lambda using above views for each chromosome           
            background.region=as.integer(
                readCoverage(obj.control)[[chr]][peakPositions[[chr]]])
            background.chr=sum(
                readCoverage(obj.control)[[chr]])/chrSize(obj.control)[chr]
            background.1k=viewSums(views.control.1k)/1000
            background.5k=viewSums(views.control.5k)/5000
            background.10k=viewSums(views.control.10k)/1e4

            # exclude background.region, background.500bp if obj.control is absent
            if (identical(readCoverage(obj.chip),readCoverage(obj.control))){
                lambda=round(
                    pmax(background.chr,background.5k,background.10k))                
            }else{
                lambda=round(
                    pmax(background.chr,background.region,
                         background.1k,background.5k,background.10k))
            }
            return(lambda)
            
        }))

        # calculate p-value, fdr and foldChange
        # (p-value, probability of observing x in the poisson distribution of the  background)  
        cat("Estimating significace of enrichement:","\n")
        cat("p-value, fdr, foldChange...","\n")
        
        p.values=NumericList(sapply(chrl,function(chr){
            ppois(peakHeights[[chr]],lambda[[chr]],lower.tail=FALSE)}))
        
        fdr=NumericList(sapply(chrl,function(chr){
            p.adjust(p.values[[chr]],method="fdr")}))       
        
        # calculate the fold change over the "average" in the control        
        foldChange <- (peakHeights +1)/(lambda+1)
        
        ##-----change the coverageViews to the final significantPeaks-----##

        # subseting coverageViews views, so it only contains significant peaks
        peakIndex=(fdr<=fdr.max&foldChange>=foldChange.min)
        fdr=fdr[peakIndex]
        foldChange=foldChange[peakIndex]
        p.values=p.values[peakIndex]

        coverageViews(obj.chip)=coverageViews(obj.chip)[peakIndex]
        peakReadCountSums=viewSums(coverageViews(obj.chip))
        peakPositions=viewWhichMaxs(coverageViews(obj.chip))
        chipPeakHeights=viewMaxs(coverageViews(obj.chip))
        controlPeakHeights=lambda[peakIndex]
        
        ## ---------------------------------------------------------------------
        ## Output 
        
        # make peaksTable
        peakList=sapply(chrl,function(chr){      
            chromosome=rep(chr,length(coverageViews(obj.chip)[[chr]]))
            start=start(coverageViews(obj.chip)[[chr]])
            end=end(coverageViews(obj.chip)[[chr]])
            numBasesCovered=sum(width(coverageViews(obj.chip)[[chr]]))
            peakPositions=peakPositions[[chr]]
            peakReadCountSums=peakReadCountSums[[chr]]
            rpkm.peakReadCountSums=.rpkm.libSize(
                peakReadCountSums,libSize(obj),numBasesCovered)
            chipPeakHeights=chipPeakHeights[[chr]]
            controlPeakHeights=controlPeakHeights[[chr]]
            p.values=p.values[[chr]]
            fdr=fdr[[chr]]
            foldChange=foldChange[[chr]]
            cbind(
                chromosome,start,end,peakPositions,peakReadCountSums,
                rpkm.peakReadCountSums,chipPeakHeights,controlPeakHeights,
                p.values,fdr,foldChange)
        })
        
        peakTable=do.call(rbind,peakList)
        
        # sort table based on fdr
        peakTable.fdr=data.frame(peakTable[order(peakTable[,10]),])
        
        #output table
        cat("Output peakTable...","\n")
        
        fileName=paste("peakTable-",.timeStamp(bamFile(obj)),sep="")
        write.csv(file=fileName,peakTable.fdr)
        
    return(obj.chip)

    })


## findPeaks(obj.chip,obj.control)
## a dispatcher when just pass in obj.chip and obj.control, given their bamFile slot are non-empty
setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="SeqData",obj.control="SeqData",
        bamFile.chip="missing",bamFile.control="missing"),
    definition=function(
        obj.chip=NULL,obj.control=NULL,
        fdr.max=1e-5,foldChange.min=2){
        
        #.checkBamSlot
        # check whether bamFile slot is empty
        if (length(bamFile(obj.chip))==0||length(bamFile(obj.control))==0) {
            stop ("Slot bamFile is empty, please fill in bamFile slot first","\n")    
        }
        
        bamFile.chip=bamFile(obj.chip)
        bamFile.control=bamFile(obj.control)
        
        obj=findPeaks(obj.chip,obj.control,bamFile.chip,bamFile.control,fdr=fdr.max,foldChange=foldChange.min)
        return(obj)

})

## findPeaks(obj.chip)
## a dispatcher when only obj.chip is available
setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="SeqData",obj.control="missing",
        bamFile.chip="missing",bamFile.control="missing"),
    definition=function(
        obj.chip=NULL,fdr.max=1e-5,foldChange.min=2){
      
      # if obj.control is missing, use obj.chip itself as its control
      obj.control=obj.chip
      obj=findPeaks(obj.chip,obj.control,fdr=fdr.max,foldChange=foldChange.min)
      return(obj)
        
    })


## findPeaks(bamFile.chip,bamFile.control)
## a dispatcher for passing in bamFile.chip and bamFile.control
setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="missing",obj.control="missing",
        bamFile.chip="character",bamFile.control="character"),
    definition=function(
        bamFile.chip=character(0),bamFile.control=character(0),
        fdr.max=1e-5,foldChange.min=2){
        
        obj.chip=new("SeqData")
        obj.control=new("SeqData")
        
        obj.chip=getReadCoverage(obj.chip,bamFile.chip)
        obj.control=getReadCoverage(obj.control,bamFile.control)
        
        obj.chip=findPeaks(
            obj.chip,obj.control,fdr=fdr.max,foldChange=foldChange.min)
        
        return(obj.chip)
       
    })


## findPeaks(bamFile.chip)
## a dispatcher for passing in just bamFile.chip
setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="missing",obj.control="missing",
        bamFile.chip="character", bamFile.control="missing"),
    definition=function(
        bamFile.chip=character(0),fdr.max=1e-5,foldChange.min=2){
        
        # if bamFile.control is missing, use bamFile.chip itself as its control
        bamFile.control=bamFile.chip
        obj=findPeaks(bamFile.chip=bamFile.chip,bamFile.control=bamFile.control,fdr=fdr.max,foldChange=foldChange.min)
        
        return (obj)

        })
        
##-----------------------------------------------------------------------------
## TODO:



# TODO: switch for either TFBS or Histone modification 
# TODO: add peakPositions.control for sample swap to remove negtive peaks, or just add a swapper which is better off. 
# TODO: add scale factor for libSizeDiff
# TODO: plot enriched.peaks vs unenriched (It's a nice pic, worth doing)
   

## performance
## chr12 with H3 as control
# $sampling.time
# [1] 91.2

