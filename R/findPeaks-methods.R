
## findPeaks-methods
## 
## 
###############################################################################
##' @name findPeaks
##' @aliases findPeaks
##' @title findPeaks
##' @rdname findPeaks-methods
##' @docType methods 
##' @description Find reads enriched regions using poission distribution model.
##' @usage findPeaks(bamFile.chip) # this roxygen directive does not working
##' @method findPeaks # this roxygen directive does not working
##' @param bamFie.chip Full path to chip bam file.
##' @param bamFile.control Full path to control bam file.
##' @param obj.chip SeqData object containing information on chip.
##' @param obj.control SeqData object containing information on control.
##' @param description Optional. User can input short desciption for the output file, it will show up in the file name of the output file. 
##' @return
##' \itemize{
##' \item{peakTable.csv} A csv file contaning information about the significantly enriched regions, including chromosome,start,end,peakPositions,peakCoverageSums,rpkm.peakCoverageSums,chipPeakHeights,controlPeakHeights,p.values,fdr,foldChange.
##' \item{coverageView slot} A RleViews object containing Views of the peak regions. 
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
##' # to view coverageView slot 
##' # coverageView(obj.chip)[[1]]
##' @details
##' possibility of a reads been random
##' possibility of peak heights (base coverage) at specific location is esitmated from the average (lambda) coverage of the peak region, 1kb region, 5kb region and 10 kb regions. 
##' This is to account for local fluctuations.

###############################################################################


setMethod(
    f="findPeaks",
    signature=c(
        obj.chip="SeqData",obj.control="SeqData",
        bamFile.chip="character",bamFile.control="character"),
    definition=function(
        obj.chip=NULL,obj.control=NULL,
        bamFile.chip=character(0),bamFile.control=character(0),
        fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        
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
        
        # peakHeight as measure
        peaks.chip=viewMaxs(islands.chip)
        peaks.control=viewMaxs(islands.control)
        
        enriched=peaks.chip>=peaks.control
        # (added "=" to the equation so that when obj.chip and obj.control are the same, they still can be processed)
        
        coverageView(obj.chip)=islands.chip[enriched]
        coverageView(obj.control)=islands.control[enriched]
        
        ## ---------------------------------------------------------------------
        ## findPeaks 
        
        cat("Finding peak summits ... ","\n")
        
        # get the postion of the maximum peak height (summit)
        peakPositions=viewWhichMaxs(coverageView(obj.chip))
        
        ## center the flat peaks
        # viewWhichMaxs returns the leftmost positon for flat peaks, to center the positon of the flat peaks
        flatPeaks = width(viewRangeMaxs(coverageView(obj.chip)))>1
        flatPeaks.increment=floor(
            width(viewRangeMaxs(coverageView(obj.chip)))[flatPeaks])/2
        peakPositions[flatPeaks]= peakPositions[flatPeaks]+flatPeaks.increment
        
        ##-----normalize libSize difference-----##
        
        # normalize libSize/sequencing depth difference
        # by sampling the bigger library to match the samll library
        # then compute coverage
        
        chipSet=.samplingNorm(obj.chip,obj.control)
        obj.chip=chipSet$chip
        obj.control=chipSet$control
        
        ## add the scale up option in previous version
        
        ##-----estimate significance of peaks from the control-----##
        
        # use peakHeights as significance estimater
        peakHeights=viewMaxs(coverageView(obj.chip))
        
 
        # a bigger loop has to be vectorized into many smaller ones
        
        # setting up different size of windows for calcualting lambda
        views.control.1k=RleViewsList(sapply(chrl,function(chr){
            Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-500,
                end=peakPositions[[chr]]+500)
        }))
        
        views.control.5k=RleViewsList(sapply(chrl,function(chr){
            Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-2500,
                end=peakPositions[[chr]]+2500)
        }))

        views.control.10k=RleViewsList(sapply(chrl,function(chr){
            Views(
                readCoverage(obj.control)[[chr]],
                start=peakPositions[[chr]]-5000,
                end=peakPositions[[chr]]+5000)
        }))
        
        views.control.region=RleViewsList(sapply(chrl,function(chr){
            Views(
                readCoverage(obj.control)[[chr]],
                start=start(islands.control[[chr]]),
                end=end(islands.control[[chr]]))
        }))

        # calculate lambda using above views for each chromosome range  
        # lambda=the total reads/length, it is the average reads within these different windows 

        background.1k=viewSums(views.control.1k)/1000
        background.5k=viewSums(views.control.5k)/5000
        background.10k=viewSums(views.control.10k)/1e4
        
        ## background.region is  the coverage sum at that region divided by its length to get the average-lambda
        
        background.region=viewSums(views.control.region)/width(views.control.region)
        
        # this is not a deep list, it is a 24 value vector, one value for each chromosome
        background.chr=sapply(chrl,function(chr){
            sum(readCoverage(obj.control)[[chr]])/chrSize(obj.control)[chr]
        })   
        
        # to use pmax, this need to have the same structure as other backgrounds, as NumericList
        names(background.chr)=chrl
        
        background.chr=NumericList(sapply(chrl,function(chr){
            rep(background.chr[[chr]],length(background.1k[[chr]]))
            
        }))
        
        ## ---------------------------------------------------------------------
        ## calculate lambda and stats
        cat("local.lambda = ",local.lambda,"\n")
        # exclude background.region, background.500bp if obj.control is absent
        if (identical(readCoverage(obj.chip),readCoverage(obj.control))){
            
            if (local.lambda==F) {
                cat("lambda estimated from chr only\n")
                lambda=round(background.chr) 
            }else{
                cat("lambda estimated from chr,5k and 10k.\n")
                lambda=round(
                    pmax(background.chr,background.5k,background.10k)) 
            }
           
        }else{
            if (local.lambda==F) {
                cat("lambda estimated from chr only\n")
                lambda=round(background.chr) 
            }else{
                cat("lambda estimated from chr,[region,1k,] 5k and 10k.\n")
                lambda=round(
                    pmax(background.chr,background.region,
                         background.1k,background.5k,background.10k))
                
            }
            
        }


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
        
        ##-----change the coverageView to the final significantPeaks-----##
        
        # subseting coverageView views, so it only contains significant peaks
        peakIndex=(fdr<=fdr.max&foldChange>=foldChange.min)
        fdr=fdr[peakIndex]
        foldChange=foldChange[peakIndex]
        p.values=p.values[peakIndex]
        
        coverageView(obj.chip)=coverageView(obj.chip)[peakIndex]
        peakCoverageSums=viewSums(coverageView(obj.chip))
        peakPositions=viewWhichMaxs(coverageView(obj.chip))
        chipPeakHeights=viewMaxs(coverageView(obj.chip))
        controlPeakHeights=lambda[peakIndex]
        
        ## ---------------------------------------------------------------------
        ## Output 
        
        # make peaksTable
        peakList=sapply(chrl,function(chr){      
            chromosome=rep(chr,length(coverageView(obj.chip)[[chr]]))
            start=start(coverageView(obj.chip)[[chr]])
            end=end(coverageView(obj.chip)[[chr]])
            numBasesCovered=sum(width(coverageView(obj.chip)[[chr]]))
            peakPositions=peakPositions[[chr]]
            peakCoverageSums=peakCoverageSums[[chr]]
            # from coverage to rpkm, peakReadCountSums/qwidth
            peakReadCountSums=peakCoverageSums/mean(qwidth(readAlignment(obj.chip)))
            rpkm=.rpkm.libSize(
                peakReadCountSums,libSize(obj.chip),numBasesCovered)
            chipPeakHeights=chipPeakHeights[[chr]]
            controlPeakHeights=controlPeakHeights[[chr]]
            p.values=p.values[[chr]]
            fdr=fdr[[chr]]
            foldChange=foldChange[[chr]]
            cbind(
                chromosome,start,end,peakPositions,peakCoverageSums,
                rpkm,chipPeakHeights,controlPeakHeights,
                p.values,fdr,foldChange)
        })
        
        peakTable=do.call(rbind,peakList)
        
        # sort table based on fdr
        peakTable.fdr=data.frame(peakTable[order(peakTable[,10]),])
        
        #output table
        cat("Output peakTable...","\n")
        
        fileName=paste("peakTable-",description,"-",.timeStamp(bamFile(obj.chip)),".csv",sep="")
        
        
        
        
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
        fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        
        #.checkBamSlot
        # check whether bamFile slot is empty
        if (length(bamFile(obj.chip))==0||length(bamFile(obj.control))==0) {
            stop ("Slot bamFile is empty, please fill in bamFile slot first","\n")    
        }
        
        bamFile.chip=bamFile(obj.chip)
        bamFile.control=bamFile(obj.control)
        
        obj=findPeaks(obj.chip,obj.control,bamFile.chip,bamFile.control,fdr=fdr.max,foldChange=foldChange.min,description=description,local.lambda=local.lambda)
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
        obj.chip=NULL,fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        
        # if obj.control is missing, use obj.chip itself as its control
        obj.control=obj.chip
        obj=findPeaks(obj.chip,obj.control,fdr=fdr.max,foldChange=foldChange.min,description=description,local.lambda=local.lambda)
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
        fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        
        obj.chip=new("SeqData")
        obj.control=new("SeqData")
        
        obj.chip=getReadCoverage(obj.chip,bamFile.chip)
        obj.control=getReadCoverage(obj.control,bamFile.control)
        
        obj.chip=findPeaks(
            obj.chip,obj.control,fdr=fdr.max,foldChange=foldChange.min,description=description)
        
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
        bamFile.chip=character(0),fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        
        # if bamFile.control is missing, use bamFile.chip itself as its control
        bamFile.control=bamFile.chip
        obj=findPeaks(bamFile.chip=bamFile.chip,bamFile.control=bamFile.control,fdr=fdr.max,foldChange=foldChange.min,description=description,local.lambda=local.lambda)
        
        return (obj)
        
    })

##-----------------------------------------------------------------------------
## TODO:



# TODO: switch for either TFBS or Histone modification 
# TODO: add peakPositions.control for sample swap to remove negtive peaks, or just add a swapper which is better off. 
# TODO: add scale factor for libSizeDiff
# TODO: plot enriched.peaks vs unenriched (It's a nice pic, worth doing)

# TODO: when libSize has huge difference, this may leads to loss of information, further normalization using MAnorm package maybe added if needed.

## performance
## chr12 with H3 as control
# $sampling.time
# [1] 91.2

