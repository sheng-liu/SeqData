
##' @title findPeaks
##' @name findPeaks
##' @description 
##' find reads enriched regions. 
##' @details poission distribution test

###############################################################################


setMethod(
    f="findPeaks",
    signature="SeqData",
    definition=function(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(1),bamFile.control=character(1)){
        
        
        ##----------------------------------------------------------------------
        ## Assign enriched Regions   
        ##
        # TODO: check if readCoverage is available
        
        ## because this function involves two obj, to avoid the complication, keep this code here. 
        
        cat("Defining read enriched regions ...","\n")
        
        #peakRegions(obj)=slice(readCoverage(obj),lower=5)
        
        
        islands.chip=slice(readCoverage(obj.chip),lower=5)
        
        if (exists("obj.control")){
            chrl=names(chrSize(obj.chip))
            
            islands.control=RleViewsList(sapply(chrl,function(chr) {
                Views(
                    readCoverage(obj.control)[[chr]],
                    start=start(islands.chip[[chr]]),
                    end=end(islands.chip[[chr]]))
            }))
            # this is what Views are, start, end and readCoverage
            
            
            peaks.chip=viewMaxs(islands.chip)
            peaks.control=viewMaxs(islands.control)
            
            enriched=peaks.chip>peaks.control
            
            peakRegions(obj.chip)=islands.chip[enriched]
            peakRegions(obj.control)=islands.control[enriched]
            
            
        }else{
            peakRegions(obj.chip)=islands.chip
        }
        
        # slicedRegions=slice(readCoverage(obj),lower=5)
        # peakRegions(obj)=slicedRegions[ranges(slicedRegions>100)]
        
        
        #obj.chip=getpeakRegions(obj.chip)
        
        #obj.chip=getpeakRegions(obj.chip=obj.chip,obj.control=obj.control) 
        #obj.control=getpeakRegions(obj.chip=obj.control)



        ##----------------------------------------------------------------------
        ##  findPeakSummit  
        ##       
        
        cat("Finding peak summits ... ","\n")
        
        
        ### TODO: switch for either TFBS or Histone modification 
        # now go with TFBS
        

        #         peaks.chip=viewMaxs(peakRegions(obj.chip))
        #         peaks.control=viewMaxs()
        #         
        #         enriched=peaks.chip>peaks.control
        #         
        #         peakRegions(obj.chip)=peaks.chip[enriched]
        #         peakRegions(obj.control)=peaks.control[enriched]
        
        ## TODO: add peakPositions.control for sample swap to remove negtive peaks, or just add a swapper which is better off. 
        
        ### findPeaks
        # get the postion of the maximum peak height (summit)
        
        peakPositions=viewWhichMaxs(peakRegions(obj.chip))
        
        ## center the flat peaks
        # viewWhichMaxs returns the leftmost positon for flat peaks, to center the positon of the flat peaks
        
        # create a logical to pick flat peaks
        flatPeaks = width(viewRangeMaxs(peakRegions(obj.chip)))>1
        
        flatPeaks.increment=floor(width(viewRangeMaxs(peakRegions(obj.chip)))[flatPeaks])/2
        
        peakPositions[flatPeaks]= peakPositions[flatPeaks]+flatPeaks.increment
        
        ### estimate significance of peaks from the control   
        
        # TODO: if (bamFile.control==NULL) {exlude 1k}
        if (!exists("obj.control")||!exists("bamFile.control")){
            bamFile.control=bamFile.chip
        } 
        
        # calculate lambda, the average
        lambda=NumericList(sapply(chrl,function(chr) {
            
            # setting up different size of windows for calcualting lambda
            views.control.200=Views(readCoverage(obj.control)[[chr]],start=peakPositions[[chr]]-100,end=peakPositions[[chr]]+100)
            
            views.control.1k=Views(readCoverage(obj.control)[[chr]],start=peakPositions[[chr]]-500,end=peakPositions[[chr]]+500)
            
            views.control.5k=Views(readCoverage(obj.control)[[chr]],start=peakPositions[[chr]]-2500,end=peakPositions[[chr]]+2500)
            
            views.control.10k=Views(readCoverage(obj.control)[[chr]],start=peakPositions[[chr]]-5000,end=peakPositions[[chr]]+5000)
            
            # calculate lambda using above views for each chromosome           
            background.peakPositions=as.integer(readCoverage(obj.control)[[chr]][peakPositions[[chr]]])
            
            background.chr=sum(readCoverage(obj.control)[[chr]])/chrSize(obj.control)[chr]
            
            background.200=viewSums(views.control.200)/200
            background.1k=viewSums(views.control.1k)/1000
            background.5k=viewSums(views.control.5k)/5000
            background.10k=viewSums(views.control.10k)/1e4
            
            lambda=round(pmax(background.peakPositions,background.chr,background.200,background.1k,background.5k,background.10k))
            
            return(lambda)
            
        }))
        
        # calculate p-value, fdr and foldChange
        # p-value, probability of observing x in the poisson distribution of the  background        
        
        
        
        
        
        # normalize libary size difference when chip has more reads
        
        peakHeights=viewMaxs(peakRegions(obj.chip))
        libSizeDiff=libSize(obj.chip)/libSize(obj.control)
        if (libSizeDiff>=1){
            peakHeights=peakHeights/libSizeDiff
            cat("libSizeDiff=",libsizeDiff,"normalizeing peakHeight to libSize","\n")
        }
        
        ## TODO: add scale factor for libSizeDiff
        
        cat("Estimating significace of enrichement:","\n")
        cat("p-value, fdr, fold change...","\n")
        
        p.values=NumericList(sapply(chrl,function(chr){ppois(peakHeights[[chr]],lambda[[chr]],lower.tail=FALSE)}))
        
        fdr=NumericList(sapply(chrl,function(chr){p.adjust(p.values[[chr]],method="fdr")}))       
        
        # calculate the fold change over the "average" in the control        
        foldChange <- (peakHeights +1)/(lambda+1)
        #         
##----------------------------------------------------------------------------
#        # output just peakTable
        
        
#         # make peaksTable
#         peakTable=do.call(cbind,sapply(chrl,function(chr){
#             data.frame(
#                 # TODO: add chromosome
#                 #index=seq(1:length(peakRegions(obj.chip)[[chr]]))
#                 chromosome=rep(chrl[which(chrl==chr)],length(peakHeights[[chr]])),
#                 #chromosome=chr==chrl,
#                 start=start(peakRegions(obj.chip)[[chr]]),
#                 end=end(peakRegions(obj.chip)[[chr]]),
#                 peakPositions=peakPositions[[chr]],
#                 chip=peakHeights[[chr]],
#                 control=round(lambda)[[chr]],
#                 p.values=p.values[[chr]],
#                 fdr=fdr[[chr]],
#                 foldChange=foldChange[[chr]]
#             )}))
#         
#         colnames(peakTable)=c("chromosome","starts","ends","peakPositions","chip","control","p.values","fdr","foldChange")
#         
# 
#         # peakTable.fdr=peakTable[order(peakTable[,8]),]
#         # somehow this gives an matrix
#         
#         peakTable.fdr=data.frame(peakTable[order(peakTable[,8]),])
#         
#         significantPeaks=subset(peakTable.fdr,fdr<=1e-3&foldChange>=3)
#         
#         cat("Output peakTable","\n")
#         write.csv(file="PeakTable.csv",significantPeaks)
#         
 ##----------------------------------------------------------------------------       
        # change the peakRegions to the final significantPeaks
        # or have findPeaks with user adjustable fdr setting ready to intake peakRegions as a quick start. no. 
        
        peakIndex=(fdr<=1e-3&foldChange>=3)
        #TODO: make fdr adjustable to user
        
        fdr=fdr[peakIndex]
        foldChange=foldChange[peakIndex]
        p.values=p.values[peakIndex]
        
        
        
        
        #subseting peakRegions views, so it only contains significant peaks, use it for calculating readCounts
        peakRegions(obj.chip)=peakRegions(obj.chip)[peakIndex]      
              
        peakReadCountSums=viewSums(peakRegions(obj.chip))
        
        peakPositions=viewWhichMaxs(peakRegions(obj.chip))
        
        chipPeakHeights=viewMaxs(peakRegions(obj.chip))

        controlPeakHeights=lambda[peakIndex]
        

        
        
        # make peaksTable
        peakTable=do.call(cbind,sapply(chrl,function(chr){
            data.frame(
                # TODO: add chromosome
                #index=seq(1:length(peakRegions(obj.chip)[[chr]]))
                #chromosome=rep(chrl[which(chrl==chr)],length(chipPeakHeights[[chr]])),
                chromosome=rep(chr,length(peakRegions(obj.chip)[[chr]])),
                #chromosome=chr==chrl,
                start=start(peakRegions(obj.chip)[[chr]]),
                end=end(peakRegions(obj.chip)[[chr]]),
                peakPositions[[chr]],
                peakReadCountSums[[chr]],
                chipPeakHeights[[chr]],
                controlPeakHeights[[chr]],
                p.values=p.values[[chr]],
                fdr=fdr[[chr]],
                foldChange=foldChange[[chr]]
                
            )}))
        
        colnames(peakTable)=c("chromosome","starts","ends","peakPositions","peakReadCountSums","chipPeakHeights","controlPeakHeights","p.values","fdr","foldChange")
       
        peakTable.fdr=data.frame(peakTable[order(peakTable[,9]),])
        
        cat("Output peakTable...","\n")
        write.csv(file="PeakTable.csv",peakTable.fdr)
        
        
       
        
        
        
        
        
##----------------------------------------------------------------------------         
        
        
        #TODO: plot enriched.peaks vs unenriched (It's a nice pic, worth doing)
        
        return(obj.chip)
        
        
    })

##------------------------------------------------------------------------------
##  findPeakSummit  dispatcher
##   

setMethod(
    f="findPeaks",
    signature="missing",
    definition=function(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(1),bamFile.control=character(1)){
        
        obj.chip=new("SeqData")
        obj.control=new("SeqData")
        
        obj.chip=getReadCoverage(obj.chip,bamFile.chip)
        obj.control=getReadCoverage(obj.control,bamFile.control)
        
        # normalize libSize difference
        libSizeDiff=libSize(obj.chip)/libSize(obj.control)
        
        cat("library size ratio (chip/control) is",libSizeDiff,"\n")
        
        if (libSizeDiff<1) {
            
            cat("Sample control data for normalizing library size difference...","\n")
            libSize.normalize=min(libSize(obj.chip),libSize(obj.control))
            readAlignment(obj.chip)=sample(readAlignment(obj.chip),libSize.normalize)
            readAlignment(obj.control)=sample(readAlignment(obj.control),libSize.normalize)
            
        }
        
        obj.chip=findPeaks(obj.chip,obj.control)
        
        
    })

## performance
## chr12 with H3 as control
# $sampling.time
# [1] 91.2

