# getReadCoverage
# 
# 
###############################################################################


## many of the dispatch here can be removed, user would only use either SeqData with bamFile filled, or a bamFile itself. nobody is going to use an empty SeqData,BamFile, or to say, the author shouldn't encourage this usage. 

# usage

# getReadCoverage(bamFile)
# getReadCoverage(obj,bamFile)
# getReadCoverage(obj)             # obj needs to have bamFile slot non-empty

# return
# SimpleRleList composed of smoothed read coverage

setMethod(
    f="getReadCoverage",
    signature=c(obj="SeqData",bamFile="character"),
    definition=function(obj,bamFile=character(1)        
    ){
        
        # set bamFile slot       
        bamFile(obj)=bamFile
        
        # check if readAlignment slot is empty, if it is then read in readAlignment first
        if(length(readAlignment(obj))==0){
            obj=getReadAlignment(obj)
        }else{
            cat("readAlignment is not empty, use existing readAlignment\n")
        }
        
        # get isize
        isize=mcols(readAlignment(obj))
        
        # estimate meanfragmentSize for pair end data
        meanFragmentSize=round(mean(abs(isize[,1]),na.rm=T))

        
        # estimate meanfragmentSize for single end data
        if (meanFragmentSize==0){
            cat("SeqData is likely to be single-end, estimating meanFragmentSize...")
            meanFragmentSize=round(.estimateFragmentSize(readAlignment(obj)))
                                                                                                 
            cat("estimated meanFragmentSize is ",meanFragmentSize,"\n")
        }else{
            cat("SeqData is pair-end, meanFragmentSize is ",
                meanFragmentSize,"\n")
        }
        
        # coerce GAlignments to GRanges
        aln <- as(readAlignment(obj), "GRanges")
        
        
        # extend the reads
        if (meanFragmentSize==0){
            cat("Extendeding fragments...\n")
            aln.extended= suppressWarnings(resize(aln, width=meanFragmentSize)) 
        }else{
            cat("Fragments are not extended\n")
            aln.extended=aln
            
        }
       
        # compute bp coverage 
        cat("Computing read coverage...","\n")        
        readCoverage(obj)=coverage(aln.extended,width=seqlengths(readAlignment(obj)))
        
        return(obj)
    })

# a dispatcher when input is SeqData
setMethod(
    f="getReadCoverage",
    signature=c(obj="SeqData",bamFile="missing"),
    definition=function(obj){
        
        # check whether bamFile slot is empty
        if (length(bamFile(obj))==0) {
            stop ("Slot bamFile is empty, please fill in bamFile slot first","\n")    
        }
        
        # get coverage
        obj=getReadCoverage(obj,bamFile(obj))
        
        return(obj)      
    })

# a dispatcher when input is bamFile
setMethod(
    f="getReadCoverage",
    signature=c(obj="missing",bamFile="character"),
    definition=function(bamFile){
        
        obj=new("SeqData")        
        bamFile(obj)=bamFile       
        obj=getReadCoverage(obj,bamFile=bamFile)         
        return(obj)
        
    })

##-----------------------------------------------------------------------------
## helperFunctions


# estimate meanFragmentSize based on cross-correlation between positive and negtive strand, the bestShift is the shift that gives the highest cross-correlation. 

.estimateFragmentSize=function(aln){
    
    # find first unempty chromosome to do the estimation    
    chr=as.character(seqnames(aln)[1])
    aln.sample= aln[seqnames(aln)==as.integer(chr)]
    
    # if want to sample the data to reduce calculation
    # sampleSize=2*10e9   # the range of integer is -2*10e9~ 2*10e9
    # aln.sample=sample(aln,sampleSize)
    
    # calculate coverage for both positive and negative strand
    posStr=strand(aln.sample)=="+"
    posCover=coverage(aln.sample[posStr])
    negCover=coverage(aln.sample[!posStr])
    
    # compute the bestShift
    shifts=seq(100,500,by=3)        
    corrs=shiftApply(shifts,negCover[[chr]],posCover[[chr]],FUN=cor)
    
    shiftsTable=data.frame(shifts,corrs)
    bestShift=shiftsTable$shifts[which.max(shiftsTable$corrs)]
    
    return(bestShift)
    
}

# for chr12
# $sampling.time
# [1] 98.44
