# getBaseScore-method
# 
# 
###############################################################################

## get inidvidual base coverage
## can export bedGraph or bigWig directly from the baseAlignment
## this funciton unify it to coverage, so it is easy to output either genome wide and small scale views

# result output to coverage slot of SeqData

# provide a formula function for user to define their way of calculating score

# bases, basesAlignment 
# FUN, function used to manipulate metacolumns of the bases (currently support two metacolumns)

# for coverage
##' @export .cov
.cov=function(x,y){x+y}

##' @export .perc
# for percent of methylation
.perc=function(x,y){  
    perc=100*x/(x+y)
    perc[is.na(perc)]=0
    perc=as.integer(round(perc)) 
}
# to facilitate viewSums, make it integer, otherwise it is numeric, viewSums only performs on integer coverage
#perc=round(perc,digits=2) # round number to 2 digits
# perc[is.na(perc)]=-1 # replace NA 
# dispatch on SeqData with baseAlignment filled
# and baseReport as well
# obj with baseAlignment filled
# default FUN=.perc
# default output="coverage"

setMethod(
    f="getBaseScore",
    signature=c(obj="SeqData"),
    definition=function(obj,bamFile=character(0),baseReport=character(0),FUN=.perc,output="coverage"){ 
        
        #check if the baseAlignment slot is filled
        if(length(baseAlignment(obj))==0) 
            stop("Please fill in baseAlignment slot first.")
        
        baseScore=.getBaseScore(bases=baseAlignment(obj),FUN=FUN,output=output)
        
        if(output=="coverage"){
            readCoverage(obj)=baseScore
        }else{
            return(baseScore)
        }
        
        return(obj)
        
    })   
   
# dispatch on baseReport
setMethod(
    f="getBaseScore",
    signature=c(baseReport="character"),
    definition=function(obj=NULL,bamFile=character(0),baseReport=character(0),FUN=.perc,output="coverage"){ 
        
        obj=getBaseAlignment(baseReport=baseReport)
        obj=getBaseScore(obj,FUN=FUN)
        return(obj)
        
        })   
        
.getBaseScore=function(bases,FUN,output=c("GRanges","coverage")){
    
    FUN=match.fun(FUN)
    
    print(FUN)
    
    baseInfo1=mcols(bases)[[1]]
    baseInfo2=mcols(bases)[[2]]
    
    score=FUN(baseInfo1,baseInfo2)
    DF=DataFrame(score)
    
    baseScore=bases
    mcols(baseScore)=DF
    # print(DF)
    # change the score into coverage format, RleList, 
    # for integration with viewCoverage method
    
    if(output=="coverage"){
        baseScore.chr=split(baseScore,seqnames(baseScore))
        # this gives genomewide 1 base coverage
        baseScore.cov=coverage(baseScore) 
        for (chr in seqlevels(bases)){   
            print(chr)
            # this converts GRanges into coverage (RleList)
            baseScore.cov[[chr]][start(baseScore.chr[[chr]])]=mcols(baseScore.chr[[chr]])[[1]] 
           
        }
            return(baseScore.cov)  
    
    }
    
    return(baseScore)
    
}

## output is integer-Rle List on all list members



