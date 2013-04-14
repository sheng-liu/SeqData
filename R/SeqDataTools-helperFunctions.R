# internalMethods
# 
# 
###############################################################################

## check and change chromosome names into standard UCSC naming convention "chr12"
.checkChrNames=function(chrNames){
    
    # convert names to character
    if(is.factor(chrNames)) chrNames=as.character(chrNames)
    
    
    # check and change chr to UCSC chr naming convention
    grepChr= grep("chr",chrNames)
    if (length(grepChr) == 0) chrNames=paste("chr",sub("MT","M",chrNames),sep="")
    # go back: chrNames=substr(chrNames,"4","6")
    
    return(chrNames)
}

## strand correlation method to determine meanFragmentSize
## calculate the coverage for both positive and negtive strand then determine the best shift that maximizes the correlation between these two coverages, the bestShift is then the estimated meanFragmentSize. 

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

.geneAnnot=function(refGene){
    # indexing exons for each gene
    # TODO: change the name convention here
    
    exon.index=c()
    for (i in 1:length(refGene$chrom)) {    		
        if (refGene$strand[i]=="+")		   		
            s=seq(from=1,to=refGene$exonCount[i])    				
        else
            s=seq(from=refGene$exonCount[i],to=1)					
        exon.index=c(exon.index,s)
    }
    
    
    
    # generate RangedData for annotation
    geneAnnot=RangedData(IRanges(
        start=as.integer(unlist(strsplit(refGene$exonStarts,split=","))),
        end=as.integer(unlist(strsplit(refGene$exonEnds,split=",")))),
                         space=rep(refGene$chrom,refGene$exonCount),
                         strand=rep(refGene$strand,refGene$exonCount),
                         exon=paste(rep(refGene$name,refGene$exonCount),exon.index,sep=":"),
                         transcript=rep(refGene$name,refGene$exonCount),
                         gene=rep(refGene$name2,refGene$exonCount))
    
    # corce RangedData to GRanges object
    cat("Setting featureAnnotation slot...","\n")
    geneAnnot=as(geneAnnot,"GRanges")
    return(geneAnnot)
}

.tssAnnot=function(refGene){
    # create tssAnnot
    # using GRanges
    tssAnnot=GRanges(
        seqnames=refGene$chrom,
        ranges=IRanges(
            start=refGene$txStart,
            end=refGene$txStart),
        strand=refGene$strand,
        EntrezID=refGene$name,
        symbol=refGene$name2)          
    #tssAnnot=resize(tssAnnot,width=2000,fix="center")
    
    
    #             # using IRanges
    #             tssRanges=IRanges(start=refGene$txStart,end=refGene$txStart)
    #             tssRanges=resize(tssRanges,width=2000,fix="center")
    #             tssAnnot=RangedData(tssRanges,
    #                                 space=refGene$chrom,
    #                                 strand=refGene$strand,
    #                                 EnrezID=refGene$name,
    #                                 symbol=refGene$name2)
    #             tssAnnot=as(tssAnnot,"GRanges") 
    return(tssAnnot)
}



