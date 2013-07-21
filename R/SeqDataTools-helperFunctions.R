## helperFunctions 
# 
###############################################################################

##-----------------------------------------------------------------------------
## .libSize

.libSize=function(bamFile=character(0),
                  param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                                      isUnmappedQuery=FALSE))){
    
    count=countBam(bamFile,param=param)
    libSize=count$records
    return(libSize)
}


##----------------------------------------------------------------------------- 
## .rpkm and .rpkm.libSize

# simple rpkm normalization method

# normalized to total read count over regions of interest
# for comparison within one feature
.rpkm=function(counts,numBasesCovered){
    # millionsMapped=sum(counts)/1e6
    
    ## counted reads / total of counted reads in millions
    # rpm=counts/millionsMapped	
    
    # geneLengthInKB=numBasesCovered/1e3
    
    ## reads per million per geneLength in Kb
    # rpkm=rpm/geneLengthInKB
    
    ## counts=countData
    ## numBasesCovered=sum(width(geneAnnot))
    
    rpkm=counts/(sum(counts)*1e-9*numBasesCovered)
    
}

# normalized to total mapped reads in this sample, libSize
# for any feature comparison in one sample.
# note rpkm is not fit for cross sample comparison. 
.rpkm.libSize=function(counts,libSize,numBasesCovered){
    rpkm=counts/(libSize*1e-9*numBasesCovered)
}


##-----------------------------------------------------------------------------
## .checkChrNames

## check and change chromosome names into standard UCSC naming convention "chrX"
.checkChrNames=function(chrNames){
    
    # convert names to character
    if(is.factor(chrNames)) chrNames=as.character(chrNames)
    
    # check and change chr to UCSC chr naming convention
    grepChr= grep("chr",chrNames)
    if (length(grepChr) == 0) chrNames=paste("chr",sub("MT","M",chrNames),sep="")
    # if want to go back: chrNames=substr(chrNames,"4","6")
    
    return(chrNames)
}

.checkChrNames.rev=function(chrNames){
    
    # convert names to character
    if(is.factor(chrNames)) chrNames=as.character(chrNames)
    
    # check and change chr to to match naming convention
    grepChr= grep("chr",chrNames)
    if (length(grepChr) != 0) chrNames=chrNames=substr(chrNames,"4","6")
    return(chrNames)
}
    
    
    


##-----------------------------------------------------------------------------
## intersectChr

# a method for intersecting chr name and number between reads (readCoverage, readAlignment) with features (featureAnnotation)

setMethod(
    f="intersectChr",
    signature=c("GAlignments","GRanges"),
    definition=function(reads,features){
      
        # check chromosome naming convention of the reads and features
        seqlevels(reads)=.checkChrNames(seqlevels(reads))
        seqlevels(features)=.checkChrNames(seqlevels(features))
        
        # matching the chromosome number in reads and features  
        intersectChrNames=sort(intersect(seqlevels(reads),seqlevels(features)))
        
        ## subsetting GRanges multiple chromosome wise (a method for looping over GRanges)
        # subsetting GRanges
        intersectReads=do.call("c",lapply(intersectChrNames,function(x) reads[seqnames(reads)==x]))
        
        # adjusting corespond seqlevels(It is important to do this step)  
        seqlevels(intersectReads)=intersectChrNames
        
        # same for annotation
        intersectFeatures=do.call("c",lapply(intersectChrNames,function(x) features[seqnames(features)==x]))
        
        seqlevels(intersectFeatures)=intersectChrNames
        seqlengths(intersectFeatures)=seqlengths(intersectReads)
        
        # return reads and features together in a list, use intersect$reads and intersect$features subsetting each of them
        intersect=list(reads=intersectReads,features=intersectFeatures)
        
        return(intersect)
        
        
    })

setMethod(
        f="intersectChr",
        signature=c("SimpleRleList","GRanges"),
        definition=function(reads,features){
            
            # check chromosome naming convention 
            names(reads)=.checkChrNames(names(reads))
            seqlevels(features)=.checkChrNames(seqlevels(features))
                        
            intersectChrNames=sort(intersect(names(reads),seqlevels(features)))
            
            # subsetting SimpleRleList
            intersectReads=do.call("c",lapply(intersectChrNames,function(x) reads[names(reads)==x]))
            
            # adjusting corespond seqlevels(It is important to do this step)  
            names(intersectReads)=intersectChrNames
            
            # same for annotation
            intersectFeatures=do.call("c",lapply(intersectChrNames,function(x) features[seqnames(features)==x]))
            
            seqlevels(intersectFeatures)=intersectChrNames
            
            #set seqlengths(intersectFeatures)
            chrl=intersectChrNames
            chrLength=c()
            chrLength=sapply(chrl,function(chr) {
                chrLen=length(intersectReads[[chr]])
                chrLength=c(chrLen,chrLength)
                
            })
            seqlengths(intersectFeatures)=chrLength

            # return reads and features together in a list, use intersect$reads and intersect$features subsetting each of them
            intersect=list(reads=intersectReads,features=intersectFeatures)

            return(intersect)
            
            
        })

##-----------------------------------------------------------------------------
## .timeStamp

# add time stamp and bam file name as a unique signature of the output file
.timeStamp=function(filename){
   
    basename=basename(filename)
    name=unlist(strsplit(basename,split="[.]"))
    fileName=paste(name[1],"-",format(Sys.time(),"%Y%m%d.%H%M%S"),".csv",sep="")
    
}



