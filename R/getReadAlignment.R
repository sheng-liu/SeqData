# getReadAlignment
#
#
###############################################################################


## this is a method to set readAlignment slot of SeqData
# -create a new SeqData objects and return the filled one, if no SeqData object is passed in
# -fill the readAlignment slot of SeqData object, if SeqData object is passed in, require bamFile slot non-empty

# function can determine whether SeqData is single end or pair end by looking at isize when read in the bamFile

# when and only when annotationFile are passed in, getReadAlignment will also set featureAnnotation Slot for record, the final output is the readAlignment at those ranges specified at feature annotation slot.

# usage
# getReadAlignment(bamFile=bamFile)
# getReadAlignment(obj,bamFile=bamFile)
# getReadAlignment(obj)         # obj needs to have bamFile slot non-empty
# getReadAlignment(bamFile=bamFile,annotationFile=annotationFile)

setMethod(
    f="getReadAlignment",
    signature=c(obj="SeqData",bamFile="character"),
    definition=function(obj,bamFile){
        
        # set bamFile slot
        bamFile(obj)=bamFile
        
        # make bam index file if it doesn't exist
        if(!file.exists(paste(bamFile,"bai",sep="."))){                
            cat("Making index file for bam...","\n")
            indexBam(bamFile)
        } 
        
        # extract information from bamHeaders	 
        bamHeader= scanBamHeader(bamFile) 
        
        # get chrSize
        chrSize=bamHeader[[1]]$targets
        names(chrSize)=.checkChrNames(names(chrSize))
        chrSize(obj)=chrSize
        
        # read in bam file
        cat("Processing bam file...","\n")
        
        # read in isize to determine its insert size also determine whether it is pair end or single end
        param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isUnmappedQuery=FALSE),what=c("isize"))
        
        # note readGAlignments/readGAlignmentsFromBam default param is
        # flag=scanBamFlag(isUnmappedQuery=FALSE)
        # explicitly stated here for calrification
        aln=readGAlignments(bamFile,index=bamFile,format="BAM",use.names=T, param=param)
        
        # check chrNames
        seqlevels(aln)=.checkChrNames(seqlevels(aln))
        
        cat("Setting up readAlignment slot...","\n")
        readAlignment(obj)=aln
        
        # set libSize slot
        libSize(obj)=length(aln)
        
        return(obj) 
        
    })

# a dispatcher when input is SeqData
setMethod(
    f="getReadAlignment",
    signature=c(obj="SeqData",bamFile="missing"),
    definition=function(obj){
        
        # check whether bamFile slot is non-empty
        bamFile=bamFile(obj) 
        if (length(bamFile)==0) {
            stop ("Slot bamFile is empty, please fill in bamFile slot first","\n")
        }
        
        obj=getReadAlignment(obj,bamFile)
        return(obj)
    })



# a dispatcher when input is bamFile
setMethod(
    f="getReadAlignment",
    signature=c(obj="missing",bamFile="character"),   
    definition=function(bamFile=character(0)){
        
        obj=new("SeqData")
        bamFile(obj)=bamFile
        
        
        obj=getReadAlignment(obj)
        
        return(obj)
        
        
    })




# get read alignment at ranges specified in annotation file 

setMethod(
    f="getReadAlignment",
    signature=c(obj="missing",bamFile="character",annotationFile="character"),
    definition=function(bamFile=character(0),annotationFile=character(0)){
        
        obj=new("SeqData")
        
        # set bamFile slot
        bamFile(obj)=bamFile
        
        # make bam index file if it doesn't exist
        if(!file.exists(paste(bamFile,"bai",sep="."))){                
            cat("Making index file for bam...","\n")
            indexBam(bamFile)
        } 
        
        # extract information from bamHeaders     
        bamHeader= scanBamHeader(bamFile) 
        
        # get chrSize
        chrSize=bamHeader[[1]]$targets
        names(chrSize)=.checkChrNames(names(chrSize))
        chrSize(obj)=chrSize
        
        
        # set featureAnnotation slot and get FeatureAnnotation
        obj=getFeatureAnnotation(obj,annotationFile=annotationFile)
        
        

        # set readAlignment slto
        readAlignment(obj)=.BamViewsAlignment(bamPaths=bamFile,bamRanges=featureAnnotation(obj))
        
        # set libSize slot
        cat("counting libSize...\n")
        libSize(obj)=.libSize(bamFile=bamFile)
        cat("libSize=",libSize(obj),"\n")
        
        return(obj)

    })


##-----------------------------------------------------------------------------
## helperFunctions

.BamViewsAlignment=function(bamPaths,bamRanges){
    #viewBam implementation
    # deal with one bam file at a time (for now and for simplicity and robusty)
    
    # check chrNames in bamFile, make chrNames in featureAnnotation match to it
    
    bamHeader= scanBamHeader(bamPaths)
    chrSize=bamHeader[[1]]$targets
    
    chrNames.bam=names(chrSize)
    chrNames.ranges=seqlevels(bamRanges)
    
    # convert chromosome names in bamRanges to match chromosome names in bam file when chromosome names in bamFiles doesn't have "chr"
    
    grepChr= grep("chr",chrNames.bam)
    
    if (length(grepChr) == 0) {
        chrNames.ranges=.checkChrNames.rev(chrNames.ranges)
        seqlevels(bamRanges)=chrNames.ranges
    }
    
    bv=BamViews(bamPaths=bamPaths,bamRanges=bamRanges)
    
    
    param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),what=c("isize"))
    aln.bv=suppressWarnings(readGAlignmentsFromBam(bv,param=param)[[1]])
    
    # check the chrosomsome name to UCSC naming convention 
    
    seqlevels(aln.bv)=.checkChrNames(seqlevels(aln.bv))
    
    
    return(aln.bv)
}


##-----------------------------------------------------------------------------
## 
## TODO: 

# add function getBamInfo 
# play with the bam file can get more information about
# mapped vs unmapped, multimatch vs unique match

# for this part, it is a function like flagstat in samtools, which is already implemented in Rsamtools as quickCountBam()

# later may add it. now the default param is 
# param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isUnmappedQuery=FALSE))
# accross the package.



# param=ScanBamParam(what=c("isize"),
# flag=scanBamFlag(isUnmappedQuery=False,isPaired=NA),
# tag="NH")





# get species which the bam file is aligned to

