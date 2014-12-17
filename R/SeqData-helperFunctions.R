## helperFunctions 
# 
###############################################################################

##-----------------------------------------------------------------------------
## .libSize
##

# param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
#                                     isUnmappedQuery=FALSE,
#                                     isNotPassingQualityControls=FALSE)

## count libSize for bamFile
## works for BamFile too

##' @export .libSize
.libSize=function(bamFile=character(0),param=ScanBamParam()){
    
    count=countBam(bamFile,param=param)
    libSize=count$records
    return(libSize)
}


## count libSize for BamFileList
##' @export .bams.libSize
.bams.libSize=function(bams=BamFileList(),param=ScanBamParam()){
    
    count=countBam(bams,param=param)
    libSize=count$records
    names(libSize)=as.character(count$file)
    return(libSize)
}

##----------------------------------------------------------------------------- 
## .rpkm and .rpkm.libSize
##

# simple rpkm normalization method

# normalized to total read count over regions of interest
# for comparison within one feature
##' @export .rpkm
.rpkm=function(counts,numBasesCovered){
    # geneAnnot is a list of exon ranges
    # geneAnnot=split(geneRanges,name)
    
    # millionsMapped=sum(counts)/1e6
    
    ## counted reads / total of counted reads in millions
    # rpm=counts/millionsMapped	
    
    # geneLengthInKB=numBasesCovered/1e3
    
    ## reads per million per geneLength in Kb
    # rpkm=rpm/geneLengthInKB
    
    ## counts=countData
    
    ## geneAnnot is a list of exon ranges
    ## numBasesCovered=sum(width(geneAnnot))
    
    rpkm=counts/(sum(counts)*1e-9*numBasesCovered)
    
}

# normalized to total mapped reads in this sample, libSize
# for any feature comparison in one sample.
# note rpkm is not a good parameter for cross sample comparison.

##' @export .rpkm.libSize
.rpkm.libSize=function(counts,libSize,numBasesCovered){
    rpkm=counts/(libSize*1e-9*numBasesCovered)
}


## normalize for count.table with libSize

##' @export .rpkm.table.libSize
.rpkm.table.libSize=function(count.table,libSize,numBasesCovered){
    
    name=colnames(count.table)
    rpkm.lst=list()
    
    for (i in 1:length(libSize)) {
        cat("processing column",i,"\n")
        rpkm.lst[[i]]=.rpkm.libSize(count.table[,i],libSize[i],numBasesCovered) 
    }
    rpkm.df=data.frame(do.call(cbind,rpkm.lst))
    names(rpkm.df)=name
    return(rpkm.df)	
}

##-----------------------------------------------------------------------------
## .checkChrNames
##

##' @export .checkChrNames
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
##

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
##

# add time stamp and bam file name as a unique signature of the output file
.timeStamp=function(filename){
   
    basename=basename(filename)
    name=unlist(strsplit(basename,split="[.]"))
    fileName=paste(name[1],"-",format(Sys.time(),"%Y%m%d.%H%M%S"),sep="")
    
}

##-----------------------------------------------------------------------------
## .samplingNorm
##

# normalize libSize/sequencing depth difference between libraries
# by sampling the bigger library to match the samll library
# then get the coverage

##' @export .samplingNorm
.samplingNorm=function(obj.chip,obj.control) {
    cat ("Normalizing library size difference...\n")
    
    libSizeDiff=libSize(obj.chip)/libSize(obj.control)
    cat("library size ratio (chip/control) is",libSizeDiff,"\n")
    
    # normalize only when libSize are different
    if (libSizeDiff!=1){
        cat("Sampling bigger library to match smaller library...","\n") 
        libSize.normalize=min(libSize(obj.chip),libSize(obj.control))
        readAlignment(obj.chip)=sample(readAlignment(obj.chip),libSize.normalize)
        readAlignment(obj.control)=sample(
            readAlignment(obj.control),libSize.normalize) 
        
        #change libSize slot of both library
        libSize(obj.chip)=length(readAlignment(obj.chip))
        libSize(obj.control)=length(readAlignment(obj.control))
        
        cat("Computing readCoverage based on normalized libSize...\n")
        obj.chip=getReadCoverage(obj.chip)
        obj.control=getReadCoverage(obj.control) 
        
    }else{
        # check whether readCoverage slot already has value 
        cat("libSize of obj.chip and obj.control are the same, compute readCoverage...\n")
        if (length(readCoverage(obj.chip))==0) obj.chip=getReadCoverage(obj.chip)
        if (length(readCoverage(obj.control))==0) obj.control=getReadCoverage(obj.control)
        
    }
    
    # return obj.chip and obj.control in a list, 
    # use chipSet$chip and chipSet$control subsetting each of them
    chipSet=list(chip=obj.chip,control=obj.control)
    return(chipSet)
    
}

##-----------------------------------------------------------------------------
## .seqlengths
##

##' @export .seqlengths
# output seqlengths based on whether one want "chr" to appear in the seqinfo
.seqlengths=function(chr=T,MT=T){
    seqlengths.mm9=c(
        197195432, 129993255, 121843856, 121257530, 120284312,
        125194864, 103494974, 98319150,  95272651,  90772031,
        61342430, 181748087, 159599783, 155630120, 152537259,
        149517037, 152524553, 131738871, 124076172, 16299, 
        166650296,  15902555)  
    mode(seqlengths.mm9)="integer"
    if (chr==T){
        if (MT==T){
            names(seqlengths.mm9)=c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15",
                "chr16","chr17", "chr18","chr19","chr2","chr3","chr4","chr5",
                "chr6","chr7","chr8", "chr9", "chrMT", "chrX", "chrY")
        }else{
            names(seqlengths.mm9)=c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15",
                "chr16","chr17", "chr18","chr19","chr2","chr3","chr4","chr5",
                "chr6","chr7","chr8", "chr9", "chrM", "chrX", "chrY")
        }
        
    }else{
        if (MT==T){
            names(seqlengths.mm9)=c(
                "1","10","11","12","13","14","15","16","17", "18","19",
                "2","3","4","5", "6","7","8", "9", "MT", "X", "Y")
        }else{
            names(seqlengths.mm9)=c(
                "1","10","11","12","13","14","15","16","17", "18","19",
                "2","3","4","5", "6","7","8", "9", "M", "X", "Y")
        }
    }  
    return(seqlengths.mm9)
}

##-----------------------------------------------------------------------------
## .convertChrNames
##


## .seqlengths is output a value for the supposed length of each chromosome
## for chr, you may only want to convert it to a proper form you want, instead of
## reset it. This function let one choose the output 
## have "chr" in or not, have use "MT" or "MT"

##' @export .convertChrNames
.convertChrNames=function(chrNames,chr=T,MT=T){
    
    ## four helper functions
    ## remove "chr"
    .rmChr=function(chrNames){substr(chrNames,"4","6")}
    ## add "chr"
    .addChr=function(chrNames){paste("chr",chrNames,sep="")}
    
    ## substitute "MT" to "M"
    .subMT=function(chrNames){sub("MT","M",chrNames)}
    ## substitute "M" to "MT"
    .subM=function(chrNames){sub("M","MT",chrNames)}
    
    
    # convert names to character
    if(is.factor(chrNames)) chrNames=as.character(chrNames)
    
    # check chr MT status
    cat("input:\n")
    grepChr= grep("chr",chrNames)
    cat("chr",length(grepChr),"\n")
    grepMT=grep("MT",chrNames)
    cat("MT",length(grepMT),"\n\n")
    
    cat("output:\n")
    chr.out=ifelse(chr,1,0)
    cat("chr.out =",chr.out,"\n")
    
    MT.out=ifelse(MT,1,0)
    cat("MT.out =",MT.out,"\n\n")
    
    # convert
    # swtich states what you want
    # if statement checks what's there in the names
    if (chr==T)
    {# chr.out=1                        
        if (MT==T)
        { # chr.out=1, MT.out=1
            if (length(grepChr)==0){
                # chr=0
                if (length(grepMT)==0){
                    # chr=0, MT=0
                    chrNames=.addChr(.subM(chrNames))
                }else{
                    # chr=0, MT=1
                    chrNames=.addChr(chrNames)
                }
            }else{
                # chr=1
                if (length(grepMT)==0){
                    # chr=1, MT=0
                    chrNames=.subM(chrNames)
                }else{
                    # chr=1, MT=1
                    chrNames
                }     
            }
            
        }else
        {# chr.out=1, MT.out=0
            
            cat("MT.out =",MT.out,"\n\n")
            
            if (length(grepChr)==0){
                # chr=0
                if (length(grepMT)==0){
                    # chr=0, MT=0
                    chrNames=.addChr(chrNames)
                }else{
                    # chr=0, MT=1
                    chrNames=addChr(.subMT(chrNames))
                }
            }else{
                # chr=1
                if (length(grepMT)==0){
                    # chr=1, MT=0
                    chrNames
                }else{
                    # chr=1, MT=1
                    chrNames=.subMT(chrNames)
                }     
            } 
        }
        
        
    }else
    { # chr.out=0
        if (MT==T)
        { # chr.out=0, MT.out=1
            if (length(grepChr)==0){
                # chr=0
                if (length(grepMT)==0){
                    # chr=0, MT=0
                    chrNames=.subM(chrNames)
                }else{
                    # chr=0, MT=1
                    chrNames
                }
            }else{
                # chr=1
                if (length(grepMT)==0){
                    # chr=1, MT=0
                    chrNames=.rmChr(.subM(chrNames))
                }else{
                    # chr=1, MT=1
                    chrNames=.rmChr(chrNames)
                }     
            }
        }else
        { # chr.out=0, MT.out=0
            if (length(grepChr)==0){
                # chr=0
                if (length(grepMT)==0){
                    # chr=0, MT=0
                    chrNames
                }else{
                    # chr=0, MT=1
                    chrNames=.subMT(chrNames)
                }
            }else{
                # chr=1
                if (length(grepMT)==0){
                    # chr=1, MT=0
                    chrNames=.rmChr(chrNames)
                    
                }else{
                    # chr=1, MT=1
                    chrNames=.rmChr(.subMT(chrNames))
                }     
            }
            
        }
        
        
    }
    
    
    return(chrNames)
    
}

## Notes on switch
## only use it as the flow control
## when use it as return, it is not as straight as if-else

##-----------------------------------------------------------------------------
## df2gr
##

## convert data.frame to GRanges based on whether one wants "chr" in the 
## chromosome name and whether one want "MT" or "M" in the chromosome name
## require seqlengths,.convertChrNames to work

##' @export df2gr
df2gr=function(df,chr=T,MT=T){
    
    
    # add strand if there isn't 
    cln=colnames(df)
    if (length(which(cln=="strand"))!=0) strand=df$strand 
    else strand=rep("*",dim(df)[1]) 
    
    # add meta columns        
    conserved.n=c("chr","start","end","strand")
    mcols.n=setdiff(cln,conserved.n)
    mcols.n=mcols.n[complete.cases(mcols.n)]
    
    ## remember here to use [[]] instead of $x or [,x]
    mcols=sapply(mcols.n,function(x){ df[[x]] })    
    
    ## chr MT
    df$chr=.convertChrNames(df$chr,chr=chr,MT=MT)
    gr=GRanges(
        seqnames=df$chr,
        ranges=IRanges(start=df$start,end=df$end),
        strand=strand,
        #mcols=mcols,
        mcols,  # this removes mcols.xxx in metacolum names
        seqlengths=.seqlengths(chr=chr,MT=MT)
    )
    
    # use suppressMessages() if don't want to see warning
    gr=trim(gr)
    cat("if warning about trimming GRanges appears, it is already trimmed\n")
    return(gr)
    
}


# the names of seqlengths(aln) is lost when called using lapply function, it is
# passed in as a second variable for lapply function
# bin.size in bp

##' @export .bins
.bins=function(bin.size,granges=T){
    SEQLEN=.seqlengths()
    bins <- GRangesList(
        lapply(SEQLEN,function(seqlen,seqlen.name) {
            GRanges(
                ranges=IRanges(breakInChunks(seqlen,bin.size)), 
                seqnames=Rle(seqlen.name[which(SEQLEN==seqlen)],length(ranges)),    
                strand=Rle("*",length(ranges)),
                seqlengths=SEQLEN[which(SEQLEN==seqlen)]
            )},names(SEQLEN))
    )
    if(granges) bins=unlist(bins,use.names=FALSE)
    return(bins)
}



