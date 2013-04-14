# countReadsOverGenes-Methods
# countReadsOverFeatures-Methods

# FIXME: use countReadsOverRanges  seems more fit and goes better with other title.

# 

##'@title countReadsOverFeatures
##' @name countReadsOverFeatures
##' @description countReadsOverFeatures
##' count reads over different genomic features, such as gene, transciripts, exons, promoter regions, or any genomic ranges of interest. It will count how many reads is mapped to each genomic feature. 
##'  we use 1 kb upstream and downstream of tss represent promoter regions for now, feature name tss. 

##@rdname countReadsOverFeatures-method
## @details this function... \itemize{\item {a} {does a} \item {annotate} {A switch for stop the pipeline at annotatePeaks()}} 
##' @details \describe{\item{features :}{ description a}}
###############################################################################

setMethod(
    f="getReadAlignment",
    signature="SeqData",
    definition=function(obj,bamFileName){
        
        
        ## TODO: add function getBamInfo 
        # play with the bam file can get more information about mappability
        # mapped vs unmapped, multimatch vs unique match
        
        # param=ScanBamParam(what=c("isize"),
        # flag=scanBamFlag(isUnmappedQuery=False,isPaired=NA),
        # tag="NH")
        
        # every function seems need the bam file 
        # set libSize
        # set readCoverage
        # set readAlignment  for countOverlaps and coverage
        
        # and other informations that can potentially come from the bam file/the alignment file.
        
        
        # TODO: get species which the bam file is aligned to
        
         
         # TODO: set a detector, if readCoverage(obj) is not empty, length()==0, then do this
     
        
#         if (pairEnd){
#             param=ScanBamParam(what="isize")
#             
#             # something to consider later: mappability
#             # param=ScanBamParam(what=c("isize"),
#             # flag=scanBamFlag(isUnmappedQuery=False,isPaired=NA),
#             # tag="NH")
#             
#             
#             
#             
#             
#         }
        
        
       # ifelse (pairEnd==TRUE,param=ScanBamParam(what="isize"),param=NULL)
       
       ### you don't need to know whether it's single end or pairEnd, just look at isize, as it can be applied to both single end and pair end
        
        
        

        # read in bam file
        # make bam index file if it doesn't exist
        if(!file.exists(paste(bamFileName,"bai",sep="."))){        		
            cat("Making index file for bam...","\n")
            indexBam(bamFileName)
        } 
        
        ### extract information from bamHeaders	 
        bamHeader= scanBamHeader(bamFileName) 
        
        # get chromSize
        # chrSize(obj)=bamHeader[[1]]$targets # don't need this slot, as no other function needs it
        chrSize=bamHeader[[1]]$targets
        names(chrSize)=.checkChrNames(names(chrSize))
        chrSize(obj)=chrSize
        
        
        # TODO: get species the bam file is aligned to
        
        ### read in bam file
        cat("Processing bam file...","\n")
        #readAlignment(obj)=readGappedAlignments(bamFileName,index=bamFileName,format="BAM",use.names=T)
        
        param=ScanBamParam(what="isize")
        
        aln=readGappedAlignments(bamFileName,index=bamFileName,format="BAM",use.names=T, param=param)
        
        # TODO: add a checker here only do it when seqlevels doesn't have "chr"
        #seqlevels(aln)=paste("chr",seqlevels(aln),sep="")
        
        
        seqlevels(aln)=.checkChrNames(seqlevels(aln))
        
        cat("Setting readAlignment slot...","\n")
        
        readAlignment(obj)=aln
        
        
        
        # set libSize slot
        libSize(obj)=length(aln)  # don't need this slot as it is a simple function of readAlignment(obj)
        
        


        return(obj)
    })

# FIXME Redundency of reading bam


setMethod(
    f="getFeatureAnnotation",
    signature="SeqData",
    definition=function(obj,annotationFileName,feature=c("gene","tss","exon","transcript")){
        
        cat("Processing annotation file...\n")
        refGene=read.csv(annotationFileName,as.is=T,header=T)
        
        # duplication check
        # sometimes exactly same item appears in UCSC refGene table, 
        # which should be removed before processing
        
        chr.start=paste(refGene$chrom,refGene$txStart,sep="_")    
        uniGene=paste(chr.start,refGene$name)
        refGene$names=uniGene
        
        refGene=refGene[!duplicated(refGene$names),]
        
        featureAnnot=switch(feature,
                            "gene"={featureAnnotation(obj)=.geneAnnot(refGene)},
                            "tss"={featureAnnotation(obj)=.tssAnnot(refGene)},
                             {stop(cat("No match was found for feature requested,please check","\n"))}
                            )

        return(obj)   
    }) 

##
##




setMethod(
    f="countReadsOverFeatures",
    signature="SeqData",
    definition=function(obj=NULL,bamFileName=character(1),annotationFileName=character(1),countingMode="Union"){
        
        ## check featureAnnotation
        if (length(featureAnnotation(obj))==0){
            stop("featureAnnotation hasn't been assigned yet")
            #obj=getFeatureAnnotation(obj,featureAnnotation)
        }
        
        ##check readAlignment
        if(length(readAlignment(obj))==0){ stop("readAlignment has not been assigned yet")
        }
        
        
        # matching names in the readAlignment and featureAnnotation
        # not only the "chr" issue
        # but also when they have different number of chromosome
        
        # check chromosome naming convention of the bam file
        seqlevels(readAlignment(obj))=.checkChrNames(seqlevels(readAlignment(obj)))
        
        # matching the chromosome number in aln and anno        
        intersectChrNames <- sort(intersect(seqlevels(readAlignment(obj)),seqlevels(featureAnnotation(obj))))
        
        ## subsetting GRanges multiple chromosome wise (a method for looping over GRanges)
        # subsetting GRanges
        intersectAlignment= do.call("c",lapply(intersectChrNames,function(x) readAlignment(obj)[seqnames(readAlignment(obj))==x]))
        
        # adjusting corespond seqlevels(It is important to do this step)                            
        seqlevels(intersectAlignment)=intersectChrNames
        
        # same for annotation       
        intersectGeneAnnot= do.call("c",lapply(intersectChrNames,function(x) featureAnnotation(obj)[seqnames(featureAnnotation(obj))==x]))                                            
        seqlevels(intersectGeneAnnot)=intersectChrNames
        
        # TODO: a switch for user to select different features: gene(default), exon, transcripts, promoter, TSS

        geneAnnot=split(intersectGeneAnnot,mcols(intersectGeneAnnot)$gene) 
        
        cat("Counting reads per gene...\n")
        
        featureHits=summarizeOverlaps(geneAnnot,intersectAlignment,mode=countingMode)
        
        countData=assays(featureHits)$counts
        
        return(countData)
        
    }
)


setMethod(
    f="countReadsOverFeatures",
    signature="missing",
    definition=function(obj=NULL,bamFileName=character(1),annotationFileName=character(1),countingMode="Union"){
        
        # TODO: set an error handler here to handle wrong input
        
        #if (length(bamFileName)==0||length(annotationFileName==0)){
        if (bamFileName==character(1)||annotationFileName==character(1)){   
            bamFileName=readline("\nPlease input the path to your bam file 
                                 or simply drag the file to this window: ")
            annotationFileName=readline("\nPlease input the path to your annotation file
                                        or simply drag the file to this window: ")   
        }
        
        
        
        obj=new("SeqData")
        
        obj=getReadAlignment(obj,bamFileName)
        
        obj=getFeatureAnnotation(obj,annotationFileName)
        
        show(obj)
        
        print(countingMode)
        
        countData=countReadsOverFeatures(obj,countingMode=countingMode)
        
        cat("output countData...\n")
        write.csv(file="countData.csv",countData[,1])
        
        cat("\nDone!\n")
        
        return(obj)
        })
