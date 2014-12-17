# getFeatureAnnotation
# 
# 
###############################################################################

# get feature annotation by features from two inbuild annotation file refGene.csv and rmasker_repClassLTR.csv. features include "gene","exon","transcript","tss","erv","5LTR"

# can also get feature annotation from  csv format annotation file downloaded from UCSC if it has the miminum five field: 'chrom','txSart','txEnd','strand' and 'name'

# feature=c("gene","exon","transcript","tss","erv","repeats","transposon","5LTR")
# defaults feature is "gene".

# usage
# getFeatureAnnotation(annotationFile=annotationFile)
# getFeatureAnnotation(obj,annotationFile=annotationFile)
# getFeatureAnnotation(obj,feature="gene")

setMethod(
    f="getFeatureAnnotation",
    signature="SeqData",
    definition=function(obj,annotationFile=character(0),feature="gene"){
        
        
        # check whether annotationFile is present
        # if yes, use user supplied file 
        # if no, use user supplied feature
        
        if (length(annotationFile!=0)){
            cat("Please make sure the annotation file has required field:'chrom','txSart','txEnd','strand','name'.\n\n")
            
            cat("Setting up annotationFile slot...\n")
            annotationFile(obj)=annotationFile
            featureAnnot=read.csv(file=annotationFile(obj),as.is=T,header=T) 
            featureAnnot=.duplicationCheck(featureAnnot)
            
            cat("Setting up featureAnnotation slot...","\n")
            featureAnnotation(obj)=.featureAnnot(featureAnnot)
        }else{
            cat("Processing annotation file...\n")
            featureAnnot=
                switch(feature,
                       "gene"={
                           annotationFile(obj)=system.file(
                               "extdata", "refGene.csv", package="SeqData")
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           featureAnnotation(obj)=.geneAnnot(featureAnnot)                                
                       },
                       "erv"={
                           annotationFile(obj)=system.file(
                               "extdata", "rmasker_repClassLTR.csv", 
                               package="SeqData")
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           featureAnnotation(obj)=
                               .featureAnnot(featureAnnot)
                       },
                       "repeats"={
                           cat("Setting up featureAnnotation slot...","\n")
                           data(repeats)
                           featureAnnotation(obj)=repeats
                       },
                       "transposon"={
                           cat("Setting up featureAnnotation slot...","\n")
                           data(transposon)
                           featureAnnotation(obj)=transposon
                       },
                       "exon"={
                           annotationFile(obj)=system.file(
                               "extdata", "refGene.csv", package="SeqData")
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           featureAnnotation(obj)=.exonAnnot(featureAnnot)                                
                       },
                       "transcript"={
                           annotationFile(obj)=system.file(
                               "extdata", "refGene.csv", package="SeqData")
                           
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           featureAnnotation(obj)=.exonAnnot(featureAnnot) 
                           
                       },
                       "tss"={
                           annotationFile(obj)=system.file(
                               "extdata", "refGene.csv", package="SeqData")
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           geneAnnot=.geneAnnot(featureAnnot)
                           featureAnnotation(obj)=
                               resize(geneAnnot,fix="start",width=1)
                       },
                       "5LTR"={
                           annotationFile(obj)=system.file(
                               "extdata", "rmasker_repClassLTR.csv", 
                               package="SeqData")
                           featureAnnot=read.csv(file=annotationFile(obj),
                                                 as.is=T,header=T) 
                           cat("Setting up featureAnnotation slot...","\n")
                           ervAnnot=.featureAnnot(featureAnnot)
                           featureAnnotation(obj)=
                               resize(ervAnnot,fix="start",width=1)
                       },                            
{stop(cat("no match was found for feature requested, please check spelling.\n"))}
                )}
return(obj)   
    }) 

# a dispatcher when input is annotationFile
setMethod(
    f="getFeatureAnnotation",
    signature="missing",
    definition=function(annotationFile=character(0),feature="gene"){
        
        obj=new("SeqData")
        
        annotationFile(obj)=annotationFile
        
        obj=getFeatureAnnotation(obj,annotationFile=annotationFile) 
        
        return(obj)
    })

##-----------------------------------------------------------------------------
## helperFunctions

.exonAnnot=function(refGene){
    
    refGene=.duplicationCheck(refGene)
    
    # indexing exons for each gene    
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
        gene=rep(refGene$name2,refGene$exonCount),
        exon=paste(rep(refGene$name,refGene$exonCount),exon.index,sep=":"),
        transcript=rep(refGene$name,refGene$exonCount))
    
    # corce RangedData to GRanges object
    exonAnnot=as(geneAnnot,"GRanges")
    return(exonAnnot)
}

# GRanges implementation
.geneAnnot=function(featureAnnot){
    featureAnnot=.duplicationCheck(featureAnnot)
    featureAnnot=GRanges(
        seqnames=featureAnnot$chrom,
        ranges=IRanges(
            start=featureAnnot$txStart,
            end=featureAnnot$txEnd),
        strand=featureAnnot$strand,
        symbol=featureAnnot$name2,
        EntrezID=featureAnnot$name)
    return(featureAnnot)     
    
}


.featureAnnot=function(featureAnnot){
    featureAnnot=.duplicationCheck(featureAnnot)
    featureAnnot=GRanges(
        seqnames=featureAnnot$chrom,
        ranges=IRanges(
            start=featureAnnot$txStart,
            end=featureAnnot$txEnd),
        strand=featureAnnot$strand,
        symbol=featureAnnot$name)
    return(featureAnnot)    
}

# duplication check
# sometimes exactly same item appears in UCSC refGene table, which should be removed before processing

.duplicationCheck=function(featureAnnot){
    chr.start=paste(featureAnnot$chrom,featureAnnot$txStart,sep="_") 
    uniGene=paste(chr.start,featureAnnot$name)
    featureAnnot$names=uniGene
    featureAnnot=featureAnnot[!duplicated(featureAnnot$names),]
    featureAnnot$names=NULL
    return(featureAnnot)
    
}

##-----------------------------------------------------------------------------
## LATER:
# getAnnotationFromUCSC, feature="tss"

## to do:

# remove the obj
# getFeatureAnnotation(obj,annotationFile=annotationFile)
# getFeatureAnnotation(obj,feature="gene")





