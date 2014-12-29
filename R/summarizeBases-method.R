# summarizeBases
# 
# 
###############################################################################
##' @name summarizeBases
##' @aliases summarizeBases
##' @title summarizeBases
##' @rdname summarizeBases-methods
##' @docType methods
##' @description count base measurement (e.g methylation) on any genomic regions. 
##' @usage 
##' summarizeBases(obj=NULL,bamFile=character(0),annotationFile=character(0),baseReport=charcter(0),description=character(0))

##' @param obj A SeqData object.
##' @param bamFile The full path to the sequencing data. The data needs to be in
##'   bam file format.

##' @param annotationFile csv annotation file for specifying ranges of interest.
##'   It must includes the miminum five field: 'chrom','txSart','txEnd','strand'
##'   and 'name'.
##' 
##' @param baseReport DNA methylaton report from methylaton aligners such as 
##'   bismark, in the form of a csv file has minimu four columns "chromosome", 
##'   "position", "strand", "count_methylated", "count_unmethylated". if strand information is not available, fill it with "*". 
##'   
##   Additionaly
##   "C_context", "trinucleotide_context". chr=chromosome,position,strand,count_methylated,count_unmethylated
##' 
##' @param description Optional. User can input short desciption for the output
##'   file, it will show up in the file name of the output file.
##' 
##' @details the method can be applied to any parameters measured along genome.
##'   essentially, the function is to group the items in bases based on region
##'   provided, then comput whatever statistics over those groups. Internally,
##'   it calls getBaseAlignment.
##'   


##' @return 
##' \itemize{
##' \item{baseTable.csv} a csv file contains the base summarization.
##' }
##' @section Usage:{summarizeBases(obj=NULL,bamFile=character(0),annotationFile=character(0),baseReport=charcter(0),description=character(0))}
##' @examples
##' # summarizeBases(obj)  
##' # summarizeBases(baseReport)
##' # summarizeBases(bamFile=bamFile,annotationFile=annotationFile)

## count base measurement (e.g methylation) on any genomic regions

## perc_meth = meth_counts/(meth_counts+unmeth_counts) within each region the
## method can be applied to any parameters measured along genome. essentially,
## the function is to group the items in bases based on region provided, then
## comput whatever statistics over those groups regions, needs to be in GRanges
## (are in GRanges) bases, needs to be in GRanges (are in data.frames) as the
## grouping is easier done by GRangs facility findOverlaps


# dispatch for baseReport
setMethod(
    f="summarizeBases",
    signature=c(baseReport="character",annotationFile="character"),
    definition=function(obj=NULL,bamFile=character(0),annotationFile,baseReport,description=character(0)){
        
        print("SummarizeBases")
        
        # alignment level function creates new obj, ie. erase obj information
        # need to put ahead
        obj=getBaseAlignment(baseReport=baseReport)
        
        #obj=SeqData(annotationFile=annotationFile)
        #obj=getFeatureAnnotation(obj)
        
        featureAnnotation(obj)=featureAnnotation(
            getFeatureAnnotation(annotationFile=annotationFile))
        ## change file level function to ouput just the slot maybe better idea
        ## avoid all these erase mistakes and long typing
        
        regions=featureAnnotation(obj)
        bases=baseAlignment(obj)
        
        baseTable=.summarizeBases(regions,bases,ignore.strand=T)
        
        #output
        cat("output baseTable...\n")
        fileName=paste("baseTable-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,baseTable,row.names=F)
        cat("\nDone!\n")
        
        return(obj)
        
    })


setMethod(
    f="summarizeBases",
    signature=c(obj="SeqData"),
    definition=function(obj,bamFile=character(0),annotationFile=character(0),baseReport=character(0),description=character(0)){
        
        
        # check if featureAnnotation and baseAlignment is filled
        if (length(featureAnnotation(obj))==0||length(baseAlignment)==0) {
            cat("Please fill in featureAnnotation or baseAlignment slot first.\n")
        }
        
        regions=featureAnnotation(obj)
        bases=baseAlignment(obj)
        
        baseTable=.summarizeBases(regions,bases,ignore.strand=T)
        
        #output
        cat("output baseTable...\n")
        fileName=paste("baseTable-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,baseTable)
        cat("\nDone!\n")
        
        return(obj)
        
    })


setMethod(
    f="summarizeBases",
    signature=c(bamFile="character",annotationFile="character"),
    definition=function(obj=NULL,bamFile,annotationFile,baseReport=character(0),description=character(0)){
        
        print("Support for bismark aligned bam file is in development, currently please use base report as input file.")
        
        
    })
    
    
.summarizeBases=function(regions,bases,ignore.strand=T){
    ##-----------------------------------------------------------------------
    ## grouping    
    # find the regions, bases overlaps
    cat("Finding overlaps between regions and bases...\n")
    query=bases
    subject=regions
    ol=suppressWarnings(findOverlaps(query,subject,ignore.strand=ignore.strand))
    
    # grouping
    cat("Grouping bases with regions...\n")
    ol.ls=split(query[queryHits(ol)],subjectHits(ol))
    
    # transition from GrangesList to data.frame (for faster calculation)
    ol.df=as.data.frame(ol.ls)
    
    ##-------------------------------------------------------------------------
    ## calculate mean of groups
    
    # colnames(ol.df)
    # R3.0          element seqnames start end width strand meth_count unmeth_count
    # R3.1 group group_name seqnames start end width strand meth_count unmeth_count    
    
    cat("Summarizing on groups...\n")
    region.summarize=ol.df %>%
        group_by(group_name) %>%
        select(group_name,width,meth_count,unmeth_count) %>%
        summarise(    					# British summarise
            num_cpg_site=sum(width), 	# num_cpg_site=table(element), # takes 11s
            methylated_count=sum(meth_count),
            unmethylated_count=sum(unmeth_count),
            perc_meth=methylated_count/(methylated_count+unmethylated_count)       	
        ) 
    
    ##------------------------------------------------------------------------
    ## associate subjectHit index with region information 
    # index regions 
    regions.df=as.data.frame(regions)
    regions.df=transform(regions.df,id=rownames(regions.df))
    
    # merge() only works on columns not rownames
    region.summarize.out=merge(regions.df,region.summarize,by.x="id",by.y="group_name")
    
    region.summarize.out$id=NULL
    names(region.summarize.out)[1]="chr"
    
    
    cat("Done!\n")
    return(region.summarize.out)
    
    
}

