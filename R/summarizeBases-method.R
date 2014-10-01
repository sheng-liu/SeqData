# summarizeBases
# 
# 
###############################################################################

## count base measurement (e.g methylation) on any genomic regions

## perc_meth = meth_counts/(meth_counts+unmeth_counts) within each region
## the method can be applied to any parameters measured along genome. essentially, the function is to group the items in bases based on region provided, then comput whatever statistics over those groups
## regions, needs to be in GRanges (are in GRanges)
## bases, needs to be in GRanges (are in data.frames)
## as the grouping is easier done by GRangs facility findOverlaps



setMethod(
    f="summarizeBases",
    signature=c(obj="SeqData"),
    definition=function(obj,description=character(0)){
        
        regions=featureAnnotation(obj)
        bases=baseAlignment(obj)
        
        baseTable=.summarizeBases(regions,bases,ignore.strand=T)
        
        #output
        cat("output baseSummarizeTable...\n")
        fileName=paste("baseTable-",description,"-",.timeStamp(bamFile(obj)),".csv",sep="")
        write.csv(file=fileName,baseTable)
        cat("\nDone!\n")
        
        return(obj)
        
    })


setMethod(
    f="summarizeBases",
    signature=c(bamFile="character",annotationFile="character"),
    definition=function(bamFile,annotationFile,description=character(0)){
        
        print("Support for bismark aligned bam file is in development, currently please use cytosine report as input file.")
        
        #sd=SeqData(bamFile=bamFile,annotationFile=annotationFile)
        #regions=featureAnnotation(getFeatureAnnotation(annotationFile=annotationFile))
        
        
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
    cat("Summarizing on groups...\n")
    region.summarize=ol.df %.%
        group_by(element) %.%
        select(element,width,meth_count,unmeth_count) %.%
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
    region.summarize.out=merge(regions.df,region.summarize,by.x="id",by.y="element")
    
    region.summarize.out$id=NULL
    names(region.summarize.out)[1]="chr"
    
    
    cat("Done!\n")
    return(region.summarize.out)
    
    
}

