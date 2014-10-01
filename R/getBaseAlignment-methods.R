# getBaseAlignment
# 
# 
###############################################################################

## baseAlignment,
#alignment at individual, separated bases. It is usually two or more different read counts.  e.g. methylated vs unmethylated, allilic difference, A vs T (single nucleotide polymorphism) which is base Alignment of A and T at that location. 


# a dispatch for methylation report
# get cytosine methylaition from aligned bamFile and form methylation report

# read in cytosine methylation from methylationFile
# add filter on read coverage, >=5, 
# or people and do it after, it is eay to sort that way

## TODO, for two couple data sets, need to find the overlap

# methylatonFile is output from DNA methylation aligners such as bismark, in the form of a csv file has below columns
# "chromosome", "position", "strand", "count_methylated", "count_unmethylated", "C_context", "trinucleotide_context"


setMethod(
    f="getBaseAlignment",
    signature=c(baseReport="character"),
    definition=function(baseReport){
        
        report.df=read.csv(file=baseReport,as.is=T,header=T)
        base.df=with(report.df,
                        data.frame(chr=chromosome,
                                   start=position,
                                   end=position,
                                   strand,
                                   meth_count=count_methylated,
                                   unmeth_count=count_unmethylated))
        
        base.gr=df2gr(base.df)  
        
        baseAlignment(obj)=base.gr
    }
    
)

## TODO:
## df2gr has a species issue, now it is all set to mouse
## may need to add in human and other species for GRanges to be correct

# species=mouse-mm9

# input:
# chr 0 
# MT 0 

# output:
# chr.out = 1 
# MT.out = 1 





# a dispatch for bismark alignment file
# it is essentially a methylation extractor
# It seems bioString, BSgenome, GAlignment package can do the work well


setMethod(
    f="getBaseAlignment",
    signature=c(bamFile="character"),
    definition=function(obj){
        print("Support for bismark aligned bam file is in development, currently please use cytosine report as input file.")
    }
    
)

setMethod(
    f="getBaseAlignment",
    signature=c(obj="SeqData"),
    definition=function(obj){
        print("Support for bismark aligned bam file is in development, currently please use cytosine report as input file.")
    }
    
)

