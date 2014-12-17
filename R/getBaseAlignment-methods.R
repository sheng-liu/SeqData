# getBaseAlignment
# 
# 
###############################################################################

## caution, it creats a new obj, which means it erases the inforamtion that was in the obj

## baseAlignment,
#alignment at individual, separated bases. It is usually two or more different read counts.  e.g. methylated vs unmethylated, allilic difference, A vs T (single nucleotide polymorphism) which is base Alignment of A and T at that location. 

# the baseReport is usually no header, in tab seperated format, it is user's repsonsibility to assign header to the table and save it as csv file

# get cytosine methylaition from aligned bamFile and form methylation report

# read in cytosine methylation from methylationFile
# add filter on read coverage, >=5, 
# or people and do it after, it is eay to sort that way

## TODO, for two couple data sets, need to find the overlap

# methylatonFile is output from DNA methylation aligners such as bismark, in the
# form of a csv file has below columns "chromosome", "position", "strand",
# "count_methylated", "count_unmethylated", "C_context", "trinucleotide_context"

# a dispatch for methylation report
setMethod(
    f="getBaseAlignment",
    signature=c(baseReport="character"),
    definition=function(obj=NULL,bamFile=character(0),baseReport){
        
        obj=SeqData()
        # this erased all information
        # none of the functions should have this statement
        # SeqData should be announced before these function, all functions only
        # add slots instead of erase, this way it is simply to use
        # getReadAlignment does the same, may change this in the future
        
        cat("Reading in baseReport...\n")
        # fill=T in case there is mising values, add NA
        report.df=read.csv(file=baseReport,as.is=T,header=T,fill=T)
        # temporarily or for new use tab
        
        cat("Filtering NAs...\n")
        # remove rows that has missing value on coverage
        cov=report.df$count_methylated+report.df$count_unmethylated
        report.df=report.df[!is.na(cov),]
        
        # complete.cases()
        # report.df=report.df[complete.cases(report.df),]
        
        
        # add strand if there isn't 
        cln=colnames(report.df)
        if (length(which(cln=="strand"))!=0) strand=report.df$strand 
        else strand=rep("*",dim(report.df)[1]) 
        
        
        cat("Format chromosome names...\n")
        # chr MT
        report.df$chromosome=.convertChrNames(
            report.df$chromosome,chr="chr",MT="MT")
        
        
        cat("Constructing data.frame...\n")
        base.df=with(report.df,
                        data.frame(chr=chromosome,
                                   start=position,
                                   end=position,
                                   strand=strand,
                                   meth_count=count_methylated,
                                   unmeth_count=count_unmethylated))
        
        cat("Constructing baseAlignment...\n")
        base.gr=df2gr(base.df)  # this step takes long
        baseAlignment(obj)=base.gr
        
        cat("Complete getBaseAlignment.\n")
        return(obj)
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


# it is essentially a methylation extractor
# It seems bioString, BSgenome, GAlignment package can do the work well

# a dispatch for bismark alignment bam file
setMethod(
    f="getBaseAlignment",
    signature=c(bamFile="character"),
    definition=function(obj=NULL,bamFile,baseReport=character(0)){
        print("Support for bismark aligned bam file is in development, currently please use cytosine report as input file.")
    })

# a dispatch for SeqData
setMethod(
    f="getBaseAlignment",
    signature=c(obj="SeqData"),
    definition=function(obj,bamFile=character(0),baseReport=character(0)){
        print("Support for bismark aligned bam file is in development, currently please use cytosine report as input file.")
    }
    
)

