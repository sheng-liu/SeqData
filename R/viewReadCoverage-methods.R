# viewReadCoverage-methods
# 
# 
###############################################################################

# Gardgets
# viewReadCoverage()
# countReadsPerGene()
# association()



# all methods
# getReadCoverage()    Calculate base pair resolution read coverage
# viewReadCoverage()    # View gene coverage    
# getReadAlignment()    # get readAlignment and set readAlignment slot.
# getGeneAnnotation()    # get geneAnnotation and set geneAnnotation
# countReadsPerGene()    # countReadsPerGene


# Naming Conventions: 
# functionName (functionVariableReceiverName), 
# variable.name 
# ConstantName, ClassName

# Naming Conventions -- Java naming coventions (pick this one for biocondutor)
# functionName, verbs
# variableName, noun 
# ClassName, PackageName
# CONSTANT

# using "." indicate belongs to, affliation, such as peaks.chip and peaks.control



# TODO: add a summmary of percentage of coverage of the genome
# make "SeqData" contain/inherit "RNAseq" class

# TODO: make viewReadCoverage() and countReadsPerGene() can accept obj, fileName, so that it can be used inside of other program, not just the start. 

## the word "coverage" in sequencing field is kind of confusing as it can mean two things, one is the how much of the genome is covered, the other is how many reads is there within a specific range. The first one is more intuitive, while the word for the second meanning should be changed to read density, how much reads in a specific range, which are more intuitive. althought when most people speak coverage, they are still talking about the reads density at a range.
# so maybe the name of this function should change to readDensity?

##'@title viewReadCoverage


#' @description this is a package for view read coverage

# @name viewReadCoverage method

#' @rdname viewReadCoverage-method
#' 
#' @param obj An {linkS4class} object
#' 
#' @param bamFileName The full path to the sequencing data in bam file format.
#' 
#' @return An S4 object. the slot readCoverage contains a contains a SimpleRleList composing of chromosome coverages as its elements.
#' 
# @section Details:  convert bam file into read density file
#' @details convert bam file into read density file
#' 
#' @author Sheng Liu
#' 
#' @seealso countReadsPerGene
#' 
#' @keywords methods
#' @examples
#' 
#' viewReadCoverage(bamFileName=  )  
#' viewReadCoverage()   
#' viewReadCoverage(obj) 
#  FIXME: the above two is not working well after changing. 




setMethod(
    f="getReadCoverage",
    signature="SeqData",
    definition=function(obj,
                        bamFileName=character(1)    	
    ){
        ## FIXME:
        ## read in as single end for both single- and pair- end data 
        ## read in as pair-end just if it is pair-end to get it's isize
        ## get mean(isize) for pair-end, and estimat.fragment.len for single end
        ## extend reads accordingly 
       
        obj=getReadAlignment(obj,bamFileName)        
        
        
        # get meanfragmentSize
        # extend the reads according to isize
        isize=mcols(readAlignment(obj))
        
        meanFragmentSize=round(mean(abs(isize[,1]),na.rm=T))
        
        if (meanFragmentSize==0){
            cat("SeqData is likely to be single-end, estimating meanFragmentSize...")
            meanFragmentSize=round(.estimateFragmentSize(readAlignment(obj)))
                                   
                                   
            cat("estimated meanFragmentSize is ",meanFragmentSize,"\n")
        }else{
            cat("SeqData is pair-end, meanFragmentSize is ",
                meanFragmentSize,"\n")
        }
        
        #fragmentExtendLength=meanFragmentSize-round(mean(width(readAlignment(obj))))
        
        aln <- as(readAlignment(obj), "GRanges")
        
        #assign("last.warning", NULL, envir = baseenv())
        #aln.extended= resize(aln, width=meanFragmentSize)
        aln.extended= suppressWarnings(resize(aln, width=meanFragmentSize))
        
        
        
        
        # compute bp coverage 
        cat("Computing read coverage...","\n")
        
        #readCoverage(obj)=coverage(readAlignment(obj),width=seqlengths(readAlignment(obj)),extend=as.integer(fragmentExtendLength))
        
        readCoverage(obj)=coverage(aln.extended,width=seqlengths(readAlignment(obj)))
                
        return(obj)
                
    }
)


setMethod(
    f="viewReadCoverage",
    signature="SeqData",
    definition=function(obj,bamFileName=character(1)){
        
        
        # get read coverage        
        obj=getReadCoverage(obj,bamFileName) 
        
        # smoothing with runing window average
        smoothedCoverage=runmean(readCoverage(obj),k=51,endrule="constant") 
        
        # export file for viewing
        cat("Output bigWig file...","\n")
        file.name=unlist(strsplit(bamFileName,split="/"))
        file.name.bigWig=paste(file.name[length(file.name)],".bigWig",sep="")    
        export(smoothedCoverage,file.name.bigWig)
        
        cat("Done.","\n\n")

        return(obj)
        
    }
)

setMethod(
    f="viewReadCoverage",
    signature="missing",
    function(obj, bamFileName){
        bamFileName=readline("\nPlease input the path to your bam file 
                          or simply drag the file to this window: ")
        # bamFileName is the path to the file plus the file name itself
        obj=new("SeqData")
        obj=viewReadCoverage(obj,bamFileName)         
        return(obj)

    })



