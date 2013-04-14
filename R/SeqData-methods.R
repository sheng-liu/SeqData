#SeqData-method

## a method for create a SeqData project which enables user to use all the functions within the package and access all the data from the object, instead of just the output file. which also facilitate or simplified the function/method writing. all the functions/method only intake one thing, that is SeqData, and SeqData has the most basic information they need, the bamFile and the annotation file.

# not tested yet. still in idea.


## not sure why most of the functions in R takes the most basic "character", it is to simplify the usage of the function for the user, as user don't need to create the SeqData objects but just simply use it.
## but here if people create one, it make so much sense as it simplifies the coding and the usage as well, as user doesn't need to know what function to run first, what do they need to use this function, they just run it, and the function does what it needs to get it done.



setMethod(
    f="SeqData",
    signature="character",
    definition=function(bamFile=character(1),annotationFile=character(1)){
        obj=new("SeqData")
        obj@bamFile=bamFile
        obj@annotationFile=annotationFile
        return(obj)
    })

setMethod(
    f="SeqData",
    signature=c("character","missing")
    definition=function(bamFile=character(1)){
        obj=new("SeqData")
        obj@bamFile=bamFile
        return(obj)
        
    })

setMethod(
    f="SeqData",
    signature=c("character","missing")
    definition=function(annotationFile=character(1)){
        
        obj=new("SeqData")
        obj@annotationFile=annotationFile
        return(obj)
    })




setMethod(
    f="SeqData",
    signature="missing",
    definition=function(obj){
        
        
        
        
    })