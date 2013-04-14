# SeqGardgets-Accessors and Setters
# 
# 
###############################################################################


setReplaceMethod(
    f="chrSize",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,chrSize=value)    
    })

setMethod(
    f="chrSize",
    signature="SeqData",
    definition=function(obj){
        obj@chrSize
    })

setReplaceMethod(
    f="libSize",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,libSize=value)
    })

setMethod(
    f="libSize",
    signature="SeqData",
    definition=function(obj){
        obj@libSize
    })

setReplaceMethod(
    f="readCoverage",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,readCoverage=value)
    })

setMethod(
    f="readCoverage",
    signature="SeqData",
    definition=function(obj){
        obj@readCoverage
    })



setReplaceMethod(
    f="readAlignment",
    signature="SeqData",
    definition=function(obj,value){   
        initialize(obj,readAlignment=value)
    })

setMethod(
    f="readAlignment",
    signature="SeqData",
    definition=function(obj){
        obj@readAlignment
    })

setReplaceMethod(
    f="featureAnnotation",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,featureAnnotation=value)
    })

setMethod(
    f="featureAnnotation",
    signature="SeqData",
    definition=function(obj){
        obj@featureAnnotation
    })

setReplaceMethod(
    f="peakRegions",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,peakRegions=value)
    })

setMethod(
    f="peakRegions",
    signature="SeqData",
    definition=function(obj){
        obj@peakRegions
    })
