# Accessors and Setters
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
    f="coverageView",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,coverageView=value)
    })

setMethod(
    f="coverageView",
    signature="SeqData",
    definition=function(obj){
        obj@coverageView
    })

setMethod(
    f="bamFile",
    signature="SeqData",
    definition=function(obj){
        obj@bamFile
    })

setReplaceMethod(
    f="bamFile",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,bamFile=value)
    })

setMethod(
    f="annotationFile",
    signature="SeqData",
    definition=function(obj){
        obj@annotationFile
    })

setReplaceMethod(
    f="annotationFile",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,annotationFile=value)
    })

setMethod(
    f="baseAlignment",
    signature="SeqData",
    definition=function(obj){
        obj@baseAlignment
    })

setReplaceMethod(
    f="baseAlignment",
    signature="SeqData",
    definition=function(obj,value){
        initialize(obj,baseAlignment=value)
    })


