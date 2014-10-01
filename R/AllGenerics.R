# SeqData - AllGenerics
# 
# 
###############################################################################




# Generics for Methods
##' @exportMethod viewCoverage
setGeneric(
    name="viewCoverage",
    def=function(obj=NULL,obj.control=NULL, bamFile=character(0),bamFile.control=character(0),weight.control=1,feature=character(0),ranges=GRanges(),annotationFile=character(0),split.metaData=F,bigWig=F,description=character(0)){
        standardGeneric("viewCoverage")
    })

##' @exportMethod getReadCoverage
setGeneric(
    name="getReadCoverage",
    def=function(obj,bamFile=character(0)){
        standardGeneric("getReadCoverage")
    })

##' @exportMethod getReadAlignment
setGeneric(
    name="getReadAlignment",
    def=function(obj,bamFile=character(0),annotationFile=character(0)){
        standardGeneric("getReadAlignment")
    })

##' @exportMethod getFeatureAnnotation
setGeneric(
    name="getFeatureAnnotation",
    def=function(obj=NULL,annotationFile=character(0),feature="gene"){
        standardGeneric("getFeatureAnnotation")
    })


##' @exportMethod countReads
setGeneric(
    name="countReads",
    def=function(obj=NULL,bamFile=character(0),annotationFile=character(0),countingMode="Union",interFeature=TRUE,feature=character(0),description=character(0)){
        standardGeneric("countReads")
    })


##' @exportMethod findPeaks
setGeneric(
    name="findPeaks",
    def=function(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(0),bamFile.control=character(0),fdr.max=1e-5,foldChange.min=2,description=character(0),local.lambda=T){
        standardGeneric("findPeaks")
    })

##' @exportMethod getBaseAlignment
setGeneric(
    name="getBaseAlignment",
    def=function(obj,bamFile=character(0),baseReport=character(0)){
        standardGeneric("getBaseAlignment")
    })

##------------------------------------------------------------------------------

# Generics for Accessors and Setters

##'@exportMethod chrSize<-
setGeneric(
    name="chrSize<-",
    def=function(obj,value){
        standardGeneric("chrSize<-")
    })

##'@exportMethod chrSize
setGeneric(
    name="chrSize",
    def=function(obj){
        standardGeneric("chrSize")
    })

##' @exportMethod libSize<-
setGeneric(
    name="libSize<-",
    def=function(obj,value){
        standardGeneric("libSize<-")
    })

##' @exportMethod libSize
setGeneric(
    name="libSize",
    def=function(obj){
        standardGeneric("libSize")
    })

##' @exportMethod readAlignment<-
setGeneric(
    name="readAlignment<-",
    def=function(obj,value){
        standardGeneric("readAlignment<-")
    })

##' @exportMethod readAlignment
setGeneric(
    name="readAlignment",
    def=function(obj){
        standardGeneric("readAlignment")
    })

##' @exportMethod readCoverage
setGeneric(
    name="readCoverage",
    def=function(obj){
        standardGeneric("readCoverage")
    })

##' @exportMethod readCoverage<-
setGeneric(
    name="readCoverage<-",
    def=function(obj,value){
        standardGeneric("readCoverage<-")
    })

##' @exportMethod featureAnnotation<-
setGeneric(
    name="featureAnnotation<-",    
    def=function(obj,value){
        standardGeneric("featureAnnotation<-")
    })

##' @exportMethod featureAnnotation
setGeneric(
    name="featureAnnotation",
    def=function(obj){
        standardGeneric("featureAnnotation")
    })


##' @exportMethod coverageView
setGeneric(
    name="coverageView",
    def=function(obj){
        standardGeneric("coverageView")
    })

##' @exportMethod coverageView<-
setGeneric(
    name="coverageView<-",
    def=function(obj,value){
        standardGeneric("coverageView<-")
    })

##' @exportMethod annotatePeaks
setGeneric(
    name="annotatePeaks",
    def=function(obj=NULL,peakTable=character(0),annotationFile=character(0),feature=character(0),distance.max=4000,description=character(0)){
        standardGeneric("annotatePeaks")
    })

##' @exportMethod bamFile
setGeneric(
    name="bamFile",
    def=function(obj){
        standardGeneric("bamFile")
    })

##' @exportMethod bamFile<-
setGeneric(
    name="bamFile<-",    
    def=function(obj,value){
        standardGeneric("bamFile<-")
    })

##' @exportMethod annotationFile
setGeneric(
    name="annotationFile",
    def=function(obj){
        standardGeneric("annotationFile")
    })

##' @exportMethod annotationFile<-
setGeneric(
    name="annotationFile<-",    
    def=function(obj,value){
        standardGeneric("annotationFile<-")
    })

##' @exportMethod coverageView
setGeneric(
    name="coverageView",
    def=function(obj,ranges){
        standardGeneric("coverageView")
    })
##' @exportMethod intersectChr
setGeneric(
    name="intersectChr",
    def=function(reads,features){
        standardGeneric("intersectChr")
    })

##' @exportMethod baseAlignment<-
setGeneric(
    name="baseAlignment<-",
    def=function(obj,value){
        standardGeneric("baseAlignment<-")
    })

##' @exportMethod baseAlignment
setGeneric(
    name="baseAlignment",
    def=function(obj){
        standardGeneric("baseAlignment")
    })


