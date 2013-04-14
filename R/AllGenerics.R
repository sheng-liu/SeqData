# SeqGardgets - Generics
# 
# 
###############################################################################

# Generics for Methods

setGeneric(
    name="viewReadCoverage",
    def=function(obj=NULL,bamFileName=character(1)){
        standardGeneric("viewReadCoverage")
    })

# TODO: set viewReadCoverage also like countReadsOverFeatures

setGeneric(
    name="getReadCoverage",
    def=function(obj,bamFileName){
        standardGeneric("getReadCoverage")
    })


setGeneric(
    name="getReadAlignment",
    def=function(obj,bamFileName){
        standardGeneric("getReadAlignment")
    })

setGeneric(
    name="getFeatureAnnotation",
    def=function(obj,annotationFileName,feature=c("gene","tss","exon","transcript")){
        standardGeneric("getFeatureAnnotation")
    })

setGeneric(
    name="countReadsOverFeatures",
    def=function(obj=NULL,bamFileName=character(1),annotationFileName=character(1),countingMode="Union"){
        standardGeneric("countReadsOverFeatures")
    })

setGeneric(
    name="findPeaks",
    def=function(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(1),bamFile.control=character(1)){
        standardGeneric("findPeaks")
    })



##------------------------------------------------------------------------------
## Generics for Accessors and Setters
##

setGeneric(
    name="chrSize<-",
    def=function(obj,value){
        standardGeneric("chrSize<-")
    })

setGeneric(
    name="chrSize",
    def=function(obj){
        standardGeneric("chrSize")
    })

setGeneric(
    name="libSize<-",
    def=function(obj,value){
        standardGeneric("libSize<-")
    })

setGeneric(
    name="libSize",
    def=function(obj){
        standardGeneric("libSize")
    })



setGeneric(
    name="readAlignment<-",
    def=function(obj,value){
        standardGeneric("readAlignment<-")
    })

setGeneric(
    name="readAlignment",
    def=function(obj){
        standardGeneric("readAlignment")
    })

setGeneric(
    name="readCoverage",
    def=function(obj){
        standardGeneric("readCoverage")
    })

setGeneric(
    name="readCoverage<-",
    def=function(obj,value){
        standardGeneric("readCoverage<-")
    })


setGeneric(
    name="featureAnnotation<-",    
    def=function(obj,value){
        standardGeneric("featureAnnotation<-")
    })

setGeneric(
    name="featureAnnotation",
    def=function(obj){
        standardGeneric("featureAnnotation")
    })



setGeneric(
    name="peakRegions",
    def=function(obj){
        standardGeneric("peakRegions")
    })

setGeneric(
    name="peakRegions<-",
    def=function(obj,value){
        standardGeneric("peakRegions<-")
    })

setGeneric(
    name="annotatePeaks",
    def=function(obj,bamFileName){
        standardGeneric("annotatePeaks")
    })