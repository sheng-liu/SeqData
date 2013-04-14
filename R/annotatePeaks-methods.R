
# this method should be the last step of peak finding, it's kind of annotation of peaks of interest

# maybe I shold name it annotatePeaks


## annotatePeaks
# annotate peaks with the provided or internal annotation table (tss for gene, 5' tss for erv. Annotation only needs tss. it's a singe point, not a range. peak postion is also a point, all these can be change, but for simplicity, this is best, at least for now). proivde the targetTable with peaks which are located close to tss within 200bps. 


# get the distance to tss and plot it.


## output

## a targetTable
## set the peakRegions(obj) to the target regions. decided not to do this, so you can annotate the same peak list with different annotations.but output a targetTable instead. 

##'@name annotatePeaks
##'
##'@title annotatePeaks



##'@description
##' get peak target genomic features, how many tss is enriched with this histone mark 
##' what are the target genes of this protein of interest.
##' 
##'this method 
##'use annotation->tss positon


##' peakPosition   tss position (subject)  overlap position
##' function "nearest" tells which subject has overlap with the query, 
##' select the positive tss positions.
##' 


setMethod(
    f="annotatePeaks",
    signature="SeqData",
    definition=function(obj){
        
        # getFeatureAnnotation(,feature="tss")
     
        # resize both to length 1?
        
        #peakRanges=IRanges(viewWhichMaxs(peakRegions(obj)),width=1)
        
        # making RleViewList peakRegions(obj) a GRanges with chromosome names
        # the only information you need is the ranges inside this RleViewList for calculating the nearst (ps. it needs to be seperated with chromosome). so first step is to linearize the list. either using IRanges with space, or using GRanges. 
        #peakRanges=as(peakRegions(obj),"GRanges")  # GRanges class
        # this however creaates lots of ranges more than actual
        obj=getFeatureAnnotation(obj,annotationFileName,feature="tss")
        
        peakIRangesList=ranges(peakRegions(obj))  #SimpleRangesList class
        
        peakGRanges=as(peakIRangesList,"GRanges")
        
        ##  you can create peakGRanges from the viewWhichMaxs(peakRegions)
        
        #peakGRanges=GRanges(viewWhichMaxs(peakRegions(obj))
        
        #peakRanges.rd=as(peakRegions(obj),"RangedData")     # RangedData class
        
        #peakRangesList=ranges(peakRanges.rd)  #SimpleRangesList class
        
        ##### coerce RangesList, which can be feedinto view to calculate read counts.
        
        
        featureGRanges=featureAnnotation(obj)  # GRanges class
        
        #featureRanges=as(featureRanges,"RangedData")   # RangedData class
         
        #featureRangesList=ranges(featureRanges)  # RangesList class
      
        # may need an lapply here to do it on multiple chromosome, no, as it is GRanges list which has chromosome already  check
        
        #suppressWarnings()
        distanceToNearest=suppressWarnings(distanceToNearest(peakGRanges,featureGRanges))
        ## GRanges::: distanceToNearest accept GRanges, which has the chromosome information
        ## IRanges:::distanceToNearest accept ranges(IRanges), not even any list, which doesn't have chromosome information, which maens need a loop to work on genome. 
        
        ####### Not sure wehter this function consider chromosome 
        ####### if not, then the whole thing needs to be inside a loop to get it work right for multiple chromosomes. The GRanges has chromosome information thought, not sure whether they have used it to seperate the ranges or just mix them all as a whole, which can be a problem when mix differen chromosome.
        
        ######## while it considers chromosome yes and no.
        ########  it provide a list with each index of  query hits. 
        # the query is what ever chromosome
        ## there might be some coordinates are same on different chromosoemes, when they are presented as GRanges, they have the seqnames(space/list name if IRanges) to seperate them. but if they are put together, there may have dupicates. I think they must have handled this. but to the 
        ## when you want to know what chromosome they are in, or use chromosome to seperate them and output them, you will need to look at the corresponding quirey hits, where are they located. 
        
        
        
        distanceToFeature=distanceToNearest$distance
        featureIndex=distanceToNearest$subjectHits
        
        #peakIndex=distanceToNearest$queryHits
        
        
#-----------------------------------        
#         

        
        

#         featureInfo=mcols(featureAnnotation(obj))
#         
#         targetFeatureNames=featureInfo$symbol[featureIndex]
#         
#         targetFeatureIDs=featureInfo$EntrezID[featureIndex]
#         
#         targetFeaturePos=start(featureAnnotation(obj))[featureIndex]
#                                
# 
#         targetFeatures=as.data.frame(cbind(targetFeatureNames,targetFeatureIDs,targetFeaturePos,distanceToFeature))
        
        ## featureIndex has no use, what's make sense is the target Index, which is a subset of feature Index, which should be output to the user.
        
        
        
        
        #peakTable needs to be subsetted 
        
        #peakTable=cbind(peakTable,targetFeatures)
        
        # targetIndex=abs(distanceToFeature)<=200  # this will be ideal
        targetIndex=abs(distanceToFeature)<=2000   # for now
        ## this parameter can kill hundres targets  
        
        
        
        #distanceToFeature[targetIndex]
        
        
        # all these need to be list to make the rest subseting work
        
        
        
        ### use the calculated ranges subset coverage
        
        # views are only in IRanges, it is subsetting coverage into small ranges.
        
        # when in IRanges, deal with chromosome usually are in lists, or one can use space. 
        
        #as(peakRanges,"RangedData")
        
        chrl=names(chrSize(obj.chip))
        
        
        #peakRangedDataList=split(peakRanges,space(peakRanges))
        
        #peakRangesList=ranges(peakRangedDataList)
        
        #Views(readCoverage(obj),peakRangesList) # worked
        # it has to be RangesList, not GRangesList
        
        # get targetRangesList
        
#         > Views(readCoverage(obj),targetGRangesList)
#         Error in RleViewsList(rleList = subject, rangesList = start) : 
#             'rangesList' must be a RangesList object
        
        
        targetGRanges=peakGRanges[targetIndex]
        #targetRanges=peakRanges.rd[targetIndex]
        
        
        
        targetGRangesList=split(targetGRanges,seqnames(targetGRanges))
        
        
        ## TODO: instead of create new views, how about just subset peakViews already have
        
        targetRegions=RleViewsList(sapply(chrl,function(chr) {
            Views(
                readCoverage(obj)[[chr]],
                start=start(targetGRangesList)[[chr]],
                end=end(targetGRangesList[[chr]]))
        }))
        
        ## peakRegions(obj)=targetRegions
        ## leave it there for now.
        
#         targetRegions.l=(sapply(chrl,function(chr) {
#             Views(
#                 readCoverage(obj)[[chr]],
#                 start=start(targetGRangesList)[[chr]],
#                 end=end(targetGRangesList[[chr]]))
#         }))
         # this gives a list, not compact as RleViewsList which can be used the same
        
       
            
            # peakRegions  a view on coverage for calculation of readCounts
            # - regions that have higher count in chip than control
            # - peakRegions after findPeakSummits
            # - targetRegions  after findTargetFeatures
            
        
        
        #targetRegions=peakRegions(obj)[ranges(peakRegions(obj))==ranges(targetRanges)]
        
        
        #x=peakRegions(obj)[["chr12"]][targetIndex]
        
        ## creating a view on the coverage is pretty fast, so no need for subseting it every time, just create a new view. 
        
#         lapply(chrl,function(chr){
#             
#             peakRegions(obj)[[chr]][targetIndex]
#             
#             
#         })
        
        # but how do you know these targetRanges or targetRegions are in the significant peak Regions?
        # because targetRanges are the results of subsetting peak regions with featureRegions, so it is inside of the peaklist. 
        
        
        # output viewSums for total reads, viewMax for peak height, viewWhichMax for peakPosition, put this with the feature name into one targetFeatureTable
        
        
        
        
        
        
        
        ## plot "distribution of peaks relative to tss"
        
        # subsetting the views to get the reads count on regions
        
        # by swapping to nearest(features, peakRegions)
        # I can answer which peakRegions are overlap/nearest to the features
        # targetRegions=unlist(ranges(peakRegions(obj)))[targetIndex]
        # it turn out the target retions are all the regions, there is no selection yet
        # until you did the subseting
        
        
        
        # summarizeOverlaps and this are answering two different questions
        # sum
        # 
        # this is you have all the peaks, and througth this, you find those peaks that are close to the gene tss within 200bp, and it gives you the read counts of the whole peak. 
        
        # summarizeOverlaps is you have the coverage and you have the ranges of interest, you use it can assay, ok, how many reads is there at the upstream downstream 200bp of TSS regions, and you can have a list of enrichement on all tss, range from huge to nothing. 
        # the difference summarize is cutting coverage with ranges. above method is keeping the peak intact but found those peaks that have a TSS nearby (targeting tss). 
        
        
        # peakRegions are RleViews

# it has the coverage information, so you don't need countReadsOverFeatures to count the reads in the regions, subletting the region with ranges of interest should give you the reads count (coverage), then you can fill it in the peakTable as RPKM (either on each peak region, or on gene region, or any feature regions,subset the regions should do that)

         # one key difference in this two is, peakRegions doesn't have overlaps, you get those by slicing the peaks, so there is no need for dealing with where to put the reads that overlaps on two peaks. this is the case for RNA seq reads, where I expect 60~80% of the gene annotations are overlaped. I presume you can use countOverlaps to count reads in peakRegions as long as you give the ranges and aln, and proved the counting mode. but it is not as simple or intuitive as just using views you already have.       
        
        
        
        
        
        
        
        
#         targetRegions=resize(featureAnnotation(obj),width=1000,fix="center")
#         
#         targetReadCount=viewSums(Views(peakRegions(obj),targetFeatureRegions))
        
          
        
        # leave peakRegions there for now as it is hard to reset back once it is set to smaller scales.
        #targetReadCount=viewSums(peakRegions(obj))
        #targetPeakPos=viewWhichMaxs(peakRegions(obj))
        
        
        targetReadCount=viewSums(targetRegions)
        targetPeakPos=viewWhichMaxs(targetRegions)
        
        
#         
#         peakTable=cbind(peakTable,targetReadCount)
        
        
        targetPeakReadCount=unlist(targetReadCount)
        targetPeakPos=unlist(targetPeakPos)
        
        
  #  ____________________________________________________________
   # these need to be in a loop as it seems doesn't out put chromosome information
        # it's a genomewide list, all the sequence information is corresponding to the query, which is the peak location. don't forget to subset them twice first with featureIndex then with targetIndex, that is the one that is picked for the output in the end.
#         chrl=names(libSize(obj))
#         lapply(chrl,function(chr){
#             
#             
#             
#             
#         }
                #featureInfo=mcols(featureAnnotation(obj))
        
        featureInfo=mcols(featureGRanges[featureIndex])
                
                targetFeatureNames=featureInfo$symbol[targetIndex]
        #featureInfo$symbol[featureIndex] returns all feature
        # targetIndex returns only the target features
                
                targetFeatureIDs=featureInfo$EntrezID[targetIndex]
                
                #targetFeaturePos=start(featureAnnotation(obj))[targetIndex]
        targetFeaturePos=start(featureGRanges[featureIndex])[targetIndex] 
                  
                peakDistanceToFeature=distanceToFeature[targetIndex]
        
 #  ____________________________________________________________            
        
                targetTable=as.matrix(cbind(targetFeatureNames,targetFeatureIDs,targetFeaturePos,targetPeakPos,peakDistanceToFeature,targetPeakReadCount))
        
          cat("Output file...","\n")
          write.csv(file="targetTable.csv",targetTable)
          return(targetTable)
        
        
    })