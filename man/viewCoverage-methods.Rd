\docType{methods}
\name{viewCoverage}
\alias{viewCoverage}
\title{viewCoverage}
\arguments{
  \item{obj}{A SeqData object.}

  \item{obj.control}{Optional. A SeqData object used as a
  control to see differential coverage with obj.}

  \item{bamFile}{The full path to the sequencing data. The
  data needs to be in bam file format.}

  \item{bamFile.control}{Optional. A bam file used as a
  control to see differential coverage with bam file.}

  \item{weight.control}{A parameter used to weight control
  data (obj.control or bamFile.control).}

  \item{feature}{interal feature annotation data set,
  includes: feature=c("gene","erv","repeats","transposon")}

  \item{ranges}{GRanges object specifying the ranges the
  coverageTable will produce.}

  \item{annotationFile}{csv annotation file for specifying
  ranges of interest. It must includes the miminum five
  field: 'chrom','txSart','txEnd','strand' and 'name'.}

  \item{split.metaData}{information in GRanges,or
  annotation file other than
  'chrom','txSart','txEnd','strand' are treated as
  metaData. if split.metaData=T, the ranges will be splited
  based on the first column of metaData, and the coverage
  will be calculated as the total coverage on all ranges
  that have the same metaData. It can be used in
  conjunction with parameter "feature","ranges" and
  "annotationFile".}

  \item{bigWig}{a switch for output bigWig file. The data
  is from readCoverage slot of SeqData.}

  \item{description}{Optional. User can input short
  desciption for the output file, it will show up in the
  file name of the output file.}
}
\value{
  \itemize{ \item{bamFileName.bigWig} a bigWig file for
  viewing the coverage of the bamFile provided
  \item{bamFileName.csv} a csv file contains the coverage
  information about the GRanges provided.
  \item{coverageView slot} if ranges are passed in, a
  RleViews object are returned and coverage corresponding
  to each of those ranges are stored in the coverageView
  slot. }
}
\description{
  Generate a bigWig file for viewing read coverage in a
  genome browser, or a coverageTable.csv for genomic ranges
  specified.
}
\details{
  Generate bigWig file from a bam file for viewing read
  coverage in a genome browser; it can also generate
  coverage views on any ranges specified, output as a
  coverageTable.csv file. coverageTable containing read
  maxPosition, maxHeight, coverageSums at specific genomic
  ranges. It also set readCoverage slot, and coverageView
  slot if ranges are provided.

  Note, parameter 'feature','ranges','annotationFile' or
  'slit.metaData' are used only with single
  SeqData/bamFile.
}
\note{
  Parameter 'feature','ranges','annotationFile' or
  'slit.metaData' are used only with single
  SeqData/bamFile.
}
\section{Usage}{
  { viewCoverage(obj=NULL,obj.control=NULL,

  bamFile=character(0),bamFile.control=character(0),weight.control=1,

  feature=character(0),ranges=GRanges(),annotationFile=character(0),

  split.metaData=F,bigWig=F,description=character(0))}
}
\examples{
# viewCoverage(obj)
# viewCoverage(obj,ranges=GRanges)
# viewCoverage(obj,annotationFile=annotationFile)

# viewCoverage(bamFile=bamFile)
# viewCoverage(bamFile=bamFile,ranges=GRanges)
# viewCoverage(bamFile=bamFile,annotationFile=annotationFile)

# viewCoverage(obj,obj.control)
# viewCoverage(bamFile,bamFile.control)
# to view coverageView slot
# coverageView(obj)[[1]]
# viewCoverage(obj,feature="transposon",split.metaData=T,description="obj_transposon")
}

