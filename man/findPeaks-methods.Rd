\docType{methods}
\name{findPeaks}
\alias{findPeaks}
\title{findPeaks}
\arguments{
  \item{bamFie.chip}{Full path to chip bam file.}

  \item{bamFile.control}{Full path to control bam file.}

  \item{obj.chip}{SeqData object containing information on
  chip.}

  \item{obj.control}{SeqData object containing information
  on control.}

  \item{description}{Optional. User can input short
  desciption for the output file, it will show up in the
  file name of the output file.}
}
\value{
  \itemize{ \item{peakTable.csv} A csv file contaning
  information about the significantly enriched regions,
  including
  chromosome,start,end,peakPositions,peakCoverageSums,rpkm.peakCoverageSums,chipPeakHeights,controlPeakHeights,p.values,fdr,foldChange.
  \item{coverageView slot} A RleViews object containing
  Views of the peak regions. }
}
\description{
  Find reads enriched regions using poission distribution
  model.
}
\section{Usage }{
  {
  findPeaks(obj.chip=NULL,obj.control=NULL,bamFile.chip=character(0),bamFile.control=character(0),fdr.max=1e-5,foldChange.min=2)
  }
}
\examples{
#findPeaks(bamFile.chip)
#findPeaks(bamFile.chip,bamFile.control)
#findPeaks(obj.chip)
#findPeaks(obj.chip,obj.control)
#findPeaks(obj.chip,obj.control,bamFile.chip,bamFile.control)
# to view coverageView slot
# coverageView(obj.chip)[[1]]
}

