% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{methods}
\name{summarizeBases}
\alias{summarizeBases}
\title{summarizeBases}
\usage{
summarizeBases(obj=NULL,bamFile=character(0),annotationFile=character(0),baseReport=charcter(0),description=character(0))
}
\arguments{
\item{obj}{A SeqData object.}

\item{bamFile}{The full path to the sequencing data. The data needs to be in
bam file format.}

\item{annotationFile}{csv annotation file for specifying ranges of interest.
  It must includes the miminum five field: 'chrom','txSart','txEnd','strand'
  and 'name'.}

\item{baseReport}{DNA methylaton report from methylaton aligners such as
  bismark, in the form of a csv file has minimu four columns "chromosome",
  "position", "strand", "count_methylated", "count_unmethylated". if strand information is not available, fill it with "*".}

\item{description}{Optional. User can input short desciption for the output
  file, it will show up in the file name of the output file.}
}
\value{
\itemize{
\item{baseTable.csv} a csv file contains the base summarization.
}
}
\description{
count base measurement (e.g methylation) on any genomic regions.
}
\details{
the method can be applied to any parameters measured along genome.
  essentially, the function is to group the items in bases based on region
  provided, then comput whatever statistics over those groups. Internally,
  it calls getBaseAlignment.
}
\section{Usage}{
{summarizeBases(obj=NULL,bamFile=character(0),annotationFile=character(0),baseReport=charcter(0),description=character(0))}
}
\examples{
# summarizeBases(obj)
# summarizeBases(baseReport)
# summarizeBases(bamFile=bamFile,annotationFile=annotationFile)
}

