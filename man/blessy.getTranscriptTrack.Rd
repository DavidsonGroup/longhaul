% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blessy-getTranscriptTrack.R
\name{blessy.getTranscriptTrack}
\alias{blessy.getTranscriptTrack}
\title{Fetch and Convert UCSC Transcript-type Track to Data Frame}
\usage{
blessy.getTranscriptTrack(genome, track)
}
\arguments{
\item{genome}{A string specifying the genome/assembly identifier (e.g., "hg38").}

\item{track}{A string specifying the table identifier of the annotation track (e.g., "wgEncodeGencodeBasicV44").}
}
\value{
A data frame in BED format containing the transcript annotation track.
}
\description{
This function fetches a transcript annotation track from the UCSC Genome Browser using UCSC.utils and
converts it into a BED-formatted data frame in R, notably the NCBI RefSeq or GENCODE tracks.
}
\examples{
# Example: Fetch the "wgEncodeGencodeBasicV44" annotation transcript annotation track
bed_df <- blessy.getTranscriptTrack("hg38", "wgEncodeGencodeBasicV44")

}
