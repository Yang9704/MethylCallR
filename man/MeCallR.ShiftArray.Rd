% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeCall.ShiftArray.R
\name{MeCallR.ShiftArray}
\alias{MeCallR.ShiftArray}
\title{Changing the version of the microarray.}
\usage{
MeCallR.ShiftArray(
  Object,
  From = "EPICv2",
  To = "EPICv1",
  type = c("rgSet", "mSet", "Matrix", "Vector")
)
}
\arguments{
\item{Object}{The object for which the version will be changed.}

\item{From}{The version of the inputted object.}

\item{To}{The version of the output object.}

\item{type}{The type of the inputted object (\code{rgSet},\code{mSet},\code{Matrix},\code{Vector}).}
}
\value{
An object with the changed version that retains the same format as the input object.
}
\description{
This function utilizes the Chain file generated by MethylCallR to change
the version of the microarray for the input object.

Specifically, this chain file is designed to retain only features
that are common across all array versions. For example, if probe X was
a type 1 probe in EPICv1, and both type 1 and type 2 probes coexist
in EPICv2, only the information related to the type 1 probe will
be preserved. This is essential because it minimizes potential bias
that can arise when using algorithms that consider differences in
probe type or channel color.

The \code{MeCall.ShiftArray()} function supports inputs in the form of
\code{minfi} objects (\code{RGChannelSet} and \code{MethylSet}), matrix, and vector.
After changing the version, it returns an object in the same form as
the inputted object.
}
\examples{
\dontrun{
# Change array version from EPICv2 to EPICv1
Shifted.rgSet <- MeCallR.ShiftArray(Object = data.Import$minfi.Set$rgSet, From = "EPICv2", 
To = "EPICv1", type = c("rgSet"))
}

}
\author{
Hyun-Ho Yang
}
