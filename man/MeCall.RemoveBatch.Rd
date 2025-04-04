% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeCall.RemoveBatch.R
\name{MeCall.RemoveBatch}
\alias{MeCall.RemoveBatch}
\title{Removing identified batch effects.}
\usage{
MeCall.RemoveBatch(
  meth,
  pd,
  batches = NULL,
  interest = "Sample_Group",
  do.sva = FALSE,
  exception = NULL
)
}
\arguments{
\item{meth}{A matrix of methylation level.}

\item{pd}{A dataframe of sample data file is contained in an object returned by \code{MeCall.ReadIdat()} or \code{MeCall.Filtering()}.}

\item{batches}{vector. The column names present in the sample data file. If set to NULL, the correlation coefficient will be calculated for all column names present in the sample data file.}

\item{interest}{The column names in the sample data file containing the phenotype of interest.}

\item{do.sva}{Logical. Perform surrogate variable analysis. If set to True, perform surrogate variable analysis (SVA) to identify and remove unknown batch effects.}

\item{exception}{It is a vector composed of column names from the sample data file. Surrogate variables that are correlated with the specified variables are excluded from the adjustment. If set to NULL, only the values of the interest variable are considered.}
}
\value{
A list object that includes a methylation level matrix with the
specified batch effect removed and a sample metadata data frame modified
according to the study design.
\itemize{methylation.level} {A methylation level matrix with batch effects removed}
\itemize{modified.pd} {Adjusted sample metadata according to the study design}
\itemize{groups} {The original names of each group assigned to the modified group}
}
\description{
The function removes identified known batch effects through \code{MeCall.FindBatch()},
or detects and removes unknown batch effects using the singular value decomposition (SVD).

Known batch effects can be effectively removed using the ComBat algorithm.
Users can specify the variable names to be removed for batch correction by
inputting them into the \code{batch} argument.

Sometimes, within the components of a known batch, there are instances
where one component has only one sample associated with it.

In such cases, the ComBat algorithm can not calculate variance properly.

Therefore, \code{MethylCallR} sets the \code{mean.only} parameter to TRUE when
performing the \code{ComBat} algorithm.

For more details, please refer to the \code{ComBat()} function in the \code{SVA} package.

When information on technical batch effects and potential confounding
factors is limited, utilizing surrogate variable analysis (SVA) to
identify latent variables can be beneficial.

Since major confounding factors and variables of interest are not
independent, considering all latent variables to adjust data may distort
the variance associated with the variables of interest.

To prevent this, \code{MethylCallR} excludes latent variables that are
significantly correlated with variables added to the \code{exception}
parameter from the correction process.
}
\examples{
\dontrun{
# Run [MeCall.RemoveBatch] with default setting.
data.rmBatch <- MeCall.RemoveBatch(meth = data.Norm$beta, pd = data.filtered$pd, 
batches=c("Slide","Array"), interest="Sample_Group", do.sva = FALSE)

# Run [MeCall.RemoveBatch] with SVA.
data.rmBatch.SVA <- MeCall.RemoveBatch(meth = data.Norm$beta, pd = data.filtered$pd, 
batches=c("Slide","Array"), interest="Sample_Group", do.sva = TRUE, exception = c("Age","Sex"))
}

}
\references{
For SVA R package
Leek, J.T. et al. (2012). The sva package for removing batch effects and other unwanted
variation in high-throughput experiments. Bioinformatics.
/url{https://doi.org/10.1093/bioinformatics/bts034}

For ComBat algorithm
Johnson, W.E. et al. (2007). Adjusting batch effects in microarray expression data using
empirical Bayes methods. Biostatistics. /url{https://doi.org/10.1093/biostatistics/kxj037}

For SVA
Leek, J.T. and Storey, J.D. (2008). A general framework for multiple testing dependence.
Proceedings of the National Academy of Sciences. /url{https://doi.org/10.1073/pnas.0808709105}
}
\seealso{
\code{\link[=MeCall.FindBatch]{MeCall.FindBatch()}}, \code{\link[sva:ComBat]{sva::ComBat()}}, \code{\link[sva:sva]{sva::sva()}}
}
\author{
Hyun-Ho Yang
}
