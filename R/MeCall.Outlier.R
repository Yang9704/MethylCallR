#' Finding outlier samples using Mahalanobis distance.
#' 
#' @import ggplot2
#' @import ggrepel
#' @import stats
#' @importFrom car ellipse
#' 
#' @description
#' The Mahalanobis distance is calculated independently for each subgroup of
#' samples based on the elements in the `group` variable.
#'
#' In each subgroup, the Mahalanobis distance is calculated `nrep` times 
#' using `n.component` principal components of the subset of `meth`.
#'
#' Since the Mahalanobis distance is sensitive to outliers, in each iteration,
#' a random sub-sample of `sample.fraction*100`% is selected to calculate 
#' the covariance. This covariance is then applied to the entire sample 
#' for Mahalanobis distance calculation.
#'
#' Samples with distance values outside the `conf` confidence interval of 
#' the chi-squared distribution with 2 degrees of freedom are designated 
#' as outliers. Samples identified as outliers in more than 20% of the 
#' total iterations are classified as "obvious outliers", while 
#' those identified in between 5% and 20% of the iterations are classified 
#' as "suggestive outliers". The detection rates are included in the 
#' returned result object.
#'
#' @param meth A matrix of methylation level.
#' @param pd A dataframe of sample data file.
#' @param group The samples are divided into groups based on the specified criteria, and outliers will be identified within each group.
#' @param nrep The number of iterations for outlier detection.
#' @param sample.fraction The proportion of subsamples for each iteration.
#' @param conf Samples with distance values beyond the given confidence interval are considered as outliers.
#' @param n.component The number of principal components to include in calculation of Mahalanobis distance.
#' @param plot Create a scatter plot using PC1 and PC2 to display outliers.
#'
#' @return The sample name of obvious outlier and suggestive outlier. (Optional) A ggplot2 object that displays outliers.
#'
#' @author Hyun-Ho Yang
#'
#' @examples
#' \dontrun{
#' # Identify outliers with default settings and generate scatter plot
#' Outliers <- MeCall.Outlier(meth = data.Norm$beta, pd = data.filtered$pd, group = "Sample_Group", 
#' nrep = 1000, sample.fraction = 0.9, conf = 0.999, n.component = 2, plot = TRUE)
#' 
#' # To call the ggplot2 object
#' Outliers$Plot
#'}
#'
#' @export
MeCall.Outlier <- function(meth = norm, pd = filtered$pd ,group = NULL ,nrep=1000, sample.fraction = 0.9, conf= 0.999, n.component =2 ,plot = FALSE){

if(is.null(meth)){stop("\n[MeCall]-!!ERROR!! : Data is not exist.")
}else if(!is.matrix(meth)){stop("\n[MeCall]-!!ERROR!! : Data format should be matrix. Please check your input data.")
}else{data <- meth}
  
if(!group %in% colnames(pd)){stop("\n[MeCall]-!!ERROR!! : There is no group variable in your pd parameter.")
}
G <- pd[,group]
names(G) <- rownames(pd)

if(!is.null(G)){
group.var <- unique(G)
message("\n[MeCall]-[notice] : Saperate data according to 'group' variable.")
message("                 ",length(group.var)," groups : [",paste0(group.var, collapse=', '),"].")
group.idx <- lapply(1:length(group.var), function(x){which(G == group.var[x])})
names(group.idx) <- group.var
s.data <- lapply(group.idx, function(x){x <- as.vector(x)
meth[,x]})
message("\n[MeCall]-[notice] : Run PCA. ",n.component," components will be returned.")
data <- lapply(s.data, function(i){t.i <- t(i)
p <- prcomp(t.i)
comp <- p$x[,1:n.component]
return(comp)})}

message("\n[MeCall]-[notice] : Estimate possible outlier using mahalanobis distance method for each group.")
outlier.probability <- lapply(data, function(i){EstimatePossibleOutlier(i, nrep, sample.fraction, conf)})

obvious.outlier <- lapply(outlier.probability, function(x){names(x[x > 20])})
suggestive.outlier <- lapply(outlier.probability, function(x){names(x[5<x & x <= 20])})

Mahalanobis.Dist <- lapply(1:length(group.var), function(x){
Temp_data <- data[[x]]
outlierinfo <- obvious.outlier[[x]]
suggestiveinfo <- suggestive.outlier[[x]]
EstimateDist(Temp_data,outlierinfo,suggestiveinfo,conf)})
names(Mahalanobis.Dist) <- group.var

if(plot){message("\n[MeCall]-[notice] : Generate plot using 2 principal component.")
if(n.component >2){message("\n[MeCall]-[notice] : This plot use only 2 principal component. Mahalanobis distance ellipse might be different with expactation.")}
M.plot=DrawMahalanobisPlot(Mahalanobis.Dist)
}else{M.plot=NULL}

result = list(Obvious.outlier = unlist(obvious.outlier), Suggestive.outlier = unlist(suggestive.outlier), Outlier.Probability = unlist(outlier.probability),Plot = M.plot)

return(result)}

