#' Finding potential confounding factors.
#'
#' @import ggplot2
#' @import stringr
#' @import stats
#' @importFrom reshape melt 
#' @importFrom scales squish
#' 
#' @description
#' Perform principal component analysis (PCA) on the input data and calculate
#' the correlations between the principal components and the variables 
#' specified in the batch parameter to assess the association between these
#' variables and the input data.
#'
#' @param meth A matrix of methylation level.
#' @param pd A dataframe of sample data file is contained in an object returned by `MeCall.ReadIdat()` or `MeCall.Filtering()`.
#' @param batches vector. The column names present in the sample data file. If set to NULL, the correlation coefficient will be calculated for all column names present in the sample data file.
#' @param n.comp Correlation coefficients are calculated using a given number of principal components. If set to 0, a maximum of 15 principal components are used.
#'
#' @return a ggplot2 object representing a heatmap that illustrates the 
#' correlation between principal components of methylation level and variables 
#' in sample data file. The correlation coefficients are displayed at the 
#' center of each tile. Blue tiles indicate significant correlations 
#' (p < 0.05), with darker shades of blue representing lower p-values.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [MeCall.ReadIdat()], [MeCall.Filtering()]
#'
#' @examples
#' \dontrun{
#' data.Heatmap <- MeCall.FindBatch(meth = data.Norm, pd = data.filtered$pd, batches=NULL, 
#' n.comp = 0)
#' 
#' # To call the ggplot2 object
#' data.Heatmap
#'
#' # To save the ggplot2 object as pdf file
#' pdf("./your_directory/data_Heatmap.pdf")
#' data.Heatmap
#' dev.off()
#' }
#'
#' @export
MeCall.FindBatch <- function(meth, pd, batches = NULL, n.comp = 0){
message("\n[MeCall]-[notice] : Finding significant batch effect in your data set.")
message("\n[MeCall]-[notice] : Please check whether the technical batch effect is factor or character type.")

if(is.null(batches)){
message("\n[MeCall]-[notice] : Estimate the batch effect for all column names existing in the pd data.")
batches <- colnames(pd)}

if(length(setdiff(batches, colnames(pd))) >0){
stop("[MeCall]-!!ERROR!! : There is no column name in pd data : ",paste(setdiff(batches, colnames(pd)),collapse = ","))}


batch.df <- pd[,batches]
batches.c <- c()
for (i in batches){if(is.character(batch.df[,i]) | is.factor(batch.df[,i])){batches.c <- c(batches.c,i)}}
CC <- lapply(batches.c, function(x){as.numeric(as.factor(batch.df[,x]))})
CC <- as.data.frame(do.call("cbind", CC))
colnames(CC) <- batches.c

for (i in batches.c){batch.df[,i] <- CC[,i]}

Pcomp <- prcomp(t(meth), scale=TRUE)
message("\n[MeCall]-[notice] : n.comp parameter is not specified. The Best number of component based on 'elbow' method will be used for estimation of correlation.")
if(n.comp == 0){n.comp <- pred.n.comp(Pcomp)}
if(n.comp > 15){
message("\n[MeCall]-[notice] : [",n.comp," PCs] are too high. The heatmap can display a maximum of 15 PCs.")
n.comp <- 15
}
message("\n[MeCall]-[notice] : [",n.comp," PCs] are selected.")
message("\n[MeCall]-[notice] : Generate heatmap.")
prcomp.mat <- Pcomp$x[,1:n.comp]

cor.pval.mat <- lapply(batches, function(a){apply(prcomp.mat,2,function(b){cor.test(b,batch.df[,a],method="spearman",exact=FALSE)$p.value})})
cor.cor.mat <- lapply(batches, function(a){apply(prcomp.mat,2,function(b){cor.test(b,batch.df[,a],method="spearman",exact=FALSE)$estimate})})
names(cor.pval.mat) <- batches
names(cor.cor.mat) <- batches
cor.pval.mat <- as.data.frame(do.call("cbind",cor.pval.mat))
cor.cor.mat <- as.data.frame(do.call("cbind",cor.cor.mat))

#library(reshape)
melt.pval <- reshape::melt(cor.pval.mat)
melt.pval <- setNames(melt.pval, c("Batch", "pval"))
melt.pval$pc <- rownames(cor.pval.mat)
melt.cor <- reshape::melt(cor.cor.mat)
melt.pval$cor <- round(melt.cor$value,2)

bi.pval <- apply(apply(cor.pval.mat, 2 , function(i){ifelse(i < 0.05, 1,0)}),2,sum)
need.adj <- names(bi.pval[bi.pval>0])
ord <- stringr::str_sort(unique(melt.pval$pc), numeric = TRUE)

heatmap <- ggplot2::ggplot(melt.pval, aes(x = factor(pc, levels = ord), y = Batch, fill = pval)) + geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = cor), color = "black", size = 4) + scale_fill_gradient(low = "#004C99", high = "white", guide = "none", limits = c(0, 0.05),oob=scales::squish) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.title = element_blank())

message("\n[MeCall]-[notice]: White tile represents not significant correlation between batch and PC. Darker blue color represents lower p-value of correlation.")
message("\n[MeCall]-[notice] : Follow elements have significant correlation with principal components of your data : ",paste(need.adj,collapse = ', '),". \n                    You may use 'ComBat' algorithm to remove technical batch effect (MeCall.RemoveBatch) or generate 'modelmatrix' to adjust covariates (MeCall.MakeModelMatrix) before differetial methylation analysis.")

return(heatmap)
}


