#' Generate model matrix for differential methylation analysis.
#'
#' @description
#' This function generates a model matrix containing information about 
#' covariates to be considered for conducting a differential methylation 
#' analysis based on the provided data.
#'
#' If the `adjustment` parameter is set to `TRUE`, the function will 
#' perform linear regression on the methylation level matrix to correct 
#' covariates instead of providing a model matrix. This function will 
#' return a list object containing the adjusted methylation level matrix 
#' and a data frame that includes covariate information.
#'
#' @param pd A data frame of modified sample data file from `MeCall.RemoveBatch()`. 
#' @param covariate Vector. The covariates identified by `MeCall.FindBatch()`.
#' @param interest The column names in the sample data file containing the phenotype of interest.
#' @param PCA.result A returned object from `MeCall.PCA()`. Set to `NULL` if you do not want to include it as covariate.
#' @param n.pc The number of principal components to include in the covariate.
#' @param Cell.comp A dataframe containing the estimated cell composition for each sample. Set to `NULL` if you do not want to include it as covariate.
#' @param celltypes A vector containing the names of the cell types to be included in the covariate from the given `Cell.comp` parameter.
#' @param adjustment Logical. Methylation levels are adjusted using linear regression for given covariates.
#' @param meth When `adjustment` is set to `TRUE`, the matrix of methylation levels to be adjusted.
#'
#' @return The design matrix including the given covariates. 
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [MeCall.RemoveBatch()], [MeCall.FindBatch()], [MeCall.PCA()], [MeCall.CellComp()]
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.MakeModelMatrix] with default setting.
#' data.design <- MeCall.MakeModelMatrix(pd = data.rmBatch$modified.pd, 
#' covariate = c("Age","Sex"), interest = "Sample_Group", PCA.result = data.PCA, n.pc = 2, 
#' Cell.comp = data.CellComp$cell.composition, celltypes = colnames(data.CellComp$cell.composition),
#' adjustment=FALSE, meth=NULL)
#' }
#'
#' @export
MeCall.MakeModelMatrix <- function(pd = data.rmBatch$modified.pd, covariate = c("Age","Sex"), interest = "Sample_Group", PCA.result = data.PCA, n.pc = 2, Cell.comp = data.CellComp$cell.composition, celltypes = colnames(data.CellComp$cell.composition), adjustment=FALSE, meth=NULL){
message("\n[MeCall]-[notice] : Make design matrix with given input parameters.")
cov = FALSE
pca = FALSE
cell = FALSE
fcov <- c(covariate,interest)

rownames(pd) <- pd$Sample_Name

if(!is.null(pd)){
        if(all.equal(sort(intersect(colnames(pd),fcov)),sort(fcov)) == TRUE){
#               if(covariate == c()){
                if(is.null(covariate)){
                        message("\n[MeCall]-[notice] : You didn't set any covariates. ['covariate' parameter = c()]")}
                else {message("\n[MeCall]-[notice] : You select ",length(covariate)," covariate. -> ",paste(covariate,collapse=', '))}
                cov <- TRUE}
        else {stop("\n[MeCall]-!!ERROR!! : Your covariate parameter elements are not included in [myLoad$pd]. Check it again.")}
        } else { stop("\n[MeCall]-!!ERROR!! : You have to set sample information table. ['pd' parameter = NULL]")}
if(!is.null(PCA.result)){
        if(n.pc > ncol(PCA.result@Components)){stop("\n[MeCall]-!!ERROR!! : Your 'n.pc' parameter has too high value for PCA.result input data. your PCA data has '",ncol(PCA.result@Components),"' PCs")}
        else { pca <- TRUE}
        } else { message("\n[MeCall]-[notice] : You didn't set PCA result. ['PCA.result' parameter = NULL]")}
if(!is.null(Cell.comp)){
        if(all.equal(intersect(colnames(Cell.comp),celltypes),celltypes) == TRUE){cell <- TRUE}
        else { stop("\n[MeCall]-!!ERROR!! : Your celltypes parameter elements are not included in [Cell.comp]. Check it again.")}
        } else {message("\n[MeCall]-[notice] : You didn't set Cell composition result. ['Cell.comp' parameter = NULL]")}
fpd <- pd
fPCA.result <- PCA.result
fCell.comp <- Cell.comp
message("\n[MeCall]-[notice] : Pickup covariate from input data.")
if(cov){temcov <- paste(covariate,collapse='+')
pdata <- fpd
} else {temcov=NULL
pdata <- NULL}
if(pca){tempc <- paste(colnames(fPCA.result@Components)[1:n.pc],collapse='+')
pcadata <- fPCA.result
} else { tempc=NULL
pcadata <- NULL}
if(cell){temcell <- paste(celltypes,collapse='+')
celldata <- fCell.comp
} else { temcell=NULL
celldata <- NULL}
message("\n[MeCall]-[notice] : Finish check input data.")
message("\n[MeCall]-[notice] : Merging input data sets...")
if(is.null(pcadata)){
A <- pdata
rownames(A) <- rownames(pdata)}
else if (is.null(pdata)){
        A <- pcadata
rownames(A) <- rownames(pcadata)}
else if(nrow(pdata) == nrow(pcadata)){
        if(all.equal(pdata$Sample_Name,rownames(pcadata))){
                   A <- cbind(pdata,pcadata)}
        else {stop("\n[MeCall]-!!ERROR!! : pd data and PCA result data have different order. you should fix it.")}
rownames(A) <- rownames(pdata)
message("\n[MeCall]-[notice] : finish pd + pca")}
else {stop("\n[MeCall]-!!ERROR!! : row number of your data is not matched. You should recheck it again.")}
if(is.null(A)){
        whole.data <- celldata
rownames(whole.data) <- rownames(celldata)}
else if (is.null(celldata)){
        whole.data <- A
rownames(whole.data) <- rownames(A)}
else if(nrow(A) == nrow(celldata)){
        if(all.equal(rownames(A),rownames(celldata))){
        whole.data <- cbind(A,celldata)}
        else {stop("\n[MeCall]-!!ERROR!! : Your cell composition data have different sample name order. you should fix it.")}
rownames(whole.data) <- rownames(A)
message("\n[MeCall]-[notice] : finish all.")}
else {stop("\n[MeCall]-!!ERROR!! : row number of your data is not matched. You should recheck it again.")}
if(is.null(whole.data)){stop("\n[MeCall]-!!ERROR!! : You didn't set any covariate. You should input data to make design matrix.")}
if(adjustment==TRUE){
        if(!is.null(meth)){
                message("\n[MeCall]-[notice] : Adjustment = TRUE, methylation level table = TRUE, Adjustment start with linear regression method.")
                message("\n[MeCall]-[notice] : Adjusted methylation level table and combined covariate table will be returned.")
                b <- as.matrix(meth)
                whole.data <- as.data.frame(whole.data)
                mtype <- ismethyltype(b)
                message("\n[MeCall]-[notice] : ",mtype,"-value is detected.")
                adj.meth <- R4A(b, whole.data, c(temcell,temcov,tempc),mtype)
                adj <- list(adjusted.meth=adj.meth, cov.table = whole.data)
                return(adj)}
                else {stop("\n[MeCall]-!!ERROR!! : There is no methylation level table. Check it again.")}
        } else {
message("\n[MeCall]-[notice] : Generating model matrix.")
frm <- as.formula(paste(paste(paste('~', paste(c(temcell,temcov,tempc), collapse='+')),'+'),interest))
design <- model.matrix(frm,data=whole.data)
message("\n[MeCall]-[notice] : model matrix will be return in your variable.")
return(design)}
}

