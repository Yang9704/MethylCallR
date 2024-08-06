#' Estimate cell composition of various tissue type.
#' 
#' @import minfi
#' @import deconvR
#' @importFrom FlowSorted.Blood.EPIC estimateCellCounts2
#' 
#' 
#' @description
#' This function involves estimating the proportions of cell types for each
#' sample based on a selected or provided reference matrix and correcting 
#' biases using linear regression that result from cell type heterogeneity 
#' between samples.
#' 
#' MethylCallR provides a function to estimate white blood cell composition
#' using the `FlowSorted.Blood.EPIC` R packages.
#' 
#' In case of solid tissue samples, estimation of the cell type proportions is 
#' conducted using reference matrices and functions implemented in 
#' the `EPISCORE` and `deconvR` R packages.
#' 
#' If a reference matrix is available, the user may want to set the `method`
#' parameter to `manual` and input the reference matrix into the `Reference`
#' parameter. Using the `deconvolute()` function implemented in the `deconvR` R
#' package, the cell type proportions will be calculated. Please ensure that
#' the order of probes in the reference matrix matches the order of probes 
#' in the methylation level matrix.
#' 
#' If the `adjust` parameter is set to `TRUE`, `MethylCallR` will automatically 
#' correct for estimated cell type heterogeneity using linear regression for
#' all inferred cell types. This function consecutively performs estimation 
#' and correction, and if `return.All` is also set to `TRUE`, the estimated 
#' cell type proportions and the corrected methylation levels will be 
#' returned together as part of the result object.
#'
#' @param data A return object from `MeCall.Filtering()` for blood tissue type.
#' @param meth A matrix of methylation level for solid tissue type or adjustment.
#' @param processMethod Option for blood tissue type. The normalization method to be applied along with the reference data ("preprocessNoob", "preprocessFunnorm", "preprocessQuantile", "preprocessSwan", "preprocessRaw").
#' @param CellType The tissue type of the sample.
#' @param method Option for solid tissue type. The reference matrix of the selected method will be applied.
#' @param Reference Data frame of reference to be used for cell type estimation. The columns are the sample names, and the rows are the CpG IDs. The `method` parameter must be set to `manual`.
#' @param adjust Logical. Methylation levels should be provided. Methylation levels are adjusted and returned based on the estimated cell composition.
#' @param return.All Logical. Both the adjusted methylation levels and a data frame of the estimated cell type composition are returned.
#' @param arraytype The version of the microarray.
#'
#' @return A data frame of estimated cell compostion or a matrix of adjusted methylation level. If `return.All` is set to TRUE, a list object containing two data frames will be returned.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [FlowSorted.Blood.EPIC::estimateCellCounts2()],  [deconvR::deconvolute()]
#'
#' @references 
#' Adjustment code comes from  
#' Jones, M.J. et al. (2017). Adjusting for Cell Type Composition in DNA Methylation Data Using a 
#' Regression-Based Approach. Methods in Molecular Biology (Clifton, N.J.). 
#' /url{https://doi.org/10.1007/7651_2015_262}
#'
#' For FlowSorted.Blood.EPIC
#' Salas, L.A. et al. (2018). An optimized library for reference-based deconvolution of 
#' whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. 
#' Genome Biology. /url{https://doi.org/10.1186/s13059-018-1448-7}
#' 
#' For EPISCORE
#' Teschendorff, A.E. et al. (2020). EPISCORE: cell type deconvolution of bulk tissue DNA 
#' methylomes from single-cell RNA-Seq data. Genome Biology. 
#' /url{https://doi.org/10.1186/s13059-020-02126-9}
#'
#' For deconvR
#' Moss, J. et al. (2018). Comprehensive human cell-type methylation atlas reveals origins of 
#' circulating cell-free DNA in health and disease. Nature Communications. 
#' /url{https://doi.org/10.1038/s41467-018-07466-6}
#'
#'
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.CellComp] with Blood tissue type.
#' data.CellComp <- MeCall.CellComp(data = data.filtered, meth = data.rmBatch$methylation.level, 
#' processMethod = "preprocessNoob", CellType = "Blood", adjust = TRUE, return.All = FALSE)
#'
#' # Run [MeCall.CellComp] with Solid tissue type.
#' data.CellComp <- MeCall.CellComp(meth = data.rmBatch$methylation.level, CellType = "Lung", 
#' method = c("EpiSCORE"), adjust = TRUE, return.All = FALSE,  arraytype = "EPICv1")
#'}
#'
#' @export
MeCall.CellComp <- function(data = Filtered, meth = beta, processMethod = "preprocessNoob", CellType = "Blood",method = c("EpiSCORE","deconvR","manual") ,Reference = NULL, adjust = TRUE, return.All = FALSE, arraytype = "EPICv1"){
message("\n[MeCall]-[notice] : Estimate Cell Count procedure start.")

Blood.types <- c("Blood", "CordBloodCombined", "CordBlood", "CordBloodNorway", "CordTissue And Blood", "DLPFC")


if(!is.null(meth)){
Ori.meth <- meth
mtype <- ismethyltype(Ori.meth)
}

# Tissue tpye - BLOOD
if(CellType %in% Blood.types){
CellType <- match.arg(CellType, Blood.types, several.ok = FALSE)
if(data$TAG == "EPICv1"){
message("\n[MeCall]-[notice] : Array type = Illumina EPICv1 array.")
ref <- "IlluminaHumanMethylationEPIC"
Input.obj <- data$minfi.Set$rgSet
CustomCpGs = NULL
}else if(data$TAG == "450K"){
message("\n[MeCall]-[notice] : Array type = Illumina 450K array.")
ref <- "IlluminaHumanMethylation450K"
Input.obj <- data$minfi.Set$rgSet
CustomCpGs = NULL
}else{
message("\n[MeCall]-[notice] : Array type = Illumina EPICv2 array.")
ref <- "IlluminaHumanMethylationEPIC"
rgSet <- data$minfi.Set$rgSet
Input.obj = MeCallR.ShiftArray(Object = rgSet, From = "EPICv2", To = "EPICv1", type = c("rgSet"))
mSet <- minfi::preprocessRaw(Input.obj)
Filtered.CpGs <- unlist(data$filtered.CpG.list)
Filtered.EPICv1 <- MeCallR.ShiftArray(Filtered.CpGs, From = "EPICv2", To = "EPICv1", type = c("Vector"))
CustomCpGs = setdiff(mSet@NAMES,Filtered.EPICv1)
}
cell.prop <- FlowSorted.Blood.EPIC::estimateCellCounts2(rgSet = Input.obj,compositeCellType=CellType,processMethod = processMethod,CustomCpGs=CustomCpGs,referencePlatform = ref,returnAll = FALSE)
cell.df <- as.data.frame(cell.prop$prop)
celltypes <- colnames(cell.df)
}


# Tissue tpye - SOLID
if(!CellType %in% Blood.types){



if(is.null(meth)){
stop("\n[MeCall]-!!ERROR!! : Methylation level matrix should be provided for solid tissue type.")
}

if(mtype == "M"){
meth <- MtoBeta(meth)
}

EpiSCORE.types <- c("Bladder", "Brain", "Breast","Colon", "Esophagus", "Heart","Kidney","Liver","Lung","Pancreas","Prostate","Skin","Olfactory epithelium")

deconvR.types <- c('Erythrocyte_progenitors', 'Adipocytes', 'Neuron', 'Hepatocytes', 'Lung', 'Pancreas', 'Vascular', 'Colon', 'Left atrium', 'Bladder', 'Breast', 'Larynx', 'Kidney', 'Prostate', 'Thyroid', 'Upper Gastrointestinal', 'Cervix')

if(arraytype == "450K"){
message("\n[MeCall]-[notice] : Array type = Illumina 450K array.")
arraytype = "450k"
}else if(arraytype == "EPIC"){
message("\n[MeCall]-[notice] : Array type = Illumina EPICv1 array.")
arraytype = "850k"
}else{
message("\n[MeCall]-[notice] : Array type = Illumina EPICv2 array.")
arraytype = "850k"
meth <- MeCallR.ShiftArray(meth, From = "EPICv2", To = "EPICv1", type= c("Matrix"))
}

#library(EpiSCORE)
if(method == "EpiSCORE"){
if(!requireNamespace("EpiSCORE", quietly = TRUE)){
stop("\n[MeCall]-!!ERROR!! : [EpiSCORE] R package is not founded in your R environment. Please load [EpiSCORE] R package to run this function using EpiSCORE method {library(EpiSCORE)}. You can download [EpiSCORE] R package at following github page : https://github.com/aet21/EpiSCORE")
}
if(!CellType %in% EpiSCORE.types){
EpiSCORE.msg <- paste(EpiSCORE.types,collapse=", ")
stop("\n[MeCall]-!!ERROR!! : ",CellType," is not supported in EpiSCORE R package. Please refer following cell types : ",EpiSCORE.msg)
}

if(CellType == "Esophagus"){
prefix <- "Eso"
}else if(CellType == "Olfactory epithelium"){
prefix <- "OE"
}else{
prefix <- CellType
}

Ref.obj <- paste0(prefix,"Ref")
data(list = Ref.obj)
Ref.m.obj <- sprintf("mref%s.m", prefix)
refDNAm.m <- get(Ref.m.obj)

message("\n[MeCall]-[notice] : Estimating cell type proportion using reference matrix implemented in EpiSCORE R package.")

avDNAm.m <- EpiSCORE::constAvBetaTSS(meth,type=arraytype)
cell.df <- as.data.frame(EpiSCORE::wRPC(avDNAm.m,refDNAm.m,useW=TRUE,wth=0.4,maxit=200)$estF)
celltypes <- colnames(cell.df)
}


#library(deconvR)
if(method == "deconvR"){
data("HumanCellTypeMethAtlas")
if(!CellType %in% deconvR.types){
deconvR.msg <- paste(deconvR.types,collapse=", ")
stop("\n[MeCall]-!!ERROR!! : ",CellType," is not supported in deconvR R package. Please refer following cell types : ",deconvR.msg)
}
message("\n[MeCall]-[notice] : Estimating cell type proportion using reference matrix implemented in deconvR R package.")
Bulk_Beta <- as.data.frame(meth)
IDs <- rownames(meth)
Bulk_Beta <- cbind(IDs, Bulk_Beta)

cell.df <- deconvR::deconvolute(reference = HumanCellTypeMethAtlas, bulk = Bulk_Beta)$proportions
celltypes <- colnames(cell.df)
}

#library(deconvR) - With manual reference
if(method == "manual"){
message("\n[MeCall]-[notice] : Estimating cell type proportion using given reference matrix.")
message("\n[MeCall]-[notice] : Deconvolution is conducted using [deconvolute] function in 'deconvR' R package.")

Bulk_Beta <- as.data.frame(meth)
IDs <- rownames(meth)
Bulk_Beta <- cbind(IDs, Bulk_Beta)

cell.df <- deconvR::deconvolute(reference = Reference, bulk = Bulk_Beta)$proportions
celltypes <- colnames(cell.df)
}}

cell_result <- cell.df


if(adjust & is.null(meth)){message("\n[MeCall]-[notice] : Methylation level data is not detected. Adjustment procedure will be skipped.")}

if(!is.null(meth)){
if(adjust){
message("\n[MeCall]-[notice] : 'adjust' parameter is TRUE. Methylation level data is detected. Type : ",mtype)
message("\n[MeCall]-[notice] : Adjust ",mtype," value with linear regression method.")

Ori.meth <- as.data.frame(Ori.meth)
a_meth <- R4A(Ori.meth, cell.df, celltypes, type=mtype)

message("\n[MeCall]-[notice] : Finish adjustment. Cell type composition adjusted ",mtype," value will be returned")

cell_result <- a_meth

if(return.All){
message("\n[MeCall]-[notice] : returnAll == TRUE. Cell composition table will be returned together as a list object.")
cell_result <- list(adjusted.meth = a_meth, cell.prop = cell.df)
}
}}

return(cell_result)}

