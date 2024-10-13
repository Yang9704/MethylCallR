#' Removing identified batch effects.
#'
#' @import sva
#' 
#' @description 
#' The function removes identified known batch effects through `MeCall.FindBatch()`,
#' or detects and removes unknown batch effects using the singular value decomposition (SVD).
#' 
#' Known batch effects can be effectively removed using the ComBat algorithm.
#' Users can specify the variable names to be removed for batch correction by 
#' inputting them into the `batch` argument.
#'
#' Sometimes, within the components of a known batch, there are instances 
#' where one component has only one sample associated with it.
#'
#' In such cases, the ComBat algorithm can not calculate variance properly.
#'
#' Therefore, `MethylCallR` sets the `mean.only` parameter to TRUE when 
#' performing the `ComBat` algorithm.
#'
#' For more details, please refer to the `ComBat()` function in the `SVA` package.
#'
#'
#' When information on technical batch effects and potential confounding 
#' factors is limited, utilizing surrogate variable analysis (SVA) to 
#' identify latent variables can be beneficial.
#'
#' Since major confounding factors and variables of interest are not 
#' independent, considering all latent variables to adjust data may distort 
#' the variance associated with the variables of interest.
#'
#' To prevent this, `MethylCallR` excludes latent variables that are
#' significantly correlated with variables added to the `exception` 
#' parameter from the correction process.
#'
#'
#' @param meth A matrix of methylation level.
#' @param pd A dataframe of sample data file is contained in an object returned by `MeCall.ReadIdat()` or `MeCall.Filtering()`.
#' @param batches vector. The column names present in the sample data file. If set to NULL, the correlation coefficient will be calculated for all column names present in the sample data file.
#' @param interest The column names in the sample data file containing the phenotype of interest.
#' @param do.sva Logical. Perform surrogate variable analysis. If set to True, perform surrogate variable analysis (SVA) to identify and remove unknown batch effects.
#' @param exception It is a vector composed of column names from the sample data file. Surrogate variables that are correlated with the specified variables are excluded from the adjustment. If set to NULL, only the values of the interest variable are considered.
#'
#' @return A list object that includes a methylation level matrix with the 
#' specified batch effect removed and a sample metadata data frame modified 
#' according to the study design.
#' \itemize{methylation.level} {A methylation level matrix with batch effects removed}
#' \itemize{modified.pd} {Adjusted sample metadata according to the study design}
#' \itemize{groups} {The original names of each group assigned to the modified group}
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [MeCall.FindBatch()], [sva::ComBat()], [sva::sva()]
#'
#' @references
#' For SVA R package 
#' Leek, J.T. et al. (2012). The sva package for removing batch effects and other unwanted 
#' variation in high-throughput experiments. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/bts034}
#'
#' For ComBat algorithm
#' Johnson, W.E. et al. (2007). Adjusting batch effects in microarray expression data using 
#' empirical Bayes methods. Biostatistics. /url{https://doi.org/10.1093/biostatistics/kxj037}
#'
#' For SVA 
#' Leek, J.T. and Storey, J.D. (2008). A general framework for multiple testing dependence. 
#' Proceedings of the National Academy of Sciences. /url{https://doi.org/10.1073/pnas.0808709105}
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.RemoveBatch] with default setting.
#' data.rmBatch <- MeCall.RemoveBatch(meth = data.Norm$beta, pd = data.filtered$pd, 
#' batches=c("Slide","Array"), interest="Sample_Group", do.sva = FALSE)
#'
#' # Run [MeCall.RemoveBatch] with SVA.
#' data.rmBatch.SVA <- MeCall.RemoveBatch(meth = data.Norm$beta, pd = data.filtered$pd, 
#' batches=c("Slide","Array"), interest="Sample_Group", do.sva = TRUE, exception = c("Age","Sex"))
#' }
#'
#' @export
MeCall.RemoveBatch <- function(meth, pd, batches=NULL, interest="Sample_Group", do.sva= FALSE, exception = NULL){
message("\n[MeCall]-[notice] : Estimate batch effect and start remove procedure based on SVA R package.")

message("\n[MeCall]-[notice] : Known batch list : ",ifelse(!is.null(batches), paste(batches ,collapse=", "), "None"))
message("                 Surrogate variable analysis to find unknown batch effect : ", do.sva)


if(!is.null(batches)){
if(FALSE %in% (batches %in% colnames(pd))){stop("\n[MeCall]-!!ERROR!! : Variable in 'batches' parameter is not exist in 'pd' data. : ", paste(batches[!batches %in% colnames(pd)],collapse = ", "))
} else {
Check.v <- lapply(batches, function(x) {table(pd[,x])})
names(Check.v) <- batches
Check.v <- unlist(lapply(Check.v, function(x) {ifelse(min(x)<2,"fail","check")}))
Check.v.f <- names(Check.v[Check.v == "fail"])
if(length(Check.v.f)>=1){
message("\n[MeCall]-[notice] : ",length(Check.v.f)," variables have 1 sample in 1 phenotype. : ", paste(Check.v.f, collapse = ", "))
message("                 Those variables will be adjusted with reference method.")
batches.ref <- batches[batches %in% Check.v.f]
batches.nref <- batches[!batches %in% Check.v.f]
}else{batches.nref <- batches
batches.ref <-c()}
message("\n[MeCall]-[notice] : ",(length(batches.ref)+length(batches.nref))," variables are confirmed.")
}}

if(ismethyltype(meth) == "beta"){
if((sum(meth <= 0) > 0) | (sum(meth >= 1) > 0)){
message("\n[MeCall]-[notice] : The beta value matrix contains values of 0 or 1, which may interfere with the ComBat algorithm. These values will be replaced with 0.00001 and 0.99999, respectively.")
meth[meth <= 0] <- 0.00001
meth[meth >= 1] <- 0.99999
}
c_meth <- BetatoM(meth)
}else{c_meth <- meth}


fpd <- CheckMultipleGroup(pd=pd, interest = interest)
CaCo <- fpd$groups
original.interest <- fpd$original.pd[,interest]
mpd <- cbind(fpd$modified.pd,original.interest)

if(!is.null(batches)){
for (i in 1:length(batches)){
if((length(batches) - i) > 0){
t.batch <- batches[i]
message("\n[MeCall]-[notice] : Apply ComBat algorithm to adjust following batch effect : ",t.batch)
l.batch <- batches[(i+1) : length(batches)]
frm <- as.formula(paste0("~", paste(c(interest,l.batch), collapse = "+")))

} else if ((length(batches) - i) == 0){
t.batch <- batches[i]
message("\n[MeCall]-[notice] : Apply ComBat algorithm to adjust following batch effect : ",t.batch)
frm <- as.formula(paste0("~", interest))

} else{stop("[MeCall]-!!ERROR!! : Please check your batches parameter.")}

modcombat = model.matrix(frm, data=mpd)
if(t.batch %in% batches.nref){
c_meth <- sva::ComBat(dat = c_meth, batch = mpd[,t.batch], mod = modcombat,par.prior=TRUE, prior.plots = FALSE)
}else{
message("\n[MeCall]-[notice] : ",t.batch," has one sample in one phenotype. Variance adjustment will be skipped. -> mean.only = TRUE")
batch.table <- table(mpd[,t.batch])
refs <- names(batch.table[which(batch.table == max(batch.table))][1])
c_meth <- sva::ComBat(dat = c_meth, batch = mpd[,t.batch], mod = modcombat,par.prior=TRUE, prior.plots = FALSE,mean.only=TRUE,ref.batch = refs)}

colnames(c_meth) <- colnames(meth)
rownames(c_meth) <- rownames(meth)
}
}else{
message("\n[MeCall]-[notice] : There is no variable to adjust using ComBat algorithm.")
}




if(do.sva){
sva.pd <- pd
sva.pd[,interest] <- mpd[,interest]
if(sum(sva.exception %in% colnames(sva.pd)) != length(sva.exception)){
stop("\n[MeCall]-!!ERROR!! : There is an invalid column name. Please check and try again.")}
message("\n[MeCall]-[notice] : Surrogate variable analysis start.")
svaNullmod <- model.matrix(~ 1,data=sva.pd)
svamod <- model.matrix(as.formula(paste0("~", interest)),data=sva.pd)
message("\n[MeCall]-[notice] : Find appropriate number of latent factors based on leek method.")
n.sv <- num.sv(c_meth,svamod,method="leek")
}else{
if(ismethyltype(meth) == "beta"){
c_meth <- MtoBeta(c_meth)
}
adjustedobj <- list(methylation.level = c_meth, modified.pd = mpd, groups = CaCo)
return(adjustedobj)
}


if(n.sv == 0){
message("\n[MeCall]-[notice] : There is no significant latent factor. MethylCallR will skip sva procudure.")
if(ismethyltype(meth) == "beta"){
c_meth <- MtoBeta(c_meth)
}
adjustedobj <- list(methylation.level = c_meth, modified.pd = mpd, groups = CaCo)
return(adjustedobj)
}else {message("\n[MeCall]-[notice] : Number of significant latent factors :",n.sv)
message("\n[MeCall]-[notice] : Estimate the surrogate variables based on latent factors.")
svobj = sva(c_meth,svamod,svaNullmod,n.sv=n.sv)
svs <- sapply(1:n.sv, function(i){paste0("SV",i)})
colnames(svobj$sv) <- svs
rownames(svobj$sv) <- colnames(c_meth)
}

message("\n[MeCall]-[notice] : Estimate correlation batween ",n.sv," surrogate variables and known variables. SVs which have no correlation with [exception] and [interest] parameter will be used for unknown batch effect correction.")

if(!is.null(sva.exception)){
Exp.var <- c(sva.exception,interest)
} else {Exp.var <- c(interest)}

remain.var <- setdiff(colnames(sva.pd), Exp.var)

message("\n[MeCall]-[notice] : Following variables will not be considered : ", paste(remain.var, collapse = ", "))

Cvs <- sva.pd[,Exp.var]
nas <- colSums(is.na(Cvs))
if(any(nas > 0)){
Failed.var <- names(nas[nas>0])
message("\n[MeCall]-[WARNING] : NAs were detected in some variables. SVA does not allow NAs. The following variables will be excluded from the exception parameter : ",paste(Failed.var, collapse = ", "))
Cvs <- as.data.frame(Cvs[, names(nas[nas==0])])
}
colnames(Cvs)[ncol(Cvs)] <- interest
c_cols <- sapply(Cvs, is.character)
Cvs[,c_cols] <- lapply(Cvs[,c_cols], as.factor)
Cvs[,c_cols] <- lapply(Cvs[,c_cols], as.numeric)

message("\n[MeCall]-[notice] : Finding significant correlation.")
cor.pval.mat <- lapply(colnames(Cvs), function(a){apply(svobj$sv,2,function(b){cor.test(b,Cvs[,a],method="spearman",exact=FALSE)$p.value})})

names(cor.pval.mat) <- colnames(Cvs)
cor.pval.mat <- as.data.frame(do.call("cbind",cor.pval.mat))
bi.pval <- t(as.matrix(apply(cor.pval.mat, 2 , function(i){ifelse(i < 0.05, 1,0)})))
bi.pval <- rowSums(bi.pval)
need.adj <- rownames(cor.pval.mat[which(bi.pval == 0),])

if(length(need.adj) == 0){
message("\n[MeCall]-[notice] : There is no SVs which have no significant correlation with known variables. Adjustment procedure will be skipped.")
if(ismethyltype(meth) == "beta"){
c_meth <- MtoBeta(c_meth)
}
adjustedobj <- list(methylation.level = c_meth, modified.pd = mpd, groups = CaCo)
return(adjustedobj)
}else{
message("\n[MeCall]-[notice] : ",length(need.adj)," SVs have no significant correlation with variables included in sample data file : ", paste(Exp.var, collapse = ", "))
SVs <- as.data.frame(svobj$sv[,need.adj])
message("\n[MeCall]-[notice] : Estimating surrogate variables Done. ",ncol(SVs)," SV objects will be adjusted using linear regression.")
svs <- need.adj
colnames(SVs) <- svs
mtype <- ismethyltype(c_meth)
c_meth_df <- as.data.frame(c_meth)
c_meth <- R4A(c_meth_df, SVs, svs, type=mtype)
}

if(ismethyltype(meth) == "beta"){
c_meth <- MtoBeta(c_meth)
}

adjustedobj <- list(methylation.level = c_meth, modified.pd = mpd, groups = CaCo)
return(adjustedobj)}
