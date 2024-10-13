#' Normalization of methylation level
#'
#' @Import ChAMP
#' @Import minfi
#' @ImportFrom ENmix preprocessENmix
#' @ImportFrom wateRmelon adjustedFunnorm adjustedDasen dasen
#'
#'
#' @description
#'
#' This function provide several method to normalize the data after 
#' removing samples and probes of poor quality.
#'
#' `MeCall.Norm()` function provides 6 normalization algorithms available 
#' across 4 different R packages (ENmix, minfi, WateRmelon, and ChAMP) :
#' * Functional normalization (`method = "Funn"`)
#' * Beta-mixture quantile normalization (BMIQ) (`method = "BMIQ"`)
#' * Subset-quantile within array normalization (SWAN) (`method = "SWAN"`)
#' * dasen (`method = "dasen"`)
#' * Exponential-Normal mixture signal intensity background correction (ENmix) (`method = "ENmix"`)
#' * Normal-exponential out-of-band normalization (Noob) (`method = "Noob"`)
#'
#' Due to variations in input formats required by each normalization algorithm,
#' we recommend to utilize the return object of `MethylCallR.filtering()` function.
#' 
#' The interpolatedXY step is an additional process that normalizes CpG sites on sex chromosomes without bias based on independently normalized autosomal CpG sites. This step can be utilized when the method parameter is set to `Funn` or `dasen`. For more details, refer to [wateRmelon::adjustedFunnorm()] and [wateRmelon::adjustedDasen()].
#'
#' @param data A return object from `MeCall.filtering()`.
#' @param method Normalization method ("BMIQ","Funn","SWAN","ENmix","dasen","Noob").
#' @param InterpolatedXY logical. If set to TRUE, the interpolatedXY method is activated, which normalizes CpG sites on sex chromosomes without bias based on independently normalized autosomal CpG sites. It is recommended to set this to TRUE if CpG sites on sex chromosomes are included in the methylation level data. This parameter applies only when the normalization method is selected as `Funn` or `dasen`.
#' @param offset Integer. A fixed value will be added when generating methylation levels. An integer value used identically to the value used in the `MeCall.ReadIdat()` function.
#' @param return.type The representation method of methylation levels (Beta-value = `B`, M-value = `M`).
#'
#' @return A matrix consisting of normalized methylation levels.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [minfi::preprocessFunnorm()], [minfi::preprocessNoob()], [minfi::preprocessSWAN()], 
#' [ChAMP::champ.norm()], [ENmix::preprocessENmix()], [wateRmelon::dasen()], [MeCall.filtering()], 
#' [wateRmelon::adjustedFunnorm()], [wateRmelon::adjustedDasen()]
#'
#' @references 
#' For Functional normalization
#' Triche, T.J., Jr et al. (2013). Low-level processing of Illumina Infinium DNA Methylation 
#' BeadArrays. Nucleic Acids Research. /url{ttps://doi.org/10.1093/nar/gkt090}
#'
#' For SWAN
#' Maksimovic, J. et al. (2012). SWAN: Subset-quantile Within Array Normalization for Illumina 
#' Infinium HumanMethylation450 BeadChips. Genome Biology. /url{https://doi.org/10.1186/gb-2012-13-6-r44}
#'
#' For Noob
#' Triche, T.J., Jr et al. (2013). Low-level processing of Illumina Infinium DNA Methylation 
#' BeadArrays. Nucleic Acids Research. /url{https://doi.org/10.1093/nar/gkt090}
#' 
#' For ENmix
#' Xu, Z. et al. (2016). ENmix: a novel background correction method for Illumina 
#' HumanMethylation450 BeadChip. Nucleic Acids Research. /url{https://doi.org/10.1093/nar/gkv907}
#' 
#' For BMIQ
#' Teschendorff, A.E. et al. (2013). A beta-mixture quantile normalization method for correcting 
#' probe design bias in Illumina Infinium 450 k DNA methylation data. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/bts680}
#'
#' For dasen
#' Pidsley, R. et al. (2013). A data-driven approach to preprocessing Illumina 450K methylation 
#' array data. BMC Genomics. /url{https://doi.org/10.1186/1471-2164-14-293}
#'
#' For interpolatedXY
#' Wang, Y. et al. (2022). InterpolatedXY: a two-step strategy to normalize DNA methylation 
#' microarray data avoiding sex bias. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/btac436}
#'
#' @examples
#' \dontrun{
#' # For [functional normalization] and [Beta-value]
#' data.Norm.FUNN <- MeCall.Norm(data = data.filtered, method = c("Funn"), offset=100) 
#'
#' # For [BMIQ] and [M-value]
#' data.Norm.BMIQ <- MeCall.Norm(data = data.filtered, method = c("BMIQ"), offset=100)
#' }
#'
#' @export
MeCall.Norm <- function(data = filtered, method=c("BMIQ","Funn","SWAN","ENmix","dasen","Noob"),InterpolatedXY = FALSE, offset = 100){

Normlist=c("BMIQ","Funn","SWAN","ENmix","dasen","Noob")

mod.data <- data

if(!all(method %in% Normlist)){
stop("\n[MeCall]-!!ERROR!! : Please check your method parameter. You may select the method whithin BMIQ, Funn, SWAN, ENmix, dasen, and Noob")}


if(length(method) >= 2){
m.ck <- method[-1]
rgs <- c("Funn","SWAN","ENmix","Noob")
if(any(m.ck %in% rgs)){
stop("\n[MeCall]-!!ERROR!! : Due to the input format (RGchannel-Set), Funn, SWAN, ENmix, and Noob should be used in first.")
}}


for (i in c(1:length(method))){

pr.method <- method[i]

if(pr.method == "BMIQ"){
mod.data <- BMIQ.f(mod.data)
}

if(pr.method == "Funn"){
mod.data <- Funn.f(mod.data,InterpolatedXY,offset)
}

if(pr.method == "SWAN"){
mod.data <- SWAN.f(mod.data,offset)
}

if(pr.method == "ENmix"){
mod.data <- ENmix.f(mod.data,offset)
}

if(pr.method == "dasen"){
mod.data <- dasen.f(mod.data,InterpolatedXY,offset)
}

if(pr.method == "Noob"){
mod.data <- Noob.f(mod.data,offset)
}
}

mod.data$beta[mod.data$beta <= 0] <- 0.00001
mod.data$beta[mod.data$beta >= 1] <- 0.99999

return(mod.data)}


# [BMIQ]------------------------------------
BMIQ.f <- function(data){
if(is.null(data$beta)){
stop("\n[MeCall]-!!ERROR!! : Beta-value matrix is not founded.")
}
message("\n[MeCall]-[notice] : [BMIQ] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) will be returned.")

if(data$TAG == "EPICv2"){
beta <- data$beta
beta[beta <= 0] <- 0.0001
beta[beta >= 1] <- 0.9999
mani <- callmanifest(data$TAG)
design.v <- subset(mani, rownames(mani) %in% rownames(beta))[rownames(beta),]$Infinium_Design_Type
names(design.v) <- rownames(beta)

suppressMessages(beta.n <- lapply(c(1:ncol(beta)), function(x){n.beta.v <- gmqn::BMIQ(beta[,x], design.v)}))
beta.n <- as.matrix(do.call(cbind,beta.n))


}else{
suppressMessages(library(ChAMPdata))
suppressMessages(beta.n <- ChAMP::champ.norm(beta=data$beta, method="BMIQ", arraytype=data$TAG))
}

colnames(beta.n) <- colnames(data$beta)
rownames(beta.n) <- rownames(data$beta)
beta.val <- beta.n
M.val <- BetatoM(beta.n)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= NULL,filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)

return(OUTobj)}

# [FUNN]------------------------------------
Funn.f <- function(data,InterpolatedXY,offset){
if(is.null(data$minfi.Set$rgSet)){
stop("\n[MeCall]-!!ERROR!! : RGchannel-Set is not founded.")
}
filtered.CpG <- unique(unlist(data$filtered.CpG.list))

message("\n[MeCall]-[notice] : [Functional normalization] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) and mSet will be returned.")

rgSet <- data$minfi.Set$rgSet
if(InterpolatedXY){
message("\n[MeCall]-[notice] : [interpolatedXY] method is activated. This function utilizes the [adjustedFunnorm] function included in wateRmelon R package. For more details, please refer to the following paper: [Yucheng Wang, Bioinformatics, 2022].")
normed <- wateRmelon::adjustedFunnorm(rgSet,nPCs=2, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = FALSE,verbose = TRUE)
} else {normed <- preprocessFunnorm(rgSet,nPCs=2, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = FALSE,verbose = TRUE)}

mset.n <- f.out.b.norm(normed,filtered.CpG)
beta.val <- toBeta(mset.n,offset)
M.val <- toM(mset.n,offset)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= list(rgSet = NULL, mSet = mset.n),filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)

return(OUTobj)}

# [SWAN]------------------------------------
SWAN.f <- function(data,offset){
filtered.CpG <- unique(unlist(data$filtered.CpG.list))
if(is.null(data$minfi.Set$rgSet)){
stop("\n[MeCall]-!!ERROR!! : RGchannel-Set is not founded.")
}

message("\n[MeCall]-[notice] : [SWAN] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) and mSet will be returned.")

rgSet <- data$minfi.Set$rgSet
mSet <- minfi::preprocessRaw(rgSet)
normed <- preprocessSWAN(rgSet,mSet)
normed <- f.out.b.norm(normed,filtered.CpG)

beta.val <- toBeta(normed,offset)
M.val <- toM(normed,offset)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= list(rgSet = NULL, mSet = normed),filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)

return(OUTobj)}

# [ENmix]------------------------------------
ENmix.f <- function(data,offset){
filtered.CpG <- unique(unlist(data$filtered.CpG.list))
if(is.null(data$minfi.Set$rgSet)){
stop("\n[MeCall]-!!ERROR!! : RGchannel-Set is not founded.")
}
rgSet <- data$minfi.Set$rgSet
message("\n[MeCall]-[notice] : [Enmix] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) and mSet will be returned.")
normed <- ENmix::preprocessENmix(rgSet, bgParaEst="oob", dyeCorr="RELIC", QCinfo=NULL, exCpG=filtered.CpG, nCores=4)

normed <- f.out.b.norm(normed,filtered.CpG)

beta.val <- toBeta(normed,offset)
M.val <- toM(normed,offset)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= list(rgSet = NULL, mSet = normed),filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)
return(OUTobj)}

# [dasen]------------------------------------
dasen.f <- function(data,InterpolatedXY,offset){
filtered.CpG <- unique(unlist(data$filtered.CpG.list))
if(is.null(data$minfi.Set$mSet)){
stop("\n[MeCall]-!!ERROR!! : Methyl-Set is not founded.")
}

message("\n[MeCall]-[notice] : [dasen] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) will be returned.")
mSet <- f.out.b.norm(data$minfi.Set$mSet,filtered.CpG)
mn<- mSet@assays@data@listData$Meth
un <- mSet@assays@data@listData$Unmeth
Slide_array <- paste0(mSet@colData@listData$Slide,"_",mSet@colData@listData$array)
names(Slide_array) <- mSet@colData@listData$Sample_Name
colnames(mn) <- Slide_array[colnames(mn)]
colnames(un) <- Slide_array[colnames(un)]
Manifest <- callmanifest(data$TAG)
onetwo <- subset(Manifest, rownames(Manifest) %in% rownames(mn))[rownames(mn),]$Infinium_Design_Type
names(onetwo) <- rownames(mn)
if(InterpolatedXY){
message("\n[MeCall]-[notice] : [interpolatedXY] method is activated. This function utilizes the [adjustedFunnorm] function included in wateRmelon R package. For more details, please refer to the following paper: [Yucheng Wang, Bioinformatics, 2022].")

if(data$TAG == "EPICv2"){CHRs <- subset(Manifest, rownames(Manifest) %in% rownames(mn))[rownames(mn),]$CHR
} else {
CHRs <- subset(Manifest, rownames(Manifest) %in% rownames(mn))[rownames(mn),]$CHR_hg38
}
names(CHRs) <- rownames(mn)
beta.n <- wateRmelon::adjustedDasen(mn = mn, un=un, onetwo=onetwo, chr = CHRs, fudge = offset, ret2=FALSE)
colnames(beta.n) <- names(colnames(beta.n))
} else {
beta.n <- dasen(mns = mn, uns = un, onetwo=onetwo, fudge = offset, ret2=FALSE)
colnames(beta.n) <- names(colnames(beta.n))
}

mset.n <- beta.n
beta.val <- toBeta(normed,offset)
M.val <- toM(normed,offset)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= list(rgSet = NULL, mSet = mset.n),filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)

return(OUTobj)}

# [Noob]------------------------------------
Noob.f <- function(data,offset){
filtered.CpG <- unique(unlist(data$filtered.CpG.list))
if(is.null(data$minfi.Set$rgSet)){
stop("\n[MeCall]-!!ERROR!! : RGchannel-Set is not founded.")
}

message("\n[MeCall]-[notice] : [Noob] method was selected for normalization procedure.")
message("                    Methylation level matrices (Beta and M) and mSet will be returned.")
normed <- preprocessNoob(rgSet=data$minfi.Set$rgSet, offset = offset, dyeCorr = TRUE, verbose = FALSE, dyeMethod="single")

mset.n <- f.out.b.norm(normed,filtered.CpG)
beta.val <- toBeta(mset.n,offset)
M.val <- toM(mset.n,offset)

OUTobj <- list(beta = beta.val, M = M.val, pd = data$pd, TAG = data$TAG, minfi.Set= list(rgSet = NULL, mSet = mset.n),filtered.Sample.list = data$filtered.Sample.list, filtered.CpG.list = data$filtered.CpG.list)

return(OUTobj)}

