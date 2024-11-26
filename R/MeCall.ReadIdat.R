#' Read the idat file included in the sample data file
#'
#' @import stringr
#' @importFrom xfun file_ext
#' @importFrom fastmatch fmatch
#' @importFrom illuminaio readIDAT
#'
#' @description
#' `MeCall.ReadIdat()` will attempt to find IDAT files located 
#' in the directory specified by `idat.dir` parameter by 
#' combining the Sentrix ID and Sentrix position information 
#' in the sample metadata (`pd.file` parameter).
#'
#' If you want to load only some of the IDAT files present in the directory,
#' delete the information of unnecessary samples from the sample metadata 
#' before running the function.
#'
#' @param idat.dir The directory path where the idat file exists.
#' @param pd.file A file containing Sentrix_ID and Sentrix_Position (csv, txt, tsv).
#' @param offset Integer. A fixed value will be added when generating methylation levels.
#' @param platform The version of the microarray to import (450K, EPICv1, EPICv2)
#'
#' @return A list object containing basic informations to process methylation microarray : 
#' \itemize{beta} {A beta-value matrix defined with columns as samples and rows as CpG sites}
#' \itemize{M} {A M-value matrix structured similarly to a beta-value matrix}
#' \itemize{intensity} {A total intensity matrix structured similarly to a beta-value matrix}
#' \itemize{meth} {A methylated intensity matrix structured similarly to a beta-value matrix}
#' \itemize{Unmeth} {An unmethylated intensity matrix structured similarly to a beta-value matrix}
#' \itemize{detP} {A detection P-value matrix structured similarly to a beta-value matrix}
#' \itemize{B.count} {A bead count matrix structured similarly to a beta-value matrix}
#' \itemize{pd} {A data frame composed of sample metadata, which includes information about each sample}
#' \itemize{TAG} {A vector consisting of strings that indicate the version of the microarray}
#' \itemize{minfi.Set} {A list object consisting of two minfi objects (RGChannelSet Set and MethylSet)}
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [minfi::MethylSet()], [minfi::RGChannelSet()], [ChAMP::champ.import()]
#'
#' @references 
#' Aryee, M.J. et al. (2014). Minfi: a flexible and comprehensive Bioconductor package for the 
#' analysis of Infinium DNA methylation microarrays. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/btu049}
#'
#' Tian, Y. et al. (2017). ChAMP: updated methylation analysis pipeline for Illumina BeadChips. 
#' Bioinformatics. /url{https://doi.org/10.1093/bioinformatics/btx513}
#'
#'
#' @examples
#' \dontrun{
#' data.Import <- MeCall.ReadIdat(idat.dir = "~/IdatDir", pd.file = "~/IdatDir/Sample_Sheet.csv", 
#' offset=100, platform = c("EPICv1"))
#' }
#'
#' @export
MeCall.ReadIdat <- function(idat.dir = NULL, pd.file=NULL,offset=100,platform=c("450K","EPICv1","EPICv2")){

taps <- c('txt','tsv')
message("\n[MeCall]-[notice] : Read sample file and IDAT file directory.")
if (file.exists(pd.file)){
ext <- xfun::file_ext(pd.file)
} else {
stop("\n[MeCall]-!!ERROR!! : The file does not exist. Please check your file again.")
}

if("csv" %in% ext){
pd <-read.csv(pd.file, header=T)
rownames(pd) <- pd$Sample_Name
}else if(ext %in% taps){
pd <-read.table(pd.file,header=T,sep='\t')
rownames(pd) <- pd$Sample_Name
} else {message("\n[MeCall]-!!ERROR!! : Not supported file format : ",ext)
stop("\n[MeCall]-!!ERROR!! : This package supports 3 kind of format : .csv / .txt / .tsv")}

message("\n[MeCall]-[notice] : Comparing IDAT file list from sample data file and IDAT file directory.")
if(!(any(colnames(pd) %in% c("Slide","Sentrix_ID")) & any(colnames(pd) %in% c("array","Sentrix_Position")))){
stop("\n[MeCall]-!!ERROR!! : Your sample file does not include slide or array information.")}
if(any(colnames(pd) %in% "Sentrix_ID")){
colnames(pd)[which(colnames(pd)=="Sentrix_ID")] <- "Slide"}
if(any(colnames(pd) %in% "Sentrix_Position")){
colnames(pd)[which(colnames(pd)=="Sentrix_Position")] <- "array"}


if(any(is.na(pd$Slide))){
message("\n[MeCall]-[WARNING] : NA is detected in Slide column. Please ensure that sufficiently specific information is provided to load the correct IDAT file.")
Mod_slide <- stringr::str_replace_na(pd$Slide, replacement="")
} else {
Mod_slide <- pd$Slide
}

if(any(is.na(pd$array))){
message("\n[MeCall]-[WARNING] : NA is detected in array column. Please ensure that sufficiently specific information is provided to load the correct IDAT file.")
Mod_array <- stringr::str_replace_na(pd$array, replacement="")
} else {
Mod_array <- pd$array
}


IDATnm <- paste(Mod_slide,Mod_array,sep='_')
Grnlist <- list.files(idat.dir,pattern='Grn')
Redlist <- list.files(idat.dir,pattern='Red')

for(i in 1:length(IDATnm)){
if(table(str_detect(Grnlist, pattern = IDATnm[i]))[2]==0){
stop("\n[MeCall]-!!ERROR!! : There is no related Grn_IDAT file in IDAT directory. Check this sample which is contained sample table file : ",IDATnm[i])}
if(table(str_detect(Redlist, pattern = IDATnm[i]))[2]==0){
stop("\n[MeCall]-!!ERROR!! : There is no related Red_IDAT file in IDAT directory. Check this sample which is contained sample table file : ",IDATnm[i])}}
message("\n[MeCall]-[notice] : Check matching Done. MethylCallR found ",length(IDATnm)," IDAT files.")

message("\n[MeCall]-[notice] : Loading IDAT files.")
if(str_sub(idat.dir,start=(str_length(idat.dir)),end=(str_length(idat.dir)))!= "/"){
        idat.dir <- paste0(idat.dir,'/')}
act.Grn.idx <- unlist(lapply(IDATnm,function(x) {grep(x,Grnlist)}))
act.Red.idx <- unlist(lapply(IDATnm,function(x) {grep(x,Redlist)}))
Gdat <- paste0(idat.dir,Grnlist[act.Grn.idx])
Rdat <- paste0(idat.dir,Redlist[act.Red.idx])
Adat <- c(Gdat, Rdat)

sample.pair <- data.frame(sample_Name = pd$Sample_Name,IDATadress = IDATnm,row.names=2)
sample.num <- nrow(sample.pair)

#if(PA){
#message("\n[MeCall]-[notice] : Parallel Mode On. ",length(Cluster)," cores will be used to read IDAT file.")
#idx <- split(Adat, sort(rep_len(1:length(Cluster), length(Adat))))
#rParts <- clusterApply(Cluster, idx, workfn.Q)
#rcoms <- clusterApply(Cluster, rParts, workfn.com)
#All <- unlist(rParts,recursive = FALSE)
#com <- Reduce("SetOp.f",rcoms)
#} else{
message("\n[MeCall]-[notice] : Single Core Mode On. 1 core will be used to read IDAT file.")
All <- lapply(Adat,All.Quant)
com <- Reduce("SetOp.f", lapply(All, function(x){rownames(x)}))
#}
message("\n[MeCall]-[notice] : ",length(com), " common probes were found.")

All.idat.mean <- lapply(All, function(x){X <- x[com,"Mean"]
return(X)})

All.idat.NBeads <- lapply(All, function(x){X <- x[com,"NBeads"]
return(X)})

All.idat.mean <- do.call(cbind,All.idat.mean)
All.idat.NBeads <- do.call(cbind,All.idat.NBeads)

Gmean <-All.idat.mean[,1:sample.num]
Rmean <-All.idat.mean[,(sample.num+1):ncol(All.idat.mean)]

GBead <-All.idat.NBeads[,1:sample.num]
RBead <-All.idat.NBeads[,(sample.num+1):ncol(All.idat.NBeads)]

#rownames(Gmean) <- com
#rownames(Rmean) <- com
colnames(Gmean) <- sample.pair[IDATnm,"sample_Name"]
colnames(Rmean) <- sample.pair[IDATnm,"sample_Name"]
#rownames(GBead) <- com
#rownames(RBead) <- com
colnames(GBead) <- sample.pair[IDATnm,"sample_Name"]
colnames(RBead) <- sample.pair[IDATnm,"sample_Name"]




if(platform == "EPICv1"){
message("\n[MeCall]-[notice] : Extract 'Negative control probe' information from EPICv1 B5 Manifest.")
rg.annotation <- c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
}else if(platform == "EPICv2"){
message("\n[MeCall]-[notice] : Extract 'Negative control probe' information from EPICv2 A1 Manifest.")
rg.annotation <- c(array="IlluminaHumanMethylationEPICv2", annotation="20a1.hg38")
}else if(platform == "450K"){
message("\n[MeCall]-[notice] : Extract 'Negative control probe' information from 450k Manifest.")
rg.annotation <- c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
}else{
stop("\n[MeCall]-!!ERROR!! : Array platform name Unknown.")
}

Control.mani <- callcontrolprobe(platform)
NCprobes <- rownames(subset(Control.mani,Control.mani$Name == "NEGATIVE"))
NCprobes <- intersect(NCprobes,com)

RGset <- RGChannelSet(Red = Rmean, Green = Gmean,annotation=rg.annotation)
RGset@colData@listData <- as.list(pd)


Gcon <- Gmean[NCprobes,]
Rcon <- Rmean[NCprobes,]
gmd <- apply(Gcon,2,median)
gsd <- apply(Gcon,2,mad)
rmd <- apply(Rcon,2,median)
rsd <- apply(Rcon,2,mad)

if(!platform %in% c("450K","EPICv1","EPICv2")){
stop("\n[MeCall]-!!ERROR!! : ",platform," : Unknown Array platform.")
}else{mani <- callmanifest(platform)}

if(platform == "EPICv1"){
typeII <- intersect(subset(mani,Color_Channel=="G+R")$AddressA_ID,com)
typeI.R.A <- intersect(subset(mani,Color_Channel=="Red")$AddressA_ID,com)
typeI.R.B <- intersect(subset(mani,Color_Channel=="Red")$AddressB_ID,com)
typeI.G.A <- intersect(subset(mani,Color_Channel=="Grn")$AddressA_ID,com)
typeI.G.B <- intersect(subset(mani,Color_Channel=="Grn")$AddressB_ID,com)
}else if(platform == "EPICv2"){

typeII <- intersect(subset(mani,Color_Channel=="")$AddressA_ID,com)
typeI.R.A <- intersect(subset(mani,Color_Channel=="Red")$AddressA_ID,com)
typeI.R.B <- intersect(subset(mani,Color_Channel=="Red")$AddressB_ID,com)
typeI.G.A <- intersect(subset(mani,Color_Channel=="Grn")$AddressA_ID,com)
typeI.G.B <- intersect(subset(mani,Color_Channel=="Grn")$AddressB_ID,com)

}else{
typeII <- intersect(subset(mani,Color_Channel=="G+R")$AddressA_ID,com)
typeI.R.A <- intersect(subset(mani,Color_Channel=="Red")$AddressA_ID,com)
typeI.R.B <- intersect(subset(mani,Color_Channel=="Red")$AddressB_ID,com)
typeI.G.A <- intersect(subset(mani,Color_Channel=="Grn")$AddressA_ID,com)
typeI.G.B <- intersect(subset(mani,Color_Channel=="Grn")$AddressB_ID,com)
}


CpG.addressA.vec <- as.character(mani$AddressA_ID)
CpG.addressB.vec <- as.character(mani$AddressB_ID)
names(CpG.addressA.vec) <- rownames(mani)
names(CpG.addressB.vec) <- rownames(mani)

M.address <- c(typeII, typeI.R.B, typeI.G.B)
M <- rbind(Gmean[typeII,],Rmean[typeI.R.B,],Gmean[typeI.G.B,])
M.CpG <- c(names(CpG.addressA.vec[(CpG.addressA.vec %in% typeII)]),names(CpG.addressB.vec[(CpG.addressB.vec %in% typeI.R.B)]),names(CpG.addressB.vec[(CpG.addressB.vec %in% typeI.G.B)]))


U.address <- c(typeII, typeI.R.A, typeI.G.A)
U <- rbind(Rmean[typeII,],Rmean[typeI.R.A,],Gmean[typeI.G.A,])
U.CpG <- c(names(CpG.addressA.vec[(CpG.addressA.vec %in% typeII)]),names(CpG.addressA.vec[(CpG.addressA.vec %in% typeI.R.A)]),names(CpG.addressA.vec[(CpG.addressA.vec %in% typeI.G.A)]))



rownames(M) <- M.CpG
rownames(U) <- U.CpG
Mset<- MethylSet(Meth=M,Unmeth=U,annotation=rg.annotation)
Mset@colData@listData <- as.list(pd)

message("\n[MeCall]-[notice] : Generate main matrix.")
# M matrix
Meth = M
# U matrix
Unmeth = U
# beta value
Bval = M/(U+M+offset)
# M value
Mval = log2((M+offset)/(U+offset))
Intensity = U + M

# detection P matrix
typeII.CpGs <- names(CpG.addressA.vec[(CpG.addressA.vec %in% typeII)])
typeI.Red.CpGs <- names(CpG.addressA.vec[(CpG.addressA.vec %in% typeI.R.A)])
typeI.Grn.CpGs <- names(CpG.addressA.vec[(CpG.addressA.vec %in% typeI.G.A)])

message("\n[MeCall]-[notice] : Calculate detection P-value.")

detP_list <- lapply(colnames(Intensity), function(x){
TypeII_detP <- pnorm(Intensity[typeII.CpGs,x], mean = rmd[x] + gmd[x], sd=rsd[x] + gsd[x])
TypeIR_detP <- pnorm(Intensity[typeI.Red.CpGs,x], mean = rmd[x]*2, sd=rsd[x]*2)
TypeIG_detP <- pnorm(Intensity[typeI.Grn.CpGs,x], mean = gmd[x]*2, sd=gsd[x]*2)
A_detP <- c(TypeII_detP, TypeIR_detP, TypeIG_detP)
return(A_detP)})

detectP <- 1 - as.matrix(do.call(cbind, detP_list))
detectP <- detectP[rownames(Intensity),]
colnames(detectP) <- colnames(Intensity)


if(sum(is.na(detectP))) {message("\n[MeCall]-[notice] : CAUSION! -- There are NAs in detection P-value matrix.")}


message("\n[MeCall]-[notice] : Generate Bead count matrix.")
Mbead <- RBead[M.address,]
Ubead <- RBead[U.address,]
Mbead[Mbead < 3 | Ubead < 3] <- NA
rownames(Mbead) <- M.CpG
colnames(Mbead) <- colnames(Intensity)

message("\n[MeCall]-[notice] : Reading Idat file Done.")
out <- list(beta = Bval, M=Mval, intensity = Intensity, Meth = Meth, Unmeth = Unmeth, detP = detectP, B.count = Mbead, pd = pd, TAG = platform, minfi.Set=list(rgSet = RGset, mSet = Mset))
return(out)}
