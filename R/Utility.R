#---------- Translate p-value -------------
Translate.Pvalue <- function(Nominal.P.list = NULL, cutoff = 0.05,Original.method = c("none","BF", "BH"), return.method = c("none","BF", "BH")){

Methods <- c("none","BF", "BH")

Ori.idx <- which(Methods == Original.method)
re.idx <- which(Methods == return.method)

ALL.Pval <- sort(Nominal.P.list,decreasing=FALSE)
BF.Pval <- p.adjust(ALL.Pval,method = "bonferroni")
BH.Pval <- p.adjust(ALL.Pval,method = "BH")

Pval.Mat <- data.frame(nominal = ALL.Pval, bonferroni = BF.Pval, Benjamini = BH.Pval)

Target <- tail(subset(Pval.Mat, Pval.Mat[,Ori.idx] < cutoff),1)

Return.pval <- Target[,re.idx]

return(Return.pval)}


#---------- Calculate detection P-value  -------------
caldetP <- function(x){
tII.detP <- 1 - pnorm(Intensity[typeII.CpGs,x], mean=rmd[x]+gmd[x], sd=rsd[x]+gsd[x])
tI.R.detP <- 1 - pnorm(Intensity[typeI.Red.CpGs,x], mean=rmd[x]*2, sd=rsd[x]*2)
tI.G.detP <- 1 - pnorm(Intensity[typeI.Grn.CpGs,x], mean=gmd[x]*2, sd=gsd[x]*2)
dd <- c(tII.detP,tI.R.detP,tI.G.detP)
return(dd)}

#---------- Speed up function for finding common probes across IDAT files -------------
SetOp.f <- function (a, b, no.dup.guaranteed = FALSE) {
  if (no.dup.guaranteed) {
    au <- a
    bu <- b
  } else {
    au <- unique(a)
    bu <- unique(b)
  }
  ind <- fastmatch::fmatch(bu, au, nomatch = 0)
  INTERSECT <- au[ind]
  return(INTERSECT)}

#---------- Read RAW IDAT files -------------
All.Quant <- function(x){Q <- illuminaio::readIDAT(x)$Quants[,c('Mean','NBeads')]
return(Q)}

#---------- Functions for readidat-parallel mode -------------
workfn.Q <- function(x){lapply(x,All.Quant)}
workfn.com <- function(i){Reduce("SetOp.f", lapply(i, function(x){rownames(x)}))}
workfn.D <- function(x){sapply(x,caldetP)}

#---------- Remove filtered out probes from mSet -------------

f.out.b.norm <- function(mSet, filtered.CpG){
mSet@assays@data@listData$Meth <- subset(mSet@assays@data@listData$Meth, !rownames(mSet@assays@data@listData$Meth) %in% filtered.CpG)
mSet@assays@data@listData$Unmeth <- subset(mSet@assays@data@listData$Unmeth, !rownames(mSet@assays@data@listData$Unmeth) %in% filtered.CpG)
mSet@NAMES <- rownames(mSet@assays@data@listData$Unmeth)
mSet@elementMetadata@nrows <- nrow(mSet@assays@data@listData$Unmeth)
return(mSet)}

#---------- Calculate Beta-value or M-value from mSet using Meth matrix and Unmeth matrix -------------

toM <- function(mSet,offset){
m.mat <- mSet@assays@data@listData$Meth
u.mat <- mSet@assays@data@listData$Unmeth
M.value <- log2((m.mat + offset) / (u.mat + offset))
return(M.value)
}

toBeta <- function(mSet,offset){
m.mat <- mSet@assays@data@listData$Meth
u.mat <- mSet@assays@data@listData$Unmeth
Beta.value <- (m.mat) / (m.mat + u.mat + offset)
return(Beta.value)
}

#---------- Convert between two types of methylation level -------------

MtoBeta <- function(M){
beta <- 2^M/(2^M + 1)
return(beta)
}

BetatoM <- function(beta){
M <- log2(beta/(1-beta))
return(M)
}

#---------- Load illumina manifest data frame -------------
callmanifest <- function(arraytype){
if(arraytype == "450K"){
#data(IlluminaHumanMethylation450K)
#data(HM450.hg38.manifest)
HM450.hg38.manifest <- HM450.hg38.manifest[IlluminaHumanMethylation450K$Name,]
manifest <- cbind(IlluminaHumanMethylation450K,HM450.hg38.manifest[,1:3])
}else if(arraytype == "EPICv1"){
#data(IlluminaHumanMethylationEPIC.V1b5.hg38)
manifest <- IlluminaHumanMethylationEPIC.V1b5.hg38
}else if(arraytype == "EPICv2"){
#data(IlluminaHumanMethylationEPIC.V2a1.hg38)
manifest <- IlluminaHumanMethylationEPIC.V2a1.hg38
}else{
stop("\n[MeCall]-!!ERROR!! : Unknown array type. Manifest file can not be loaded.")
}
return(manifest)}


callcontrolprobe <- function(arraytype){
#data(Control_probe_chain_file)
colnms <- c("Name","color","control.type")
if(arraytype == "450K"){
subs <- Control_probe_chain_file[which(!is.na(Control_probe_chain_file[,"Methyl450K_Name"])),]
Cp.manifest <- subs[,c(6:8)]
colnames(Cp.manifest) <- colnms
rownames(Cp.manifest) <- subs[,5]

}else if(arraytype == "EPICv1" | arraytype == "EPICv2"){
subs <- Control_probe_chain_file[which(!is.na(Control_probe_chain_file[,"EPIC_Name"])),]
Cp.manifest <- subs[,c(2:4)]
colnames(Cp.manifest) <- colnms
rownames(Cp.manifest) <- subs[,1]

}else{
stop("\n[MeCall]-!!ERROR!! : Unknown array type. Manifest file can not be loaded.")
}
return(Cp.manifest)}



#---------- Prediction of appropriate principal component for estimating batch effect size ------------
# In FindBatches function
pred.n.comp <- function(pca){
cum.var <- (pca$sdev)^2/sum((pca$sdev)^2)
deg <- abs(cum.var[2:length(cum.var)] - cum.var[1:length(cum.var)-1])
ddeg <- deg[1:length(deg)-1] / deg[2:length(deg)]
n <- which.max(ddeg) + 1
return(n)
}

#---------- Taking Case Control setting input from user  -------------
# In CheckMultipleGroup[useful] function
YorN <- function(multiple.group){
if(NA %in% multiple.group){
stop("\n[MeCall]-!!ERROR!! : There is [NA] in your interest variable. Please check your data again.")
}
result <- NA
if(length(multiple.group)==2){
c <- c("y","n")
while (is.na(result)|| !result %in% c){
result <- readline('\n[MeCall]-[question]   Is it correct? (y/n) : ')
}
}else{
while(any(is.na(result)) || !all(result %in% multiple.group)){
result <- readline('\n[MeCall]-[question]  select case group : ')
result <- unlist(strsplit(result,split=','))}
result.case <- result
result <- NA

while(any(is.na(result)) || !all(result %in% multiple.group) || any(result %in% result.case)){
result <- readline('\n[MeCall]-[question]   select control group : ')
result <- unlist(strsplit(result,split=','))}
result.con <- result
result <- list(case = result.case, control = result.con)
}
return(result)}


#----------- Checking current case control setting and modifying setting-------------
# In RemoveBatch function
CheckMultipleGroup <- function(pd, interest){
Groups <- unique(pd[,interest])
if(length(Groups) > 2){message("\n[MeCall]-[notice] : There are multiple phenotypes.")
message("\n[MeCall]-[notice] : **!!CAUTION!!** - Use only comma(,) to separate groups.  <<ex) Case1,Case2,Case3>>")
message("\n[MeCall]-[notice] : **!!CAUTION!!** - Don't use other character to separate groups.")
message("\n[MeCall]-[notice] : Phenotypes in your data -> *[",paste(Groups,collapse=', '),"]*")
CaseCon <- YorN(multiple.group=Groups)
Case <- CaseCon$case
Con <- CaseCon$control
message("\n[MeCall]-[notice] : You select ",paste(Case,collapse=', ')," as Case group.")
message("\n[MeCall]-[notice] : You select ",paste(Con,collapse=', ')," as Control group.")
message("\n[MeCall]-[notice] : Code will automatically update sample group status according to your selection.")
fpd <- pd
rownames(fpd) <- pd$Sample_Name
for(i in 1:length(Case)){
        fpd[which(fpd[,interest]==Case[i]),interest] <- "Case"}
for(i in 1:length(Con)){
        fpd[which(fpd[,interest]==Con[i]),interest] <- "Base"}
Caserow <- rownames(fpd[which(fpd[,interest]=="Case"),])
Conrow <- rownames(fpd[which(fpd[,interest]=="Base"),])
select <- c(Caserow,Conrow)
message("\n[MeCall]-[notice] : Total sample number : ",length(select))
message("\n[MeCall]-[notice] : Sample number of selected case group : ",length(Caserow))
message("\n[MeCall]-[notice] : Sample number of selected control group : ",length(Conrow))
}else if(length(Groups) == 2){
message("\n[MeCall]-[notice] : Detected phenotype : ",paste(Groups,collapse=', '))
temptest <- as.factor(pd[,interest])
fpd <- pd
message("\n[MeCall]-[notice] : Check your Case-Control setting. Control Group should be printed first.")
message("\n[MeCall]-[notice] : Group names : Control / Case -> " ,paste(names(summary(temptest)),collapes=' / '))
message("\n[MeCall]-[notice] : Number of sample for each Group : ",paste(summary(temptest),collapes=' / '))
CaseCon <- YorN(multiple.group = Groups)
if(CaseCon == "n"){message("\n[MeCall]-[notice] : Code will fix your input value.")
Case <- names(summary(temptest))[1]
Con <- names(summary(temptest))[2]
for(i in 1:length(Case)){
        fpd[which(fpd[,interest]==Case[i]),interest] <- "Case"}
for(i in 1:length(Con)){
        fpd[which(fpd[,interest]==Con[i]),interest] <- "Base"}
}else{
Case <- names(summary(temptest))[2]
Con <- names(summary(temptest))[1]
}
} else {stop("\n[MeCall]-!!ERROR!! : There is only one group in interest variable.")}
G <- c(paste(Case,collapse=', '), paste(Con,collapse=', '))
names(G) <- c("Case","Control")
pdobj <- list(modified.pd = fpd, original.pd = pd, groups = G)
return(pdobj)}


#----------- Check methylation level type -------------
# In cellcomp function
ismethyltype <- function(meth){
m <- as.matrix(meth)
u0 <- length(which(m < 0))
u1 <- length(which(m > 1))
if((u0 + u1) ==0){t <- "beta"
}else{t <- "M"}
return(t)
}

#----------- Regression for adjustment using linear regression method-------------
# In RemoveBatch / cellcomp  function
R4A <- function(b, pd, variable, type="beta"){
beta.lm <- apply(b, 1, function(x){
frm <- as.formula(paste('x',paste('~', paste(c(variable), collapse='+'))))
dem <- as.data.frame(pd[colnames(b),])
colnames(dem) <- colnames(pd)
lm(frm,data=dem)})
residuals <- t(sapply(beta.lm,function(x){residuals(summary(x))}))
colnames(residuals)<-colnames(b)
adjusted <- residuals+matrix(apply(b, 1, mean),nrow=nrow(residuals), ncol=ncol(residuals))
if(type == "beta"){
mins <- apply(adjusted,2,function(c){min(c[c>0])})
maxs <- apply(adjusted,2,function(c){max(c[c<1])})
adjusted <- lapply(1:length(mins), function(i){replace(adjusted[,i],which(adjusted[,i] <= 0),mins[i])})
adjusted <- lapply(1:length(maxs), function(i){replace(adjusted[[i]],which(adjusted[[i]] >= 1),maxs[i])})
adjusted <- do.call(cbind,adjusted)
colnames(adjusted) <- colnames(b)
}
return(adjusted)}



#---------- Taking Case Control setting input from user -------------
detectOutlier <- function(x,data,conf){
D <- data[x,]
med <- apply(D[,1:2], 2, median)
shape <- cov(D[,1:2])
dist <- stats::mahalanobis(D[,1:2], med, shape)
cutoff <- qchisq(p=conf, df=2)
outlier <- dist[dist > cutoff]
return(outlier)}

#---------- Taking Case Control setting input from user -------------
EstimateDist <- function(data,outlierinfo,Suggesive,conf){
Outlier <- rep("Inlier",times=nrow(data))
names(Outlier) <- rownames(data)
if(length(outlierinfo) != 0){
temp_D <- subset(data, !rownames(data) %in% outlierinfo)
Outlier[outlierinfo] <- "Outlier"
}else{temp_D <- data}
med <- apply(temp_D[,1:2], 2, median)
shape <- cov(temp_D[,1:2])
dist <- mahalanobis(data[,1:2], med, shape)
cutoff <- qchisq(p=conf, df=2)
ellipse <- car::ellipse(center = med, shape = shape, radius = sqrt(cutoff), segments = 150, draw = FALSE)
ellipse <- as.data.frame(ellipse)
colnames(ellipse) <- colnames(data[,1:2])
if(length(Suggesive) != 0){
Outlier[Suggesive] <- "Suggestive"}
f.data <- as.data.frame(cbind(data, Outlier))
for (i in 1:ncol(data)){f.data[,i] <- as.numeric(f.data[,i])}
MahalaDistList <- list(median=med, shape = shape, dist = dist,ellipse = ellipse, Data = f.data)
return(MahalaDistList)}


#---------- Taking Case Control setting input from user -------------
EstimatePossibleOutlier <- function(data, nrep, fraction, conf){

frac <- round(nrow(data)*fraction,0)

indeces <- lapply(1:nrep, function(x){sort(sample(nrow(data), frac))})

outlier.proportion <- lapply(indeces, function(x){detectOutlier(x,data,conf)})

outlier.proportion <- table(unlist(lapply(outlier.proportion,names)))
outlier.proportion <- outlier.proportion/(nrep/100)
return(outlier.proportion)}


#---------- Taking Case Control setting input from user -------------
# Mahalanobis.Dist
DrawMahalanobisPlot <- function(Mahalanobis.Dist){   
if(is.list(Mahalanobis.Dist)){
L.data <- lapply(1:length(Mahalanobis.Dist), function(x){ Group <- rep(names(Mahalanobis.Dist)[x],times=nrow(Mahalanobis.Dist[[x]]$Data))
Mahalanobis.Dist[[x]]$Data <- cbind(Mahalanobis.Dist[[x]]$Data,Group)})

L.ellipse <- lapply(1:length(Mahalanobis.Dist), function(x){ Group <- rep(names(Mahalanobis.Dist)[x],times=nrow(Mahalanobis.Dist[[x]]$ellipse))
Mahalanobis.Dist[[x]]$ellipse <- cbind(Mahalanobis.Dist[[x]]$ellipse,Group)})

Collapsed.data <- do.call(rbind,L.data)
Collapsed.ellipse <- do.call(rbind,L.ellipse)
}else{Group <- rep("1",times=nrow(Mahalanobis.DistData))
Collapsed.data <-cbind(Mahalanobis.Dist$Data, Group)
Collapsed.ellipse<-cbind(Mahalanobis.Dist$ellipse, Group)}
O.O <- subset(Collapsed.data, Outlier == "Outlier")
S.O <- subset(Collapsed.data, Outlier == "Suggestive")
X <- colnames(Collapsed.data)[1]
Y <- colnames(Collapsed.data)[2]
figure <- ggplot(data = Collapsed.data, aes(x= .data[[X]], y= .data[[Y]],color=Group ,Group = Group)) + geom_point(size=1.5) + geom_polygon(data = Collapsed.ellipse,mapping = aes(x = PC1, y = PC2,group = Group, colour = Group),fill = "gray80", alpha = 0.3)+ labs(title = "Outlier detection", subtitle = "Ellipses repersents boundary of outlier", x= colnames(Collapsed.data)[1],y= colnames(Collapsed.data)[2]) + geom_point(data = O.O, mapping = aes(x=.data[[X]], y=.data[[Y]]),colour="red",size = 1.5) + geom_point(data = S.O, mapping = aes(x=.data[[X]], y=.data[[Y]]),colour="blue",size = 1.5) + ggrepel::geom_label_repel(data = O.O, mapping = aes(x= .data[[X]], y= .data[[Y]],label=rownames(O.O)),colour = "red",show.legend = FALSE) + ggrepel::geom_label_repel(data = S.O, mapping = aes(x= .data[[X]], y= .data[[Y]],label=rownames(S.O)),colour = "blue",show.legend = FALSE)
return(figure)}

#---------- Remove failed sample from original minfi set -------------
rmFailedSample <- function(minfiSet,Remained){

minfiSet@colData <- subset(minfiSet@colData,rownames(minfiSet@colData) %in% Remained)

if(sum(methods::is(minfiSet) == "RGChannelSet") == 1){
minfiSet@assays@data@listData$Green <- minfiSet@assays@data@listData$Green[,Remained]

minfiSet@assays@data@listData$Red <- minfiSet@assays@data@listData$Red[,Remained]
}

if(sum(methods::is(minfiSet) == "MethylSet") == 1){
minfiSet@assays@data@listData$Meth <- minfiSet@assays@data@listData$Meth[,Remained]

minfiSet@assays@data@listData$Unmeth <- minfiSet@assays@data@listData$Unmeth[,Remained]
}

return(minfiSet)}


# ---------Returns a reconstructed RGChannelSet without the information from the input data.-------------------------
# The input foramt of Import.obj parameter is a return object of "MeCall.ReadIdat()".
# The input foramt of preselected.obj parameter is a subset of the Illumina manifest.
ReShapeImportObj <- function(Import.obj, preselected.obj){

Ori.rgSet <- Import.obj$minfi.Set$rgSet
ReShape.rgSet <- ReShapergSetObj(Ori.rgSet, preselected.obj)

ReShape.mSet <- minfi::preprocessRaw(ReShape.rgSet)

RemainRows <- ReShape.mSet@NAMES

Re.beta <- Import.obj$beta[RemainRows,]
Re.M <- Import.obj$M[RemainRows,]
Re.intensity <- Import.obj$intensity[RemainRows,]
Re.Meth <- Import.obj$Meth[RemainRows,]
Re.Unmeth <- Import.obj$Unmeth[RemainRows,]
Re.detP <- Import.obj$detP[RemainRows,]
Re.B.count <- Import.obj$B.count[RemainRows,]

Re.minfi.Set <- list(rgSet = ReShape.rgSet, mSet = ReShape.mSet)

ReducedObj <- list(beta = Re.beta, M = Re.M, intensity = Re.intensity, Meth = Re.Meth, Unmeth = Re.Unmeth, detP = Re.detP, B.count = Re.B.count, pd = Import.obj$pd, TAG = Import.obj$TAG, minfi.Set = Re.minfi.Set)
return(ReducedObj)}


ReShapergSetObj <- function(rgSet, preselected.obj){

ReShape.rgSet <- rgSet

Green.ADR <- rownames(rgSet@assays@data@listData$Green)
Red.ADR <- rownames(rgSet@assays@data@listData$Red)

RmADRs <- unique(c(as.character(preselected.obj$AddressA_ID),as.character(preselected.obj$AddressB_ID)))

Reduced.Green.ADR <- setdiff(Green.ADR,RmADRs)
Reduced.Red.ADR <- setdiff(Red.ADR,RmADRs)


ReShape.rgSet.Green <- rgSet@assays@data@listData$Green[Reduced.Green.ADR,]
ReShape.rgSet.Red <- rgSet@assays@data@listData$Red[Reduced.Red.ADR,]

ReShape.rgSet@assays@data@listData$Green <- ReShape.rgSet.Green
ReShape.rgSet@assays@data@listData$Red <- ReShape.rgSet.Red
ReShape.rgSet@NAMES <- Reduced.Green.ADR
ReShape.rgSet@elementMetadata@nrows <- length(Reduced.Green.ADR)

return(ReShape.rgSet)}

