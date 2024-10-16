#' Perform differential methylation analysis to identify differentially methylated regions (DMRs).
#'
#' @import bumphunter
#' @import GenomicRanges
#' @import limma
#' @import DMRcate
#'
#' @description
#' This function provides multiple useful algorithms for identifying 
#' differentially methylated regions.
#'
#' `MeCall.DMR()` function supports 5 algorithms :
#' * Bumphunter (`method = "bumphunter"`)
#' * seqlm (`method = "seqlm"`)
#' * dmrff (`method = "dmrff"`)
#' * DMRcate (`method = "DMRcate"`)
#' * combined-pvalues (`method = "combP"`)
#'
#' Unlike other algorithms, `combined-pvalues` is a python-based algorithm.
#' Instead of returning DMR results, `MeCall.DMR()` will return the basic 
#' file format (BED format) and a guide message required for using the 
#' algorithm. For detailed instructions, please refer to GitHub : 
#' https://github.com/brentp/combined-pvalues
#'
#' @param meth A matrix of methylation level.
#' @param pd A data frame of modified sample data file from `MeCall.RemoveBatch()`.
#' @param interest The column name in the sample data file containing the phenotype of interest.
#' @param covariate Column names of covariates to be considered when identifying differentially methylated regions.
#' @param arraytype Version of the methylation microarray.
#' @param method The name of the algorithm to be used for identifying differentially methylated regions.
#' @param adjPvalDmr The p-value threshold for identifying significant differentially methylated regions.
#' @param nProbes The minimum number of probes included within a significant differentially methylated region.
#' @param probe_gap The minimum distance between probes to be considered as different regions. If set to NULL, the recommended or default values for each algorithm will be automatically applied. (bumphunter : maxGap; seqlm : max_dist; dmrff : maxgap; DMRcate : lambda) 
#' @param cutoff (Bumphunter) Threshold values to identify candidate regions. Upper and lower threshold values can be provided individually.
#' @param pickCutoff (Bumphunter) logical. If set to `TRUE`, an appropriate cutoff is selected using the permutation distribution.
#' @param pickCutoffQ (Bumphunter) The quantile value for the cutoff to be selected from the permutation distribution.
#' @param smooth (Bumphunter) logical. If set to `TRUE`, Smoothing will be performed to the estimated profile using the smoothFunction.
#' @param smoothFunction (Bumphunter) The name of the function to be applied for smoothing.
#' @param useWeights (Bumphunter) logical. If set to `TRUE` and smoothing function is set to `loessByCluster`, the standard errors of the point-wise estimates of the profile function will be used as weights.
#' @param B (Bumphunter) The number of resamples for generating the null distribution.
#' @param nullMethod (Bumphunter) The method for generating the null distribution.
#' @param max_block_length (seqlm) The maximum number of CpGs that can be included in one region.
#' @param dmp.fdr (DMRcate) The FDR threshold (Benjamini-Hochberg method) for identifying significant differential methylation probes.
#' @param C (DMRcate) Scaling factor for bandwidth. Lambda/C is defined as sigma, and a Gaussian kernel is computed accordingly.
#'
#' @return A data frame containing the genomic locations, associated genes, 
#' and major statistical information of the identified DMR.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [bumphunter::bumphunter()], [DMRcate::dmrcate()]
#'
#' @references
#' Default setting of algorithms based on
#' Mallik, S. et al. (2019). An evaluation of supervised methods for identifying differentially 
#' methylated regions in Illumina methylation arrays. Briefings in Bioinformatics. 
#' /url{https://doi.org/10.1093/bib/bby085}
#' 
#' For bumphunter
#' Jaffe, A.E. et al. (2012). Bump hunting to identify differentially methylated regions in 
#' epigenetic epidemiology studies. International Journal of Epidemiology. 
#' /url{https://doi.org/10.1093/ije/dyr238}
#'
#' For dmrff (BioRxiv)
#' Suderman, M. et al. (2018). dmrff: identifying differentially methylated regions efficiently 
#' with power and control. /url{https://doi.org/10.1101/508556}
#'
#' For DMRcate
#' Peters, T.J. et al. (2015). De novo identification of differentially methylated regions in the 
#' human genome. Epigenetics & Chromatin. /url{https://doi.org/10.1186/1756-8935-8-6}
#'
#' For seqlm
#' Kolde, R. et al. (2016). seqlm: an MDL based method for identifying differentially methylated 
#' regions in high density methylation array data. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/btw304}
#'
#' For comb-P
#' Pedersen, B.S. et al. (2012). Comb-p: software for combining, analyzing, grouping and correcting
#' spatially correlated P-values. Bioinformatics. https://doi.org/10.1093/bioinformatics/bts545
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.DMR] - Perform Bumphunter with default setting.
#' data.DMRs.Bump <- MeCall.DMR(meth = data.rmBatch$methylation.level, pd = data.rmBatch$modified.pd
#' , interest="Sample_Group", covariate = NULL, arraytype = c("EPICv2"), method = c("bumphunter"), 
#' adjPvalDmr=1, nProbes=2, probe_gap = 250, cutoff = NULL, pickCutoff = TRUE, pickCutoffQ = 0.95, 
#' smooth=TRUE, smoothFunction=loessByCluster, useWeights=FALSE, B=100, nullMethod="bootstrap")
#'
#' # Run [MeCall.DMR] - Perform seqlm with default setting.
#' data.DMRs.seq <- MeCall.DMR(meth = data.rmBatch$methylation.level, pd = data.rmBatch$modified.pd,
#' interest="Sample_Group", covariate = NULL, arraytype = c("EPICv2"), method = c("seqlm"), 
#' adjPvalDmr=1, nProbes=2, probe_gap = 1000, max_block_length=50)
#'
#' # Run [MeCall.DMR] - Perform dmrff with default setting.
#' data.DMRs.dmrff <- MeCall.DMR(meth = data.rmBatch$methylation.level, 
#' pd = data.rmBatch$modified.pd, interest="Sample_Group", covariate = NULL, 
#' arraytype = c("EPICv2"), method = c("dmrff"), adjPvalDmr=1, nProbes=2, probe_gap = 500)
#'
#' # Run [MeCall.DMR] - Perform DMRcate with default setting.
#' data.DMRs.DMRcate <- MeCall.DMR(meth = data.rmBatch$methylation.level, 
#' pd = data.rmBatch$modified.pd, interest="Sample_Group", covariate = NULL, 
#' arraytype = c("EPICv2"), method = c("DMRcate"), adjPvalDmr=1, nProbes=2, probe_gap = 500, 
#' dmp.fdr=0.05, C=5)
#'
#' # Run [MeCall.DMR] - Perform combP with default setting.
#' data.DMRs.combP.BEDfile <- MeCall.DMR(meth = data.rmBatch$methylation.level, 
#' pd = data.rmBatch$modified.pd, interest="Sample_Group", covariate = NULL, 
#' arraytype = c("EPICv2"), method = c("combP"))
#'}
#'
#' @export
MeCall.DMR <- function(meth = data.rmBatch$methylation.level, pd = data.rmBatch$modified.pd, interest="Sample_Group", covariate = NULL, arraytype = c("EPICv2","EPICv1","450K"), method = c("bumphunter", "seqlm", "dmrff", "DMRcate", "combP"), adjPvalDmr=1, nProbes=2, probe_gap = NULL, cutoff = NULL, pickCutoff = TRUE, pickCutoffQ = 0.95, smooth=TRUE, smoothFunction=loessByCluster, useWeights=FALSE, B=100, nullMethod="bootstrap", max_block_length=50, dmp.fdr=0.05, C=5){

mani <- callmanifest(arraytype)

mtype <- ismethyltype(meth)
if(mtype == "beta"){
meth <- replace(meth,which(meth <= 0),0.0001)
meth <- replace(meth,which(meth >= 1),0.9999)
}

method <- match.arg(method, c("bumphunter", "seqlm", "dmrff", "DMRcate", "combP"), several.ok = FALSE)


if(method == "bumphunter"){ 
message("------------------------------------------------------------------")
message("\n[MeCall]-[NOTICE] : Bumphunter algorithm Start.")
message("\n[MeCall]-[NOTICE] : This code based on DMR function in ChAMP.")
message("\n[MeCall]-[NOTICE] : If you want to know more detail information about code, you should visit ChAMP documentation.")
#message("\n[MeCall]-[NOTICE] : Bumphunter use ",cores," cores to calculate DMR.")

if(is.null(probe_gap)){
probe_gap = 250
}

CHR.levels <- paste0("chr",c(1:22,"X","Y","M"))
Anno.idx <- intersect(rownames(meth),rownames(mani))
mani <- mani[Anno.idx,]


if(arraytype =="EPICv2"){
Anno.Info <- mani[!is.na(mani$MAPINFO),]
Anno.Info <- Anno.Info[!is.na(Anno.Info$CHR),]
Anno.Info <- Anno.Info[!grepl("chr0",Anno.Info$CHR),]
Anno.Info$CHR <- factor(Anno.Info$CHR, levels = CHR.levels)
Anno.Info <- Anno.Info[order(Anno.Info$CHR,Anno.Info$MAPINFO),]
Anno.CHR <- as.character(Anno.Info$CHR)
Anno.Pos <- Anno.Info$MAPINFO
Anno.Info$Name <- rownames(Anno.Info)

}

if(arraytype == "EPICv1" | arraytype == "450K"){
Anno.Info <- mani[!is.na(mani$Start_hg38),]
Anno.Info <- Anno.Info[!is.na(Anno.Info$CHR_hg38),]
Anno.Info$CHR_hg38 <- factor(Anno.Info$CHR_hg38, levels = CHR.levels)
Anno.Info <- Anno.Info[order(Anno.Info$CHR_hg38,Anno.Info$Start_hg38),]
Anno.CHR <- as.character(Anno.Info$CHR_hg38)
Anno.Pos <- (Anno.Info$Start_hg38) + 1
}

names(Anno.CHR) <- Anno.Info$Name
names(Anno.Pos) <- Anno.Info$Name
#CpG.idx <- intersect(rownames(meth),rownames(Anno.Info))

probe.cluster <- bumphunter::clusterMaker(Anno.CHR, Anno.Pos, maxGap=probe_gap)
names(probe.cluster) <- Anno.Info$Name
bumphunter.idx <- Anno.idx[which(probe.cluster %in% names(which(table(probe.cluster)>=nProbes)))]

Bump.Anno.CHR <- Anno.CHR[bumphunter.idx]
Bump.Anno.Pos <- Anno.Pos[bumphunter.idx]
Bump.probe.cluster <- probe.cluster[bumphunter.idx]

Y <- meth[bumphunter.idx,]
X <- MeCall.MakeModelMatrix(pd = pd, covariate = covariate, interest = interest, PCA.result = NULL, Cell.comp = NULL, adjustment=FALSE)

if(mtype == "beta"){
Y <- BetatoM(Y)
}

Bumps <- bumphunter::bumphunter(Y, design=X, chr=Bump.Anno.CHR, pos=Bump.Anno.Pos, cluster=Bump.probe.cluster,coef = ncol(X) ,cutoff=cutoff, pickCutoff=pickCutoff, pickCutoffQ=pickCutoffQ, smooth=smooth, smoothFunction=smoothFunction, useWeights=useWeights, permutations=NULL, verbose=TRUE, B=B, nullMethod=nullMethod)

DMRs <- subset(Bumps$table, p.valueArea <= adjPvalDmr)
DMRs <- DMRs[order(DMRs$p.valueArea),]
rownames(DMRs) <- paste("DMR",1:nrow(DMRs),sep="_")
DMRobjs <- data.frame(DMRs[,1:3],width=DMRs[,3]-DMRs[,2],DMRs[,4:14])
colnames(DMRobjs)[1:3] <- c("seqnames","start","end")
message("\n[MeCall]-[NOTICE] : Bumphunter algorithm finished.")
}



if(method == "seqlm"){
if(!requireNamespace("seqlm", quietly = TRUE)){
stop("\n[MeCall]-!!ERROR!! : [seqlm] R package is not founded in your R environment. Please load [seqlm] R package to run this function using seqlm method {library(seqlm)}. You can download [seqlm] R package at following github page : https://github.com/raivokolde/seqlm")
}
message("\n[mkk]-[NOTICE] : Start seqlm algorithm.")

if(!is.null(covariate)){
message("\n[MeCall]-[NOTICE] : Methylation level will be adjusted using linear regression model before applying seqlm algorithm.")
meth <- MeCall.MakeModelMatrix(pd = pd, covariate = covariate, interest = interest, PCA.result = NULL, Cell.comp = NULL, adjustment=TRUE, meth = meth)$adjusted.meth
}

if(is.null(probe_gap)){
probe_gap = 1000
}

if(arraytype =="EPICv2"){
mani <- mani[!is.na(mani$MAPINFO),]
mani[which(mani$Strand_FR == "F"),]$Strand_FR <- "+"
mani[which(mani$Strand_FR == "R"),]$Strand_FR <- "-"
mani[which(mani$Strand_FR == "0"),]$Strand_FR <- "*"

GRange.colnm.idx1 <- c(5,7,8) 
GRange.colnm.idx2 <- c(9,12,10,1) 
GRange.df <- mani[,GRange.colnm.idx1]
GRange.sup.info <- mani[,GRange.colnm.idx2]
colnames(GRange.sup.info) <- c("Genomic.featrue","CGI.feature","gene","CpGid")

seqlm_anno <- GenomicRanges::makeGRangesFromDataFrame(GRange.df, seqnames.field=c("CHR"), start.field = c("MAPINFO"), end.field = c("MAPINFO"), strand.field = c("Strand_FR"))
seqlm_anno@elementMetadata@listData <- as.list(GRange.sup.info)
}

if(arraytype == "EPICv1" | arraytype == "450K"){
mani <- mani[!is.na(mani$CHR_hg38),]
mani[which(mani$Strand == "F"),]$Strand <- "+"
mani[which(mani$Strand == "R"),]$Strand <- "-"
mani[which(mani$Strand == ""),]$Strand <- "*"
Pos <- mani[,'Start_hg38'] + 1
seqlm_anno <- GenomicRanges::GRanges(seqnames=Rle(mani[,'CHR_hg38']),ranges = IRanges(Pos),strand=Rle(mani[,'Strand']),Genomic.featrue=mani[,'UCSC_RefGene_Group'],CGI.feature=mani[,'Relation_to_UCSC_CpG_Island'],gene=mani[,'UCSC_RefGene_Name'],CpGid=mani[,'Name'])
names(seqlm_anno) <- seqlm_anno$CpGid
}

Group.Var <- pd[colnames(meth),interest]
names(Group.Var) <- colnames(meth)

message("\n[MeCall]-[NOTICE] : Calculate seqlm segmentation ")
segments = seqlm::seqlm(values = meth, genome_information = seqlm_anno, annotation = Group.Var, max_block_length=max_block_length, max_dist=probe_gap)
seqlmDMR <- as.data.frame(segments@elementMetadata@listData)
DMRobjs <- seqlmDMR[which(seqlmDMR$fdr <= adjPvalDmr & seqlmDMR$length >= nProbes),]
rownames(DMRobjs) <- paste("DMR",1:nrow(DMRobjs),sep="_")
message("\n[MeCall]-[NOTICE] : Seqlm algorithm Done.")

}

if(method == "dmrff"){
if(!requireNamespace("dmrff", quietly = TRUE)){
stop("\n[MeCall]-!!ERROR!! : [dmrff] R package is not founded in your R environment. Please load [dmrff] R package to run this function using dmrff method {library(dmrff)}. You can download [dmrff] R package at following github page : https://github.com/perishky/dmrff")
}
message("\n[MeCall]-[NOTICE] : dmrff algorithm Start.")


if(is.null(probe_gap)){
probe_gap = 500
}

message("\n[MeCall]-[NOTICE] : dmrff algorithm will be calculate with following parameters.")
message("\n[MeCall]-[NOTICE] : Maxgap = ",probe_gap," / adjusted P-value cutoff = ",adjPvalDmr," / number of probes in region > ",nProbes)


common.cpg <- intersect(rownames(mani),rownames(meth))
mani <- mani[common.cpg,]
design.mat <- MeCall.MakeModelMatrix(pd = pd, covariate = covariate, interest = interest, PCA.result = NULL, Cell.comp = NULL, adjustment=FALSE)

fit1 <- limma::lmFit(meth,design.mat)
fit1 <- limma::eBayes(fit1)
stats <- data.frame(estimate=fit1$coefficients[,ncol(design.mat)], se=sqrt(fit1$s2.post) * fit1$stdev.unscaled[,ncol(design.mat)], p.value=fit1$p.value[,ncol(design.mat)])

stats <- stats[match(common.cpg, rownames(stats)),]
meth <- meth[match(common.cpg, rownames(meth)),]
stats <- cbind(stats, mani)

if(arraytype =="EPICv2"){
CHRs <- stats$CHR
Pos <- stats$MAPINFO
names(CHRs) <- rownames(stats)
names(Pos) <- rownames(stats)
}


if(arraytype == "EPICv1" | arraytype == "450K"){
CHRs <- stats$CHR_hg38
Pos <- (stats$Start_hg38) + 1
names(CHRs) <- rownames(stats)
names(Pos) <- rownames(stats)
}

dmrs <- dmrff::dmrff(estimate=stats$estimate, se=stats$se, p.value=stats$p.value, methylation=meth, chr=CHRs, pos=Pos, maxgap=probe_gap, verbose=T)

message("\n[MeCall]-[NOTICE] : dmrff algorithm Done.")

DMRobjs <- dmrs[which(dmrs$p.adjust <= adjPvalDmr & dmrs$n >= nProbes),]
DMRobjs <- DMRobjs[order(DMRobjs$p.value),]
rownames(DMRobjs) <- paste("DMR",1:nrow(DMRobjs),sep="_")
}

if(method == "DMRcate"){

message("\n[MeCall]-[NOTICE] : DMRcate algorithm Start.")
message("\n[MeCall]-[NOTICE] : This code based on DMR function in ChAMP.")
message("\n[MeCall]-[NOTICE] : If you want to know more detail information about code, you may want visit DMRcate R documentation.")

design.mat <- MeCall.MakeModelMatrix(pd = pd, covariate = covariate, interest = interest, PCA.result = NULL, Cell.comp = NULL, adjustment=FALSE)

if(mtype == "beta"){
meth <- BetatoM(meth)
}

if(is.null(probe_gap)){
probe_gap = 500
}


myannotation <- DMRcate::cpg.annotate("array", object=meth, what = "M", arraytype = arraytype, epicv2Remap = TRUE, analysis.type="differential", design=design.mat, coef=ncol(design.mat), fdr=dmp.fdr)

myannotationGlist <- as(myannotation@ranges,"GRangesList")
dmrcoutput <- DMRcate::dmrcate(myannotation, min.cpgs = nProbes, lambda=probe_gap, C=C, pcutoff ='fdr', betacutoff=NULL)

DMR <- as.data.frame(DMRcate::extractRanges(dmrcoutput, genome = "hg38"))
rownames(DMR) <- paste("DMR",1:nrow(DMR),sep="_")
DMRobjs <- DMR[which(DMR$HMFDR <= adjPvalDmr & DMR$n >= nProbes),]

}

if(method == "combP"){
message("\n[MeCall]-[NOTICE] : combined-pvalues algorithm is Python based algorithm.")
message("\n[MeCall]-[NOTICE] : Download combined-pvalues algorithm in here -> https://github.com/brentp/combined-pvalues.")
message("\n[MeCall]-[NOTICE] : This code will be return *.BED* format matrix object to run combined-pvalues.")
message("\n[MeCall]-[NOTICE] : BED file header => 1.chromosome 2.start 3.end 4.P-value.")

design.mat <- MeCall.MakeModelMatrix(pd = pd, covariate = covariate, interest = interest, PCA.result = NULL, Cell.comp = NULL, adjustment=FALSE)

DMPs <- MeCall.DMP(meth = meth, pd = pd, interest = interest, design = design.mat, cutoff.P=1, multi.P = 'BH', arraytype=arraytype)@Main.info

CHR <- as.character(DMPs[,"Chr"])
if(!str_detect(CHR[1],'chr')){
CHR <- sapply(CHR,function(x){
x <- paste0('chr',x)})
}

names(CHR) <- NULL


START <- (as.numeric((DMPs[,"Position"])))-1
END <- (as.numeric((DMPs[,"Position"])))
Pval <- as.numeric(DMPs[,"p"])

BED <- cbind(CHR,START,END,Pval)
BED <- as.data.frame(BED)

CHR.levels <- paste0("chr",c(1:22,"X","Y","M"))
BED$CHR <- factor(BED$CHR, levels = CHR.levels)
BED <- BED[c(order(BED$CHR,BED$START)),]
BED$CHR <- as.character(BED$CHR)

DMRobjs <- BED
message("\n[MeCall]-[NOTICE] : Please save return object as a .BED file.")
message("\n[MeCall]-[NOTICE] : To run combined-pvalues algorithm, user may refer the code below.")
message("\n[MeCall]-[NOTICE] : 1st step > $xx/acf.py -d 1:500:50 -c 4 [BED file].bed > /xx/acf.txt")
message("\n[MeCall]-[NOTICE] : 2nd step > $xx/slk.py --acf /xx/acf.txt -c 4 [BED file].bed > /xx/Pval_acf.bed")
message("\n[MeCall]-[NOTICE] : 3rd step > $xx/peaks.py --dist 750 --seed 0.05 /xx/Pval_acf.bed > /xx/region.bed")
message("\n[MeCall]-[NOTICE] : 4th step > $xx/region_p.py -p /xx/Pval_acf.bed -r /xx/region.bed -s 50 -c 5 > /xx/sig_region.bed")
message("\n[MeCall]-[NOTICE] : Please check github (https://github.com/brentp/combined-pvalues) for detail informations.")
}
return(DMRobjs)}

