#' Perform the copy number variation (CNV) analysis
#'
#' @import stats
#' @import DNAcopy
#' @importFrom preprocessCore normalize.quantiles
#'
#' @description
#' Copy number variation between two groups is identified using intensity
#' values. `MethylCallR` provides an automated function that performs 
#' data normalization, technical batch effect adjustment, and segmentation.
#' Finally, `MeCall.CNV()` returns a list object containing input files 
#' for the `GISTIC2.0` algorithm, which is used to identify significantly 
#' recurring focal alterations. By entering a directory in the 
#' `outdir` parameter, the function assumes that the user agrees to 
#' automatic file creation and saves the Segment and Marker files 
#' in the specified directory : 
#' * Segment file : `/outdir/MethylCallR_GISTIC_Segment_file.txt`
#' * Marker file : `/outdir/MethylCallR_GISTIC_Marker_file.txt`
#'
#'
#' @param data A return object from `MeCall.filtering()`.
#' @param interest The column name in the sample data file containing the phenotype of interest.
#' @param batch Before estimating CNVs, the given batch effect will be corrected using `ComBat` algoritm. The variable name must be included in the column names of the sample data file.
#' @param arraytype Version of the methylation microarray.
#' @param outdir Directory path to save the results in `GISTIC2.0` input format.
#' @param intercept logical. If set to `TRUE`, the intercept will also be considered. (for regression method)
#'
#' @return A list object containing the segment results per sample and the marker information.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [ChAMP::champ.CNA()], [DNAcopy::CNA()]
#'
#' @references
#'
#' Venkatraman, E.S. and Olshen, A.B. (2007). A faster circular binary segmentation algorithm 
#' for the analysis of array CGH data. Bioinformatics. 
#' /url{https://doi.org/10.1093/bioinformatics/btl646}
#'
#' For fitting method
#' Tian, Y. et al. (2017). ChAMP: updated methylation analysis pipeline for Illumina BeadChips.
#' Bioinformatics. /url{https://doi.org/10.1093/bioinformatics/btx513}
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.CNV] with default setting.
#' data.CNV <- MeCall.CNV(data = data.filtered, interest = "Sample_Group",
#' batch = c("Slide","Array"), arraytype=c("EPICv1"), outdir="./", intercept = TRUE)
#' }
#'
#' @export
MeCall.CNV <- function(data = data.filtered, interest = "Sample_Group", batch = c("Slide","Array"), arraytype="EPIC", outdir=NULL, intercept = TRUE){

message("\n[MeCall]-[NOTICE] : Copy Number Variation (CNV) analysis start.")

if(!is.null(batch)){
if(sum(batch %in% colnames(data$pd)) != length(batch)){
stop("\n[MeCall]-!!ERROR!! : There is an invalid batch effect name. Please check and try again.")
}}


Ori.pd <- data$pd
Intensity <- data$intensity

mani <- callmanifest(data$TAG)
CNA.mani.idx <- intersect(rownames(Intensity), rownames(mani))

message("\n[MeCall]-[NOTICE] : Fitting the signal intensity of the case group to the control group using the [mean] method. This method is written based on the source code of the [champ.CNA] function implemented in the [ChAMP R package].")
intsqn <- preprocessCore::normalize.quantiles(as.matrix(Intensity))

if(!is.null(batch)){
message("\n[MeCall]-[NOTICE] : ",length(batch)," batch effects have been detected. Known batch effects will be corrected using the ComBat algorithm.")
B.Inten.Obj <- MeCall.RemoveBatch(meth = intsqn, pd = Ori.pd, batches=batch, interest = interest, do.sva = FALSE)

Intens <- B.Inten.Obj$methylation.level
m.pd <- B.Inten.Obj$modified.pd
}else{
Intens <- intsqn
m.pd <- CheckMultipleGroup(Ori.pd, interest)$modified.pd
}

Con <- names(table(m.pd[,interest]))[1]
Case <- names(table(m.pd[,interest]))[2]

Case.Idx <- rownames(m.pd[which(m.pd[,interest] == Case),])
Control.Idx <- rownames(m.pd[which(m.pd[,interest] == Con),])

intsqnlog <- log2(Intens)
colnames(intsqnlog) <- colnames(Intensity)

Queries <- intsqnlog[,Case.Idx]
Controls <- intsqnlog[,Control.Idx]

CNA.Value <- apply(Queries,2,function(x){x - rowMeans(as.data.frame(Controls))})



rownames(CNA.Value) <- rownames(Queries)

CHR.levels <- paste0("chr",c(1:22,"X","Y","M"))

if(arraytype == "EPICv2"){
CHR.f <- factor(mani[CNA.mani.idx,"CHR"], levels = CHR.levels)
MAPINFO <- mani[CNA.mani.idx,"MAPINFO"]
} else {
CHR.f <- factor(mani[CNA.mani.idx,"CHR_hg38"], levels = CHR.levels)
MAPINFO <- mani[CNA.mani.idx,"Start_hg38"] + 1
}

CHR <- as.numeric(CHR.f)
names(CHR) <- CNA.mani.idx
names(MAPINFO) <- CNA.mani.idx

message("\n[MeCall]-[NOTICE] : Finding copy number variation using circular binary sementation algorithm.")

CNA.Result <- list()
for (i in 1:ncol(CNA.Value)){
CNA.object <- DNAcopy::CNA(cbind(CNA.Value[,i]), chrom = CHR, maploc = MAPINFO ,data.type = "logratio", sampleid = colnames(CNA.Value)[i])

smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)

segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)

seg <- segment.smoothed.CNA.object$output
CNA.Result[[colnames(CNA.Value)[i]]] <- seg}

Total.Seg.data <- do.call(rbind,CNA.Result)

Markers <- data.frame(CHR, MAPINFO)
rownames(Markers) <- CNA.mani.idx
Markers <- Markers[order(Markers[,1], Markers[,2]),]

message("\n[MeCall]-[NOTICE] : Segmentation Finished. MethylCallR will return a list object containing segment data and marker information.")

CNA.Return.Obj <- list(Segment.data = Total.Seg.data, Marker.data = Markers)

if(!is.null(outdir)){

if(str_sub(outdir, start= -1) != "/"){
outdir <- paste0(outdir,"/")
}

if(!dir.exists(outdir)){
message("\n[MeCall]-[NOTICE] : Invalid directory path. MethylCallR will return result object without auto save.")
return(CNA.Return.Obj)
}

message("\n[MeCall]-[NOTICE] : MethylCallR will save the return object as a text file in the following directory : ", outdir)

Seg.Paths <- paste0(outdir,"MethylCallR_GISTIC_Segment_file.txt")
Mark.Paths <- paste0(outdir,"MethylCallR_GISTIC_Marker_file.txt")

write.table(CNA.Return.Obj$Segment.data, file = Seg.Paths, sep = "\t", quote = FALSE,row.names=FALSE, col.names = FALSE)
write.table(CNA.Return.Obj$Marker.data, file = Mark.Paths, sep = "\t", quote = FALSE,row.names=FALSE, col.names = FALSE)
}

return(CNA.Return.Obj)}
