#' Perform differential methylation analysis to identify differentially methylated probes (DMPs).
#' 
#' @import limma
#' @importFrom methods as new
#' 
#' @description
#' Identify differentially methylated probes between two groups using 
#' preprocessed methylation level data.
#'
#' Only probes with p-values adjusted for multiple testing using 
#' the `multi.P` method and less than or equal to `cutoff.P` are 
#' included in the result object.
#'
#' By specifying the appropriate `arraytype` parameter (`450L`, `EPICv1`,
#' `EPICv2`), accurate positional information and characteristics for 
#' each probe can be provided. MethylCallR provides positional information
#' based on the GRCh38/hg38 human reference as a default, while 
#' the GRCh37/hg19 human reference is included as additional information.
#'
#' @param meth A matrix of methylation level.
#' @param pd A data frame of modified sample data file from `MeCall.RemoveBatch()`.
#' @param interest The column names in the sample data file containing the phenotype of interest.
#' @param design A design matrix from `MeCall.MakeModelMatrix()`
#' @param cutoff.P Probes with p-value above the specified threshold, as determined by the method specified in the `multi.P` parameter, are excluded from the results. If set to `1`, you can obtain results for all probes.
#' @param multi.P Set the method for multiple testing correction ("none","BH","BF").
#' @param arraytype Obtain annotation information for the given array type.
#'
#' @return "DMPresultlist" S4 class object consisted of :
#' \itemize{Main.info} {The data frame contains statistical values and major information for each probe. This information includes chromosome, position, functional region, CGI region, and gene name. For a detailed information of the statistical values, refer to the [limma::topTable()] function.}
#' \itemize{Other.info} {The data frame contains supplementary details from the Illumina manifest, excluding the major information.}
#' \itemize{accessory} {The vector includes the array type, the method used for multiple testing correction, the p-value threshold, and the names of the comparison groups.}
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [limma::topTable()]
#'
#' @references
#' Ritchie, M.E. et al. (2015). limma powers differential expression analyses for RNA-sequencing 
#' and microarray studies. Nucleic Acids Research. /url{https://doi.org/10.1093/nar/gkv007}
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.DMP] with default setting.
#' data.DMPs <- MeCall.DMP(meth = data.CellComp$adjusted.meth, 
#' pd = data.rmBatch$modified.pd, interest = "Sample_Group", design = data.design, 
#' cutoff.P=1, multi.P = 'BH', arraytype="EPICv1")
#' }
#'
#' @export
MeCall.DMP <- function(meth = data.CellComp$adjusted.meth, pd = data.rmBatch$modified.pd, interest = "Sample_Group", design = data.design, cutoff.P=1, multi.P = 'BH', arraytype=c("450k","EPICv1","EPICv2")){
message("\n[MeCall]-[notice] : Find DMPs with input data.")
adjs <- c("none","BH","BF")
if(!multi.P %in% adjs){
stop("\n[MeCall]-!!ERROR!! : You set wrong correction method. Choose correctly. : 'none' or 'BF'(bonferroni) or 'BH(benjamini-hochberg)'")
}
probe.features <- callmanifest(arraytype)
m.pd <- pd$modified.pd
CaCo <- pd$groups
fit <- limma::lmFit(meth,design)
fit.e <- limma::eBayes(fit)
IV <- colnames(fit$coefficients)[ncol(fit$coefficients)]
if(multi.P == "none" | multi.P =="BH"){
message("\n[MeCall]-[notice] : ",multi.P," adjustment method is selected.")
DMP <- limma::topTable(fit.e, coef=IV, adjust.method=multi.P, sort.by="P", number=Inf)
} else{DMP <- limma::topTable(fit.e, coef=IV, adjust.method="none", sort.by="P", number=Inf)
message("\n[MeCall]-[notice] : Bonferroni adjustment is selected.")
DMP$adj.P.Val <- p.adjust(DMP$P.Value,method = "bonferroni")
}
com.idx <- intersect(rownames(DMP),rownames(probe.features))
DMPo<- data.frame(DMP[com.idx,],probe.features[com.idx,])
DMPo <- subset(DMPo, adj.P.Val <= cutoff.P)
DMP <- Set.DMPresultlist(O = DMPo, arraytype = arraytype, adjmethod =multi.P,cutoff = cutoff.P, pd = m.pd,group = CaCo)
return(DMP)}
