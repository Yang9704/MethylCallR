#' Perfome functional enrichment analysis.
#'
#' @importFrom minfi getAnnotation
#' @importFrom missMethyl gometh
#'
#' @description
#' This function selects significant differentially methylated probes 
#' and tests which gene sets from gene ontology and KEGG database are 
#' significantly enriched based on these probes. Furthermore, the results 
#' will be restructured to allow users to view all outcomes at once.
#'
#' The `gometh()` function implemented in the `missMethyl` R package 
#' performs functional enrichment analysis using Wallenius' noncentral
#' hypergeometric distribution to minimize biases resulting from 
#' multiple relationships between CpG sites and genes.
#'
#' When conducting functional enrichment analysis using EPICv2, 
#' a matched manifest object is required due to its unique probe ID.
#' Since `MethylCallR` generates an `RGChannelSet` when loading IDAT files, 
#' the annotation object required by `gometh()` could be easily created.
#' The `MeCall.enrichment()` function will automatically perform this 
#' process by passing the return object of the `MeCall.enrichment()` function
#' to the `minfi.object.EPICv2` parameter.
#'
#' @param DMP.object A returned object from `MeCall.DMP()`. Please note that the results should include information for all probes (set to `cutoff.P = 1`).
#' @param DMP.cutoff The p-value cutoff for DMPs to be included in functional enrichment analysis.
#' @param DMP.multi Filter DMPs using the p-values adjusted with the selected multiple testing method ("none","BF", "BH").
#' @param Delta.Beta Filter DMPs using the specified delta-beta cutoff.
#' @param genomic.features Functional enrichment analysis is performed using only CpG sites located within the restricted functional region.
#' @param Enrich.cutoff The p-value cutoff for filtering significant descriptions.
#' @param Enrich.multi Filter significant descriptions using the p-values adjusted with the selected multiple testing method ("none","BF", "BH").
#' @param min.count Consider results significant only if the description includes the selected number or more genes.
#' @param minfi.object.EPICv2 To perform functional enrichment analysis using results from EPICv2, input the data object from `MeCall.ReadIdat()` function to obtain the list of genes associated with CpGs.
#'
#' @return A dataframe consisting of significantly enriched GO terms and KEGG pathways.
#'
#' @author Hyun-Ho Yang
#'
#' @seealso [missMethyl::gometh()]
#'
#' @references
#' Maksimovic, J. et al. (2021). Gene set enrichment analysis for genome-wide DNA methylation data.
#' Genome Biology. url{https://doi.org/10.1186/s13059-021-02388-x}
#'
#' @examples
#' \dontrun{
#' # Run functional enrichment analysis with default setting.
#' data.Enrichment <- MeCall.Enrichment(DMP.object = data.DMPs, DMP.cutoff = 0.05, DMP.multi = "BH",
#' Delta.Beta = c(-0.02,0.02), genomic.features = c("ALL"), Enrich.cutoff = 0.05, 
#' Enrich.multi = c("none"),min.count = 2 ,minfi.object.EPICv2 = NULL)
#'
#' # Run functional enrichment analysis using EPICv2 data.
#' data.Enrichment.EPICv2 <- MeCall.Enrichment(DMP.object = data.DMPs, DMP.cutoff = 0.05, 
#' DMP.multi = "BH", Delta.Beta = c(-0.02,0.02), genomic.features = c("ALL"), Enrich.cutoff = 0.05,
#' Enrich.multi = c("none"),min.count = 2 ,minfi.object.EPICv2 = data.Import)
#'}
#'
#' @export
MeCall.Enrichment <- function(DMP.object = data.DMPs, DMP.cutoff = 0.05, DMP.multi = "BH", Delta.Beta = c(-0.02,0.02), genomic.features = c("ALL"), Enrich.cutoff = 0.05, Enrich.multi = c("none","BF", "BH"),min.count = 2 ,minfi.object.EPICv2 = NULL){


if(DMP.object@accessory[1] == "450K"){
array.type = c("450K")
} else {
array.type = c("EPIC")
}

#### select Significant DMPs
Tr.Pvalue <- Translate.Pvalue(Nominal.P.list = DMP.object@Main.info$p, cutoff = DMP.cutoff, Original.method = DMP.multi, return.method = "BH")

subs.DMP <- subset(DMP.object@Main.info, adjusted.p < Tr.Pvalue)
Sig.DMPs <- subset(subs.DMP, deltabeta <= Delta.Beta[1] | deltabeta >= Delta.Beta[2])

SigCpGs <- rownames(Sig.DMPs)
allCpGs <- rownames(DMP.object@Main.info)

#### Set manifest for EPICv2
if(!is.null(minfi.object.EPICv2)){
annos <- getAnnotation(minfi.object.EPICv2$minfi.Set$rgSet)
com.idx <- intersect(rownames(DMP.object@Main.info), rownames(annos))
annos <- annos[com.idx, ]

annos$Name <- unlist(lapply(annos$Name, function(x){strsplit(x,split="_")[[1]][1]}))
SigCpGs <- unlist(lapply(SigCpGs, function(x){strsplit(x,split="_")[[1]][1]}))
allCpGs <- unlist(lapply(allCpGs, function(x){strsplit(x,split="_")[[1]][1]}))
}else{annos <- NULL}


#### perform GOmeth
GO.Results <- missMethyl::gometh(sig.cpg = SigCpGs, all.cpg = allCpGs, collection = c("GO"), array.type = array.type, plot.bias = FALSE, prior.prob = TRUE, anno = annos, equiv.cpg = TRUE, fract.counts = TRUE, genomic.features = genomic.features, sig.genes = TRUE)

KEGG.Results <- missMethyl::gometh(sig.cpg = SigCpGs, all.cpg = allCpGs, collection = c("KEGG"), array.type = array.type, plot.bias = FALSE, prior.prob = TRUE, anno = annos, equiv.cpg = TRUE, fract.counts = TRUE, genomic.features = genomic.features, sig.genes = TRUE)


#### Filtering enrichment analysis result
if(Enrich.multi == "BH"){
col.idx <- 5
}else{
col.idx <- 4
}

if(Enrich.multi == "BF"){
GO.cutoff <- Enrich.cutoff / nrow(GO.Results)
KE.cutoff <- Enrich.cutoff / nrow(KEGG.Results)
}else{
GO.cutoff <- Enrich.cutoff
KE.cutoff <- Enrich.cutoff
}

Sig.GO.Results <- subset(GO.Results, DE >= min.count)
Sig.GO.Results <- Sig.GO.Results[which(Sig.GO.Results[,col.idx+1] < GO.cutoff),]

Sig.KEGG.Results <- subset(KEGG.Results, DE >= min.count)
Sig.KEGG.Results <- Sig.KEGG.Results[which(Sig.KEGG.Results[,col.idx] < KE.cutoff),]

#### Return object
Results.obj <- Combine_Funtional_result(Sig.GO.Results,Sig.KEGG.Results)
return(Results.obj)}


#### Inner function #1 Combine_Funtional_result
Combine_Funtional_result <- function(GO,KEGG){

GO.counts <- table(GO$ONTOLOGY)
KEGG.counts <- nrow(KEGG)

if(is.na(GO.counts["BP"])){
BP.GO.ReCon <- data.frame(Category = c(), ID = c(), Description = c(), Gene_in_description = c(),Gene_Count = c(), p_value = c(), FDR = c(), Significant_gene = c())
}else{BP.GO <- subset(GO, ONTOLOGY == "BP")
BP.GO$ONTOLOGY <- "GO-Biological Process"
BP.GO <- BP.GO[order(BP.GO$P.DE),]

BP.GO.ReCon <- data.frame(Category = BP.GO$ONTOLOGY, ID = rownames(BP.GO), Description = BP.GO$TERM, Gene_in_description = BP.GO$N,Gene_Count = BP.GO$DE, p_value = BP.GO$P.DE, FDR = BP.GO$FDR, Significant_gene = BP.GO$SigGenesInSet)
}

if(is.na(GO.counts["CC"])){
CC.GO.ReCon <- data.frame(Category = c(), ID = c(), Description = c(), Gene_in_description = c(),Gene_Count = c(), p_value = c(), FDR = c(), Significant_gene = c())
}else{CC.GO <- subset(GO, ONTOLOGY == "CC")
CC.GO$ONTOLOGY <- "GO-Cellular Components"
CC.GO <- CC.GO[order(CC.GO$P.DE),]

CC.GO.ReCon <- data.frame(Category = CC.GO$ONTOLOGY, ID = rownames(CC.GO), Description = CC.GO$TERM, Gene_in_description = CC.GO$N,Gene_Count = CC.GO$DE, p_value = CC.GO$P.DE, FDR = CC.GO$FDR, Significant_gene = CC.GO$SigGenesInSet)
}


if(is.na(GO.counts["MF"])){
MF.GO.ReCon <- data.frame(Category = c(), ID = c(), Description = c(), Gene_in_description = c(),Gene_Count = c(), p_value = c(), FDR = c(), Significant_gene = c())
}else{MF.GO <- subset(GO, ONTOLOGY == "MF")
MF.GO$ONTOLOGY <- "GO-Molecular Function"
MF.GO <- MF.GO[order(MF.GO$P.DE),]

MF.GO.ReCon <- data.frame(Category = MF.GO$ONTOLOGY, ID = rownames(MF.GO), Description = MF.GO$TERM, Gene_in_description = MF.GO$N,Gene_Count = MF.GO$DE, p_value = MF.GO$P.DE, FDR = MF.GO$FDR, Significant_gene = MF.GO$SigGenesInSet)
}

if(KEGG.counts == 0){
KEGG.ReCon <- data.frame(Category = c(), ID = c(), Description = c(), Gene_in_description = c(),Gene_Count = c(), p_value = c(), FDR = c(), Significant_gene = c())
}else{KEGG <- KEGG[order(KEGG$P.DE),]

KEGG.ReCon <- data.frame(Category = "KEGG Pathway", ID = rownames(KEGG), Description = KEGG$Description, Gene_in_description = KEGG$N,Gene_Count = KEGG$DE, p_value = KEGG$P.DE, FDR = KEGG$FDR, Significant_gene = KEGG$SigGenesInSet)
}


Functional_object <- rbind(BP.GO.ReCon, CC.GO.ReCon, MF.GO.ReCon, KEGG.ReCon)
return(Functional_object)}

