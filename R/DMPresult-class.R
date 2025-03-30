DMPresultlist <- methods::setClass("DMPresultlist", slots=c(Main.info = "data.frame",Other.info = "data.frame",accessory = "character"))


Set.DMPresultlist <-function(O,arraytype = "EPICv1",adjmethod,cutoff,pd=NULL,group = NULL){
if(arraytype == "EPICv1" | arraytype == "450K"){
f.region <- unlist(lapply(O$UCSC_RefGene_Group,function(x){t <- str_split(x,";")[[1]][1]}))
g.name <- unlist(lapply(O$UCSC_RefGene_Name,function(x){t <- str_split(x,";")[[1]][1]}))
CGIs <- unlist(lapply(O$Relation_to_UCSC_CpG_Island,function(x){t <- str_split(x,";")[[1]][1]}))

maininfo = data.frame(CpGid=O$Name, Chr=O$CHR_hg38,Position = (O$End_hg38 - 1),strand = O$Strand_hg38,gene = g.name,functional.region = f.region,CGI.region = CGIs,deltabeta=O$logFC, t=O$t, p=O$P.Value, adjusted.p=O$adj.P.Val,AveMeth = O$AveExpr, B = O$B,SNP_ID=O$SNP_ID, SNP_DISTANCE=O$SNP_DISTANCE)

Others = data.frame(CpGid=O$Name, AddressA_ID=O$AddressA_ID, AddressB_ID = O$AddressB_ID, Infinium_Design_Type = O$Infinium_Design_Type, Color_Channel=O$Color_Channel,Chr_hg37=O$CHR, Position_hg37 = O$MAPINFO, strand_hg37=O$Strand, UCSC_RefGene_Name=O$UCSC_RefGene_Name, UCSC_RefGene_Group=O$UCSC_RefGene_Group, UCSC_CpG_Islands_Name=O$UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island=O$Relation_to_UCSC_CpG_Island, Methyl27_Loci=O$Methyl27_Loci, Methyl450_Loci=O$Methyl450_Loci, SNP_MinorAlleleFrequency=O$SNP_MinorAlleleFrequency, MFG_Change_Flagged=O$MFG_Change_Flagged)

}else{
f.region <- unlist(lapply(O$UCSC_RefGene_Group,function(x){t <- str_split(x,";")[[1]][1]}))
g.name <- unlist(lapply(O$UCSC_RefGene_Name,function(x){t <- str_split(x,";")[[1]][1]}))
CGIs <- unlist(lapply(O$Relation_to_UCSC_CpG_Island,function(x){t <- str_split(x,";")[[1]][1]}))

maininfo = data.frame(CpGid=O$Name, Chr=O$CHR,Position = O$MAPINFO,strand = O$Strand_FR,gene = g.name,functional.region = f.region,CGI.region = CGIs, deltabeta=O$logFC, t=O$t, p=O$P.Value, adjusted.p=O$adj.P.Val,AveMeth = O$AveExpr, B = O$B,SNP_ID=O$SNP_ID, SNP_DISTANCE=O$SNP_DISTANCE)


Others = data.frame(CpGid=O$Name, AddressA_ID=O$AddressA_ID, AddressB_ID = O$AddressB_ID, Strand_FR = O$Strand_FR, Infinium_Design_Type = O$Infinium_Design_Type, Color_Channel=O$Color_Channel, UCSC_RefGene_Name=O$UCSC_RefGene_Name,  UCSC_RefGene_Group=O$UCSC_RefGene_Group, UCSC_CpG_Islands_Name=O$UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island=O$Relation_to_UCSC_CpG_Island, SNP_MinorAlleleFrequency=O$SNP_MinorAlleleFrequency)
}

rownames(maininfo) <- rownames(O)
rownames(Others) <- rownames(O)

A = c(array = arraytype, Multi.correction = paste0(adjmethod,"-method"),Cutoff = cutoff,Group = group)
Object <- DMPresultlist(Main.info =maininfo, Other.info = Others, accessory = A)
return(Object)
}

