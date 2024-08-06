#' Generate various plots based on ggplot2 R package.
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import ggnewscale
#' @importFrom cowplot plot_grid
#'
#' @description
#' This function generates the selected plot based on the input data 
#' and returns a ggplot2 object.
#'
#' For a `Manhattan` plot, users may want to submit a named list object 
#' to the `highlight` parameter to simultaneously highlight probes 
#' from multiple subgroups. Through this, users can freely designate 
#' three or more levels of highlighting. To label probes with their 
#' associated gene names, users can simply input a vector composed 
#' of probe IDs into the `gene.label` argument.
#'
#'
#' @param data The `DimReductionlist` object from `MeCall.PCA()` for `scree` and `pca` plot. The `DMPresultlist` object from `MeCall.DMP()` for `Manhattan`,`Volcano`, `functional.region`, and `CGI.region` plot.
#' @param type The type of plot to create.
#' @param pd A data frame of modified sample data file from `MeCall.RemoveBatch()`. (pca plot)
#' @param batch.name Highlight samples with different colors of dots based on the given variable. (pca plot) 
#' @param PC.x The principal component to display on the x-axis. (pca plot)
#' @param PC.y The principal component to display on the y-axis.(pca plot)
#' @param Pval.cut The cutoff for the p-value regardless of the multiple testing correction method. (Manhattan, Volcano, region plot)
#' @param gene.label A vector consisting of probe IDs for which the user wants to indicate associated gene annotations. (Manhattan, Volcano plot)
#' @param highlight A list object with named slots. Each slot contains a vector consisting of probe IDs that the user wants to highlight. (Manhattan plot)
#' @param highlight.color Highlight probes in each slot with dots of different specified colors. (Manhattan plot)
#' @param delta.cut The threshold for significant delta-beta values. (Volcano, region plot)
#'
#' @return A ggplot2 object
#'
#' @author Hyun-Ho Yang
#'
#' @examples
#' \dontrun{
#' # Run [MeCall.Plot] to generate scree plot.
#' Scree_plot <- MeCall.Plot(data = data.PCA, type= c("scree"))
#'
#' # Run [MeCall.Plot] to generate pca plot.
#' pca_plot <- MeCall.Plot(data = data.PCA, type= c("pca"), batch.name="Sample_Group", 
#' PC.x = 1, PC.y = 2)
#'
#' # Run [MeCall.Plot] to generate Manhattan plot.
#' Highlight.list = list(H1 = c("cg00000001","cg00000002"), H2 = c("cg00000003","cg00000004"), 
#' H3 = c("cg00000005","cg00000006"))
#' Manhattan_plot <- MeCall.Plot(data = data.DMPs, type= c("Manhattan"), Pval.cut = 0.05, 
#' highlight = Highlight.list, highlight.color = c("red","green"),gene.label = c("cg00000001",
#' "cg00000002","cg00000003"))
#'
#' # Run [MeCall.Plot] to generate Volcano plot.
#' Volcano_plot <- MeCall.Plot(data = data.DMPs, type= c("Volcano"), Pval.cut = 0.05, 
#' delta.cut = c(-0.02,0.02), gene.label = c("cg00000001","cg00000002"))
#'
#' # Run [MeCall.Plot] to generate functional/CGI region plot.
#' Region_plot <- MeCall.Plot(data = data.DMPs, type= c("functional.region"), Pval.cut = 0.05, 
#' delta.cut = c(-0.02,0.02))
#'}
#'
#' @export
MeCall.Plot <- function(data = data, type=c("scree","pca","Manhattan","Volcano","functional.region","CGI.region"), pd=pd, batch.name=NULL, PC.x=1, PC.y=2, delta.cut = c(0,0), Pval.cut = 0.05, highlight=NULL,highlight.color = NULL,gene.label = NULL){

if(is(data) == "DMPresultlist"){
man_data <- data@Main.info
cutoff <- Pval.cut
}


if(type == c("scree")){
message("\n[MeCall]-[notice] : Generating scree plot.")
if(is(data) != "DimReductionlist"){
stop("\n[MeCall]-!!ERROR!! : Input data format is not a [MethylCallR::DimReductionlist]. Please check again")}
prop.var <- data@cum.var
df <- data.frame(components = colnames(data@Components), proportion = prop.var)
f.levels <- paste0("PC",1:length(df$components))
df$components <- factor(df$components, levels = f.levels)
out.plot <- ggplot(data=df,aes(x=components, y= proportion,group = 1)) + geom_point() + geom_text(label=round(df$proportion,digits=3),hjust=-0.3, vjust=-0.5) + geom_line(linewidth=0.5, color = "Dark gray", linetype = "longdash") + theme_gray()
}

if(type == c("pca")){
if(is(data) != "DimReductionlist"){
stop("\n[MeCall]-!!ERROR!! : Input data format is not a [MethylCallR::DimReductionlist]. Please check again")}
message("\n[MeCall]-[notice] : Generating pca plot.")
df <- as.data.frame(data@Components)
PC.x <- colnames(df)[PC.x]
PC.y <- colnames(df)[PC.y]
if(length(batch.name) > 1){
stop("\n[MeCall]-!!ERROR!! : Only 1 batch name is allowed.")

}else if(!is.null(batch.name)){
if(!(batch.name %in% colnames(pd))){
stop("\n[MeCall]-!!ERROR!! : There is no [",batch.name,"] in pd file.")

}else{df <- cbind(df, pd[,batch.name])

out.plot <- ggplot(data=df, aes(x=.data[[PC.x]],y=.data[[PC.y]], color = batch.name)) + geom_point() + geom_hline(yintercept = 0, linetype='longdash',linewidth=0.3,color="Dark gray") + geom_vline(xintercept = 0,linetype='longdash',linewidth=0.3,color="Dark gray") + xlab(PC.x) + ylab(PC.y)+ labs(color = batch.name)}
}}


if(type == "Volcano"){
message("\n[MeCall]-[notice] : Volcano plot is selected.")
if(is(data) != "DMPresultlist"){
stop("\n[MeCall]-!!ERROR!! : Input data format is not a [MethylCallR::DMPresultlist]. Please check again")}
n.cut <- tail(subset(man_data, adjusted.p <= cutoff),1)$p
man_data[,"status"] <- "not significant"
man_data[,"status"][man_data$deltabeta < delta.cut[1] & -log10(man_data$p) >= -log10(n.cut)] <- "Hypo-methylated"
man_data[,"status"][man_data$deltabeta > delta.cut[2] & -log10(man_data$p) >= -log10(n.cut)] <- "Hyper-methylated"
mycolors <- c("orangered3","deepskyblue3","dark grey")
names(mycolors) <- c("Hypo-methylated","Hyper-methylated","not significant")
numbers <- c(length(which(man_data[,"status"] == "Hyper-methylated")), length(which(man_data[,"status"] == "Hypo-methylated")), length(which(man_data[,"status"] == "not significant")))
mylabels <- paste0(names(mycolors)," (",numbers,")")


#### Volcano plot BACK BONE
out.plot <- ggplot(data = man_data, aes(x=deltabeta, y=-log10(p),col = status)) + geom_point()+ ylab("-log10(P)") + geom_hline(yintercept = -log10(n.cut), linetype=10, col="red", size=1.5) + geom_vline(xintercept=delta.cut[1],linetype=8,col="grey30",size=1.5)+ geom_vline(xintercept=delta.cut[2],linetype=8,col="grey30",size=1.5)+ scale_colour_manual("Number of CpG sites", values=mycolors,breaks = c("Hypo-methylated", "Hyper-methylated", "not significant"), labels = mylabels)

#### Volcano plot set theme
out.plot <- out.plot + theme_minimal() + theme(legend.title=element_text(size=16), legend.text=element_text(size=14), legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid',linewidth = 1), axis.title.x = element_text(size=18,face='bold'), axis.text.x = element_text(size=17), axis.line.x = element_line(linewidth=0.7), axis.title.y = element_text(size=18,face='bold'), axis.text.y = element_text(size=17),axis.line.y = element_line(linewidth=0.7),panel.grid.major.x = element_line(colour = "grey75",linewidth=0.8), panel.grid.minor.x = element_line(colour = "grey75", linewidth=0.3), panel.grid.major.y = element_line(colour = "grey75",linewidth=0.8), panel.grid.minor.y = element_line(colour = "grey75", linewidth=0.3), legend.box="vertical", legend.margin=margin(5,5,5,5)) + guides(fill = guide_legend(override.aes = aes(label = "")))

#### Volcano plot gene labeling
if(!is.null(gene.label)){
label_don <- subset(man_data, CpGid %in% gene.label)
out.plot <- out.plot + ggrepel::geom_label_repel(data = label_don, aes(label = gene), size = 6, fill= "white",colour = "black", box.padding = 1, max.overlaps= Inf, segment.color = "black")}
}




if(type == "Manhattan"){
message("\n[MeCall]-[notice] : Manhattan plot is selected.")
n.cut <- tail(subset(man_data, adjusted.p <= cutoff),1)$p
f.levels <- paste0("chr",c(1:22,"X","Y"))
f.levels <- intersect(f.levels,unique(man_data$Chr))
Chr_color <- rep(c("moccasin","plum3"), 15)[1:length(f.levels)]
man_data$Chr <- factor(man_data$Chr, levels = f.levels)
don <- man_data %>% dplyr::group_by(Chr) %>% dplyr::summarise(chr_len=max(Position)) %>% dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% dplyr::select(-chr_len) %>% dplyr::left_join(man_data, ., by=c("Chr"="Chr")) %>% dplyr::arrange(Chr, Position) %>% dplyr::mutate(BPcum=Position+tot)

axisdf = don %>% dplyr::group_by(Chr) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf$Chr <- as.numeric(axisdf$Chr)

#### Manhattan plot BACK BONE
out.plot <- ggplot(don,aes(x=BPcum, y=-log10(p),color=Chr)) +geom_point(alpha=0.8, size=1.3) +scale_color_manual("Chr",values = Chr_color,guide="none") +scale_x_continuous(labels = axisdf$Chr, breaks= axisdf$center) +xlab("Chromosome") +scale_y_continuous(limits = c(0, ceiling(-log10(min(don$p))))) +geom_hline(yintercept=-log10(n.cut), linetype='dashed', color='red', size=1)

#### Manhattan plot Highlight 
if(!is.null(highlight)){
don$highlight <- "-"
for (x in 1:length(highlight)) {don[don$CpGid %in% highlight[[x]],]$highlight <- names(highlight)[x]
}
high_don <- subset(don, !highlight == "-")
counts <- unlist(lapply(names(highlight), function(x) {nrow(subset(don, highlight == x))}))
lbs <- paste0(names(highlight)," (",counts,")")
names(highlight.color) <- names(highlight)

out.plot <- out.plot + ggnewscale::new_scale_color()  +geom_point(data =high_don, aes(x=BPcum, y=-log10(p),color =highlight)) +scale_color_manual(name= NULL, labels = lbs, values = highlight.color,breaks = names(highlight.color))
}

#### Manhattan plot set theme
out.plot <- out.plot + theme_bw() + theme(axis.title.x = element_text(size=18,face='bold'), axis.text.x = element_text(size=16),axis.line.x = element_line(linewidth=0.7), axis.title.y = element_text(size=18,face='bold'), axis.text.y = element_text(size=17),axis.line.y = element_line(linewidth=0.7),legend.title=element_text(size=18), legend.text=element_text(size=14), legend.position="right",legend.background = element_rect(fill = "White", color = 'Black',linewidth=1),panel.border = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(colour = "grey75",linewidth=0.8), panel.grid.minor.y = element_line(colour = "grey75", linewidth=0.5),legend.box="vertical",legend.margin=margin(5,5,5,5))

### Manhattan plot Gene Labeling
if(!is.null(gene.label)){
label_don <- subset(don, CpGid %in% gene.label)
out.plot <- out.plot + ggrepel::geom_label_repel(data = label_don, aes(x = BPcum, y=-log10(p), label = gene), size = 6, fill= "white", colour ="black", box.padding = 1,max.overlaps= Inf)
}
}
if(type == "functional.region"){
message("\n[MeCall]-[notice] : Polar bar plots will be generated to represent the distribution of functional regions.")
man_data$featuref <- factor(man_data$functional.region,levels=c("TSS1500","TSS200","5'UTR","1stExon","ExonBnd","Body","3'UTR","IGR"))
hyper <- subset(man_data, adjusted.p <= cutoff & deltabeta > delta.cut[2])
hypo <- subset(man_data, adjusted.p <= cutoff & deltabeta < delta.cut[1])


hypo <- as.data.frame(table(hypo$featuref))
hyper <- as.data.frame(table(hyper$featuref))
colnames(hypo) <- c("feature","count")
colnames(hyper) <- c("feature","count")


hypo_p <- ggplot(hypo,aes(x=feature,y=count)) + geom_bar(stat="identity") +ggtitle("Distribution of functional region of hypo-methylated CpGs")+geom_col(aes(x=feature,y=count,fill=count),position = "dodge2",show.legend = TRUE,alpha =0.9) +scale_fill_gradientn("Total count of each feature",colours = c("#6C5B7B","#C06C84","#F67280","#F8B195"))+guides(fill = guide_colorsteps(barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5)) +theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(color = "gray12", size = 12),legend.position = "bottom")+coord_polar()

hyper_p <- ggplot(hyper,aes(x=feature,y=count)) + geom_bar(stat="identity") +ggtitle("Distribution of functional region of hyper-methylated CpGs")+geom_col(aes(x=feature,y=count,fill=count),position = "dodge2",show.legend = TRUE,alpha =0.9) +scale_fill_gradientn("Total count of each feature",colours = c("#6C5B7B","#C06C84","#F67280","#F8B195"))+guides(fill = guide_colorsteps(barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5)) + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(color = "gray12", size = 12),legend.position = "bottom")+coord_polar()


out.plot <- cowplot::plot_grid(hypo_p,hyper_p, ncol=2)
}
if(type == "CGI.region"){
message("\n[MeCall]-[notice] : Pie plots will be generated to represent the distribution of CGI regions.")
man_data$featuref <-man_data$CGI.region
man_data$featuref[man_data$featuref == ""] <- "Opensea"
man_data$featuref[man_data$featuref == "N_Shelf" | man_data$featuref == "S_Shelf"] <- "Shelf"
man_data$featuref[man_data$featuref == "N_Shore" | man_data$featuref == "S_Shore"] <- "Shore"
hypo <- subset(man_data, adjusted.p <= cutoff & deltabeta < delta.cut[1])
hyper <- subset(man_data, adjusted.p <= cutoff & deltabeta > delta.cut[2])
hypo <- as.data.frame(table(hypo$featuref))
hyper <- as.data.frame(table(hyper$featuref))
colnames(hypo) <- c("cgi","count")
colnames(hyper) <- c("cgi","count")

hypo_p <- ggplot(hypo,aes(x="",y=count,fill=cgi))+ geom_bar(width=1,stat="identity")+ coord_polar("y",start=0)+ scale_fill_manual(values=c("darkorange1","skyblue1","Wheat1","tan1"))+ theme_bw() + theme(axis.text.x=element_blank())+ geom_text(aes(label = scales::percent(count/sum(count))),position=position_stack(vjust=0.5),size=5)+guides(fill = guide_legend(title = "")) + labs(title="Hypo-methylated CpG location") + theme_void()+ theme(plot.title = element_text(hjust=0.5),legend.position="bottom",axis.title = element_blank(),axis.ticks = element_blank())

hyper_p <- ggplot(hyper,aes(x="",y=count,fill=cgi))+ geom_bar(width=1,stat="identity")+ coord_polar("y",start=0)+ scale_fill_manual(values=c("darkorange1","skyblue1","Wheat1","tan1"))+ theme_bw() + theme(axis.text.x=element_blank())+ geom_text(aes(label = scales::percent(count/sum(count))),position=position_stack(vjust=0.5),size=5)+guides(fill = guide_legend(title = "")) + labs(title="Hyper-methylated CpG location") + theme_void()+ theme(plot.title = element_text(hjust=0.5),legend.position="bottom",axis.title = element_blank(),axis.ticks = element_blank())

out.plot <- cowplot::plot_grid(hypo_p,hyper_p, ncol=2)
}

if(!type %in% c("scree","pca","Manhattan","Volcano","functional.region","CGI.region")){
stop("\n[MeCall]-!!ERROR!! : Please check plot type parameter again.")
}

return(out.plot)}


