################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){
  install.packages("ggplot2")}
if("rrcov" %in% rownames(installed.packages()) == FALSE){
  install.packages("rrcov")}
#if("MSstatsTMT" %in% rownames(installed.packages()) == FALSE){
#  BiocManager::install("MSstatsTMT")}
if("limma" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("limma")}
if("edgeR" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("edgeR")}
if("readxl" %in% rownames(installed.packages()) == FALSE){
  install.packages("readxl")}
if("bcv" %in% rownames(installed.packages()) == FALSE){
  install.packages("bcv")}
if("FactoMineR" %in% rownames(installed.packages()) == FALSE){
  install.packages("FactoMineR")}
library(FactoMineR)
if("factoextra" %in% rownames(installed.packages()) == FALSE){
  install.packages("factoextra")}
library(factoextra)
if("pcaMethods" %in% rownames(installed.packages()) == FALSE){
  install.packages("pcaMethods")}
library(pcaMethods)
library(factoextra)
if("ggplot2" %in% rownames(installed.packages()) == FALSE){
  install.packages("ggplot2")}
library(ggplot2)
if("ggrepel" %in% rownames(installed.packages()) == FALSE){
  install.packages("ggrepel")}
library(ggrepel)
if("biomaRt" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("biomaRt")}

################################################################################
#Functions.
################################################################################
require(ggplot2)
brandontheme=theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  plot.title = element_text(color="black", size=14,
                            face="plain",hjust = 0.5),
  axis.text = element_text(color="black", size=12,
                           face="plain",hjust = 0.7),
  axis.title = element_text(color="black", size=12, face="plain",
                            hjust = 0.5))
norm_dist_plot= function(x_values,channel=NA){
  ggplot(data = data.frame(x = x_values), aes(x)) + #not necessary
    stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) + ylab("") +
    scale_y_continuous(breaks = NULL) +
    brandontheme +
    scale_x_continuous(expand = c(0, 0.01), breaks =  seq(from=-8,
                                                          to= max(x_values)*1.01,
                                                          by = 1)) +
    labs(y="Density of distribution",x="Log2 intensity",
         title = paste0("Distribution of intensities for channel", channel)
    ) + 
    brandontheme
}
pseudo_counts= function(a_list){
  list_no_invalid= na.omit(a_list[a_list != 0])
  list_no_invalid= list_no_invalid[!is.infinite(list_no_invalid)]
  a_list[a_list == 0 | is.na(a_list)]= min(list_no_invalid)
  a_list
}
require(limma)
log_frame= function(treated_matrix,CTRL_matrix,genes){
  output= cbind(CTRL_matrix,treated_matrix)
  output= output[!duplicated(genes),]
  rownames(output)=genes[!duplicated(genes)]
  output= log2(output)
  design= cbind(overlap=1,treated=c(0,0,0,1,1,1))
  data= output
  data= as.matrix(data)
  fit= lmFit(data, design)
  fit= eBayes(fit)
  topTable(fit, coef=2,number = nrow(data))
}
valid_cols= function(a_list, thershold= 0.7){
  sum(!is.na(a_list)) > round(length(a_list), 0)*thershold
}
################################################################################
#Read excel
################################################################################
TMT= readxl::read_xlsx(
  "~/Documents/Rhys_Proteomics/Data/Rhys Morgan 120122 6Plex TMT_KG1 RM.xlsx")
#bool= as.logical(apply(TMT[,21:26],1,valid_cols,thershold= 0.7)) #using columns 15 to 20 scaled intensities- i.e. normalised already.
#TMT= TMT[which(bool), ] #removing proteins with fewer than 70% of samples present (actually 66.67% due to rounding)
#imputed based on 10.1038/s41598-017-19120-0 assuming random NAs.
#bool= as.logical(apply(TMT[,21:26],1,valid_cols,thershold= 0.7)) #using columns 21 to 26 normalised intensities- i.e. normalised already.
TMT= TMT[!as.logical(apply(TMT[,21:26], 1, filter_missing_vals)), ] #getting rid of peptides lacking an intensity in 3 replicates

################################################################################
#PCA graph
################################################################################
PCA_df= t(TMT_matrix)
PCA_df= na.omit(PCA_df)
rownames(PCA_df) = gsub("_"," ",rownames(PCA_df))
res.pca <- PCA(PCA_df, ncp = 2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
#
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map of channels- 1-3 KOs")
################################################################################
#Log frame and Volcano plot set up
################################################################################
log_frame_output= log_frame(treated_matrix= TMT[,24:26], 
                            CTRL_matrix= TMT[,21:23], TMT$Accession)
log_frame_output= as.data.frame(log_frame_output)
log_frame_output$Expression="UNCHANGED"
log_frame_output$Expression[log_frame_output$logFC >= log2(2) &
                              log_frame_output$P.Value <= 0.05]= "UP"
log_frame_output$Expression[log_frame_output$logFC <= log2(0.5) &
                              log_frame_output$P.Value <= 0.05] <- "DOWN"
length(unique(log_frame_output$Expression))== 3 # Checks if all three data types present.
log_frame_output$label <- "" #Prevents errors.
mycolours <- c("blue", "red", "black") #Changing colours for plot.
names(mycolours) <- c("DOWN", "UP", "UNCHANGED")
mycolours <- c("blue", "black", "red") #Changing colours for plot.
################################################################################
#Getting list of plasma membrane proteins
################################################################################
require(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
plasma_membrane_list_full= getBM(attributes= c("hgnc_symbol","uniprot_gn_id"), 
                                 filters=c("go"),
                                 values=list("GO:0005886"),
                                 mart=ensembl)
plasma_membrane_list= plasma_membrane_list_full$uniprot_gn_id
log_frame_output$membrane= F
log_frame_output$membrane[rownames(log_frame_output) %in% 
                            plasma_membrane_list]= T
log_frame_output= log_frame_output[order(abs(log_frame_output$logFC),
                                         decreasing =T),]

################################################################################
#Adding labels to plot
################################################################################
log_frame_output_membrane_only= log_frame_output[which(
  log_frame_output$membrane ==T & log_frame_output$P.Value <0.05), ]
log_frame_output_membrane_only= log_frame_output_membrane_only[abs(log_frame_output_membrane_only$logFC) > 1,]
top_membrane_proteins= head(rownames(
  log_frame_output_membrane_only[order(log_frame_output_membrane_only$logFC),])
  ,20)
top_membrane_proteins= c(top_membrane_proteins, head(rownames(
  log_frame_output_membrane_only[order(log_frame_output_membrane_only$logFC,
                                       decreasing= T),])
  ,20))
top_membrane_proteins_in_english= plasma_membrane_list #frame of proteins names in uniprot format
log_frame_output[top_membrane_proteins,"label"]=
  as.character(lapply(top_membrane_proteins,function(x){
    y= plasma_membrane_list_full$hgnc_symbol[plasma_membrane_list_full$uniprot_gn_id == x]
    y[1]
  }))

################################################################################
#Plotting the data
################################################################################
ggplot(data=log_frame_output, aes(x=logFC, y=-log10(P.Value),col=Expression,label=label)) +
  geom_point(size=2) +
  brandontheme +
  #geom_vline(xintercept=c(-1, 1), col="black",cex=1.1) + #Addes axes to identify points of Expression.
  #geom_hline(yintercept=-log10(0.05), col="black",cex=1.1) +
  scale_colour_manual(values = mycolours) + #Changes points' colours to mycolors scheme.
  labs(y= "-log10(p value)",x="Fold change",title ="Volcano plot for CHIR treated KG1 with top 20 upregulated and downregulated proteins highlighted") +
  scale_x_continuous(expand = c(0, 1), breaks =  seq(from= round(min(log_frame_output$logFC)*1.2,0),
                                                     to= round(max(log_frame_output$logFC)*1.2,0),
                                                     by = 1)) +
  scale_y_continuous(expand = c(0, 1), breaks =  seq(from= 0,
                                                     to= round(max(-log10(log_frame_output$P.Value))*1.2,0),
                                                     by = 2)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf)
ggsave("~/sRNA_KG1_volcano.tiff",width=5000,height= 2000,units = "px")

################################################################################
#Heatmap
################################################################################
if("gplots" %in% rownames(installed.packages()) == FALSE){
  install.packages("gplots")}
require(gplots)
heatmap_matrix= as.matrix(TMT[,21:26])
heatmap_matrix= na.omit(heatmap_matrix)
colnames(heatmap_matrix)= c(paste("CTRL", 1:3), paste("shRNA", 1:3))
tiff("~/heatmap_KG1_chir.tiff",width=1200,height= 800,units = "px")
heatmap.2(heatmap_matrix,main="KG1 Sample heatmap with dendrogram",cexRow = 1,density.info="none",
          trace="none",key=T,srtCol=30,labRow = NA, col = "bluered",
          labCol = colnames(matrix),dendrogram = 'both',scale = "row",
          ylab="",xlab="Samples")
dev.off()