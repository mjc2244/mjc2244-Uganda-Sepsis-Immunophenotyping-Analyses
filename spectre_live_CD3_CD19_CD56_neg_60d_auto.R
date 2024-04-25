#Clear R environment
rm(list = ls())

library("Spectre")

#Set PrimaryDirectory
getwd()
PrimaryDirectory <- getwd()

#Set 'input' directory
setwd("C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

#Set 'metadata' directory
setwd(PrimaryDirectory)
setwd("C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg\\Live_myeloid_metadata_CD3_CD19_CD56_neg")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

#Create output directory
dir.create("Output_Spectre", showWarnings = FALSE)
setwd("C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

#Import data
setwd(InputDirectory)
list.files(InputDirectory, ".fcs")

data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".fcs",
                                 do.embed.file.names = TRUE)

warnings()

#Check the data
check <- do.list.summary(data.list)
check$name.table # Review column names and their subsequent values
check$ncol.check # Review number of columns (features, markers) in each sample
check$nrow.check # Review number of rows (cells) in each sample
data.list[[2]]

#Merge the data
cell.dat <- Spectre::do.merge.files(dat = data.list)
cell.dat <- cell.dat[,c(7:9,12:17,19:23)]
cell.dat

#Read in metadata
setwd(MetaDirectory)
meta.dat <- fread("pbmc_live_lymph_cd3_cd19_cd56_neg_metadata.csv")
meta.dat

#Data transformation
setwd(OutputDirectory)
dir.create("Output 1 - transformed plots")
setwd("Output 1 - transformed plots")
#Arcsinh transformation
as.matrix(names(cell.dat))

to.asinh <- names(cell.dat)[c(1:13)]
to.asinh

cofactor <- 150

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

#Add metadata
meta.dat

sample.info <- meta.dat[,c(1:4)]
sample.info

cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)

cell.dat

cell.dat <- na.omit(cell.dat)

#Rename asinh-transformed columns
colnames(cell.dat)[colnames(cell.dat) == "FJComp-APC-A_BDCA2_asinh"] = "BDCA2"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-APC-Cy7-A_HLADR_asinh"] ="HLA-DR"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-Alexa Fluor 700-A_CD16_asinh"] ="CD16"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BUV395-A_CD56_asinh"] ="CD56"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BUV737-A_CD15_asinh"] ="CD15"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BV421-A_CD80_asinh"] ="CD80"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BV510-A_CD33_asinh"] ="CD33"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BV605-A_CD123_asinh"] ="CD123"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BV711-A_CD11c_asinh"] ="CD11c"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-BV785-A_CD14_asinh"] ="CD14"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-PE-A_CD155_asinh"] ="CD155"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-PE-Cy7-A_PD-L1_asinh"] ="PD-L1"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-PE-Texas Red-A_CD86_asinh"] ="CD86"
colnames(cell.dat)[colnames(cell.dat) == "FJComp-PerCP-Cy5-5-A_CD40_asinh"] ="CD40"

#Define cellular columns
as.matrix(names(cell.dat))

cellular.cols <- names(cell.dat)[c(15:27)]
as.matrix(cellular.cols)

cluster.cols <- names(cell.dat)[c(15:27)]
as.matrix(cluster.cols)

sample_id.col <- "Sample"
Group_2.col <- "Vital Status"

#Subsample targets per group
data.frame(table(cell.dat[["Group_2"]]))

sub.targets <- c(4000, 4000) # target subsample numbers from each group
sub.targets

#FlowSOM clustering
setwd(OutputDirectory)
dir.create("Output - clustering")
setwd("Output - clustering")
#Clustering
cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = "auto", clust.seed = 56789, meta.seed = 56789)
fwrite(cell.dat, "clustered.data.csv")

#Dimensionality reduction
cell.sub <- do.subsample(cell.dat, sub.targets, "Group_2")
cell.sub

cell.sub <- run.umap(cell.sub, cluster.cols)
fwrite(cell.sub, "clustered.data.DR.csv")

#DR plots
setwd("C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg")
#By FlowSOM metacluster
plot1 <- make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor',
                          add.label = FALSE, dot.size = 0.05, blank.axis = TRUE)
plot1 <- plot1 + scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#91D1C2FF", "#7E614EFF", "#8491B4FF", "#EFC000FF", "#3B3B3BFF")) 
plot1 <- plot1 + theme(plot.title = element_text(color="black", size=10), legend.text = element_text(color = "black", size = 14)) 
plot1 <- plot1 + guides(colour = guide_legend(override.aes = list(size=1.5)))
plot1 <- plot1 + ggtitle("") 
plot1 <- plot1 + theme(legend.position = "right") + theme(legend.key=element_blank())
plot1

#By marker expression
plot2 <- make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols, colours = "jet", blank.axis = TRUE, dot.size = 0.05)
plot2 <- plot2 + theme(plot.title = element_text(color="black", size=10), legend.text = element_text(color = "black", size = 14))
plot2 

#By 60d vital status - density                      
cell.sub$Group <- as.factor(cell.sub$Group)
cell.sub$Group_2 <- as.factor(cell.sub$Group_2)

plot3 <- make.colour.plot(dat = subset(cell.sub, Group_2 == "Alive"), x.axis = "UMAP_X",y.axis = "UMAP_Y", dot.size = 0.05, blank.axis = TRUE, colours = "inferno")
plot3 <- plot3 + theme(plot.title = element_text(color="black", size=10), legend.text = element_text(color = "black", size = 14))
plot3 <- plot3 + ggtitle("") 
plot3 <- plot3 + theme(legend.position = "right")
plot3

plot4 <- make.colour.plot(dat = subset(cell.sub, Group_2 == "Dead"), x.axis = "UMAP_X",y.axis = "UMAP_Y", dot.size = 0.05, blank.axis = TRUE, colours = "inferno")
plot4 <- plot4 + theme(plot.title = element_text(color="black", size=10), legend.text = element_text(color = "black", size = 14))
plot4 <- plot4 + ggtitle("")
plot4 <- plot4 + theme(legend.position = "right")
plot4

#By 60d vital status - factor
plot5 <- make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", col.axis = "Group_2", col.type = 'factor', dot.size = 0.05, add.label = FALSE, blank.axis = TRUE)
plot5 <- plot5 + scale_colour_manual(values = c("#0099b4", "#925e9f"))
plot5 <- plot5 + theme(plot.title = element_text(color="black", size=10), legend.text = element_text(color = "black", size = 14))
plot5 <- plot5 + guides(colour = guide_legend(override.aes = list(size=1.5)))
plot5 <- plot5 + ggtitle("") 
plot5 <- plot5 + theme(legend.position = "right")
plot5

blankplot <- ggplot() + theme_void()

#Combine plots
library(ggpubr)
flowsomplot1 <- ggarrange(plot3, plot4, plot1, blankplot, hjust = c(-0.015,-0.03,-0.01,-0.01), ncol = 2, nrow = 2)
flowsomplot1

#Expression heatmap by cluster
exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cutree_rows = 1, cellular.cols, standard.colours = "Blues", plot.title = "", cell.size = 10, dendrograms = "none")

#Summary data 
setwd(OutputDirectory)
dir.create("Output 4 - summary data")
setwd("Output 4 - summary data")

#Select columns to measure proportions
as.matrix(cellular.cols)

#Create summary tables
sum.dat <- create.sumtable(dat = cell.dat,
                           sample.col = "FileName",
                           pop.col = "FlowSOM_metacluster",
                           annot.cols = c("Group_2"))

#Review summary data
sum.dat

#NA to zero
sum.dat[is.na(sum.dat)] <- 0

#Export sum.dat for further analysis
write.csv(sum.dat,"C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg\\sumdat.csv")

#Convert sum.dat to matrix
as.matrix(names(sum.dat))

#Make table of cluster proportions by vital status
library(tableone)
#Create a variable list which we want in Table
listVars <- c("Percent of sample -- 1",                      
              "Percent of sample -- 2",                     
              "Percent of sample -- 3",                    
              "Percent of sample -- 4",                   
              "Percent of sample -- 5",                  
              "Percent of sample -- 6",
              "Percent of sample -- 7",                    
              "Percent of sample -- 8",                   
              "Percent of sample -- 9")

catVars <- c()

tablecd3_cd19_cd56_neg_prop <- CreateTableOne(vars = listVars, 
                           data = sum.dat, 
                           factorVars = catVars,
                           strata = "Group_2",
                           includeNA = TRUE,
                           addOverall = TRUE)
tablecd3_cd19_cd56_neg_prop
tablecd3_cd19_cd56_neg_prop

tablecd3_cd19_cd56_neg_prop <- print(tablecd3_cd19_cd56_neg_prop, nonnormal = c("Percent of sample -- 1",                      
                                                        "Percent of sample -- 2",                     
                                                        "Percent of sample -- 3",                    
                                                        "Percent of sample -- 4",                   
                                                        "Percent of sample -- 5",                  
                                                        "Percent of sample -- 6",
                                                        "Percent of sample -- 7",                    
                                                        "Percent of sample -- 8",                   
                                                        "Percent of sample -- 9"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tablecd3_cd19_cd56_neg_prop), file = 'table_cd3_cd19_cd56_neg_prop.html')

#Rename columns
colnames(sum.dat)[3] ="Percent_of_sample_1"
colnames(sum.dat)[5] ="Percent_of_sample_3"
colnames(sum.dat)[6] ="Percent_of_sample_4"
colnames(sum.dat)[10] ="Percent_of_sample_8"

#Stripchart plots of cluster abundance by 60d vital status 
library(ggpubr)
library(ggprism)
#Myeloid/NK Cluster 1
myeloid_dc_nk_cluster1_plot <- ggstripchart(sum.dat, x = "Group_2", y = "Percent_of_sample_1", 
                                      title = "",
                                      color = "Group_2", fill = "Group_2", alpha = 1, size = 4, palette =  c("#3b4992","#ee0000"),
                                      order = c("Alive", "Dead"), repel = TRUE) +
  theme(legend.position = "none", text = element_text(size = 12)) 
myeloid_dc_nk_cluster1_plot  <- myeloid_dc_nk_cluster1_plot   + stat_summary(fun.y = median, 
                                                                               fun.ymin = median, 
                                                                               fun.ymax = median, 
                                                                               geom = "crossbar", 
                                                                               width = 0.5)
myeloid_dc_nk_cluster1_plot  <- myeloid_dc_nk_cluster1_plot  + theme(legend.position = "none") + stat_compare_means(label.y = 45, label.x = 1.5, size=4.5)
myeloid_dc_nk_cluster1_plot  <- myeloid_dc_nk_cluster1_plot  + labs(x ="", y = expression("% live CD3" ^ "neg"*"CD19" ^ "neg"*"CD56" ^ "neg"* " cells"), title = expression("CD56" ^ "neg"*"CD16" ^ "bright"*" NK cells"*" (1)"))
myeloid_dc_nk_cluster1_plot <- myeloid_dc_nk_cluster1_plot + theme(plot.title = element_text(size = 12))
myeloid_dc_nk_cluster1_plot

#Myeloid/NK Cluster 3
myeloid_dc_nk_cluster3_plot <- ggstripchart(sum.dat, x = "Group_2", y = "Percent_of_sample_3", 
                                      title = "",
                                      color = "Group_2", fill = "Group_2", alpha = 1, size = 4, palette =  c("#3b4992","#ee0000"),
                                      order = c("Alive", "Dead"), repel = TRUE) +
  theme(legend.position = "none", text = element_text(size = 12)) 
myeloid_dc_nk_cluster3_plot  <- myeloid_dc_nk_cluster3_plot   + stat_summary(fun.y = median, 
                                                                 fun.ymin = median, 
                                                                 fun.ymax = median, 
                                                                 geom = "crossbar", 
                                                                 width = 0.5)
myeloid_dc_nk_cluster3_plot  <- myeloid_dc_nk_cluster3_plot  + theme(legend.position = "none") + stat_compare_means(label.y = 11, label.x = 1.5, size=4.5)
myeloid_dc_nk_cluster3_plot  <- myeloid_dc_nk_cluster3_plot  + labs(x ="", y = expression("% live CD3" ^ "neg"*"CD19" ^ "neg"*"CD56" ^ "neg"* " cells"), title = expression("CD16" ^ "bright"*"PD-L1" ^ "bright"*" PMN-MDSCs"*" (3)"))
myeloid_dc_nk_cluster3_plot <- myeloid_dc_nk_cluster3_plot + theme(plot.title = element_text(size = 12))
myeloid_dc_nk_cluster3_plot

#Myeloid/NK Cluster 4
myeloid_dc_nk_cluster4_plot <- ggstripchart(sum.dat, x = "Group_2", y = "Percent_of_sample_4", 
                                      title = "",
                                      color = "Group_2", fill = "Group_2", alpha = 1, size = 4, palette =  c("#3b4992","#ee0000"),
                                      order = c("Alive", "Dead"), repel = TRUE) +
  theme(legend.position = "none", text = element_text(size = 12)) 
myeloid_dc_nk_cluster4_plot  <- myeloid_dc_nk_cluster4_plot   + stat_summary(fun.y = median, 
                                                                 fun.ymin = median, 
                                                                 fun.ymax = median, 
                                                                 geom = "crossbar", 
                                                                 width = 0.5)
myeloid_dc_nk_cluster4_plot  <- myeloid_dc_nk_cluster4_plot  + theme(legend.position = "none") + stat_compare_means(label.y = 52, label.x = 1.5, size=4.5)
myeloid_dc_nk_cluster4_plot  <- myeloid_dc_nk_cluster4_plot  + labs(x ="", y = expression("% live CD3" ^ "neg"*"CD19" ^ "neg"*"CD56" ^ "neg"* " cells"), title = expression("CD16" ^ "neg"*"PD-L1" ^ "dim"*" PMN-MDSCs"*" (4)"))
myeloid_dc_nk_cluster4_plot <- myeloid_dc_nk_cluster4_plot + theme(plot.title = element_text(size = 12))
myeloid_dc_nk_cluster4_plot

#Myeloid/NK Cluster 8
myeloid_dc_nk_cluster8_plot <- ggstripchart(sum.dat, x = "Group_2", y = "Percent_of_sample_8", 
                                      title = "",
                                      color = "Group_2", fill = "Group_2", alpha = 1, size = 4, palette =  c("#3b4992","#ee0000"),
                                      order = c("Alive", "Dead"), repel = TRUE) +
  theme(legend.position = "none", text = element_text(size = 12)) 
myeloid_dc_nk_cluster8_plot  <- myeloid_dc_nk_cluster8_plot   + stat_summary(fun.y = median, 
                                                                 fun.ymin = median, 
                                                                 fun.ymax = median, 
                                                                 geom = "crossbar", 
                                                                 width = 0.5)
myeloid_dc_nk_cluster8_plot  <- myeloid_dc_nk_cluster8_plot  + theme(legend.position = "none") + stat_compare_means(label.y = 9, label.x = 1.5, size=4.5)
myeloid_dc_nk_cluster8_plot  <- myeloid_dc_nk_cluster8_plot  + labs(x ="", y = expression("% live CD3" ^ "neg"*"CD19" ^ "neg"*"CD56" ^ "neg"* " cells"), title = expression("Plasmacytoid DCs"*" (8)"))
myeloid_dc_nk_cluster8_plot <- myeloid_dc_nk_cluster8_plot + theme(plot.title = element_text(size = 12))
myeloid_dc_nk_cluster8_plot

#Combine plots 
library(ggpubr)
myeloid_dc_nk_clusterplots <- ggarrange(myeloid_dc_nk_cluster3_plot, 
          myeloid_dc_nk_cluster4_plot,
          myeloid_dc_nk_cluster8_plot, myeloid_dc_nk_cluster1_plot,
          ncol = 4, nrow = 1)
myeloid_dc_nk_clusterplots

#Apply MEM
library(MEM)
library(cytoMEM)
library(gplots)
library(ggplot2)
library(hexbin)
library(viridis)
library(ggExtra)
library(FlowSOM)
library(flowCore)
library(Biobase)
library(tidyverse)
library(uwot)
library(RColorBrewer)
library(dplyr)

setwd("C:\\Users\\mjcum\\OneDrive\\Desktop\\Live_myeloid_pregate_CD3_CD19_CD56_neg\\MEM")
getwd()

cell.dat.mem <- cell.dat[,c(15:27,32)]
View(cell.dat.mem)

cell.dat.mem <- cell.dat.mem %>% rename("cluster" = "FlowSOM_metacluster" )
cell.dat.mem <- as.matrix(cell.dat.mem)
MEMvalues <- MEM(cell.dat.mem,
                 transform=FALSE,
                 choose.markers=FALSE,
                 markers="all",
                 rename.markers=FALSE,
                 choose.ref=FALSE,
                 zero.ref=FALSE,
                 file.is.clust = FALSE,
                 add.fileID = FALSE,
                 IQR.thresh = NULL)

str(MEMvalues)
MEMvalues$MAGpop
MEMvalues$MAGref
MEMvalues$IQRpop
MEMvalues$IQRref
MEMvalues$MEM_matrix
build_heatmaps(
  MEMvalues,
  cluster.MEM = "none",
  cluster.medians = "none",
  cluster.IQRs = "none",
  display.thresh = 2,
  output.files = TRUE,
  labels = FALSE,
  only.MEMheatmap = TRUE)



