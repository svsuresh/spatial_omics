library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(rhdf5)
library(Matrix)
library(sctransform)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Hmisc)
library(ReactomePA)
library(mygene)
library(stringr)
library(clusterProfiler)

####################################################################################
##############################Total transcripts and gene coutnts. R
##script##########
####################################################################################
##read in the coordinates of points lying on top of the tissue.position.txt is
##generated from matlab script "Pixel_identification.m".
location <-
  read.table(
    "external_input/position.txt",
    sep = ",",
    header = FALSE,
    dec = ".",
    stringsAsFactors = F
  )
x <- as.character(location[1,])
x = x[-1]
##read expression matrix and generate the Filtered_matrix.tsv, which contains only the useful pixels
my_data <-
  read.table(
    file = "results/npg_test_stdata.tsv",
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
  )
data_filtered <- my_data[my_data$X %in% x,]
write.table(
  data_filtered,
  file = 'r_results/Filtered_matrix.tsv',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

##calculate the total UMI count and Gene count
count <- rowSums(data_filtered[, 2:ncol(data_filtered)])
data_filtered_binary <-
  data_filtered[, 2:ncol(data_filtered)] %>%
  mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

##UMI Count
region <-
  2500  #change the x axis maximum, need to adjust based on different sample
df <- data.frame(number = 1, c = count)
pdf(
  file = paste("r_results/UMI.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(df, aes(x = c), color = 'blue', xlab = "Gene") +
  geom_histogram(
    aes(y = ..density..),
    binwidth = region / 20,
    color = "black",
    fill = "white",
    size = 1
  ) +
  geom_density(
    alpha = .2,
    fill = "#FF6666",
    size = 1,
    color = "red"
  ) +
  scale_x_continuous(name = "UMI", limits = c(0, region)) +
  scale_y_continuous(name = "Density", expand = c(0, 0)) +
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(colour = "black", size = 20),
    axis.title = element_text(
      colour = "black",
      size = 25,
      face = "bold"
    ),
    legend.text = element_text(colour = "black", size = 20),
    legend.title = element_text(
      colour = "black",
      size = 20,
      face = "bold"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    )
  )
dev.off()

##Gene Count
df <- data.frame(number = 1, c = gene_count)
region = 1500 #change the x axis maximum, need to adjust based on different sample
pdf(
  file = paste("r_results/Gene.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(df, aes(x = c), color = 'blue', xlab = "Gene") +
  geom_histogram(
    aes(y = ..density..),
    binwidth = region / 20,
    color = "black",
    fill = "white",
    size = 1
  ) +
  geom_density(
    alpha = .2,
    fill = "#FF6666",
    size = 1,
    color = "red"
  ) +
  scale_x_continuous(name = "Gene", limits = c(0, region)) +
  scale_y_continuous(name = "Density", expand = c(0, 0)) +
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(colour = "black", size = 20),
    axis.title = element_text(
      colour = "black",
      size = 25,
      face = "bold"
    ),
    legend.text = element_text(colour = "black", size = 20),
    legend.title = element_text(
      colour = "black",
      size = 20,
      face = "bold"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'black',
      size = 0.5,
      linetype = 'solid'
    )
  )
dev.off()

#imported_raster=OpenImageR::readImage("ventricle.jpg")     #if you want the
#microscope image under the heatmap, then uncomment this line.
imported_raster = OpenImageR::readImage("external_input/FFPE-2.jpg")
g <-
  rasterGrob(
    imported_raster,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = FALSE
  )

#UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to
#be close to the maximum number.
test <- data_filtered %>%
  separate(X, c("A", "B"),  sep = "x")

pdf(
  file = paste("r_results/UMI_heatmap.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) +
  #scale_color_gradientn(colours = c("black", "green")) +
  scale_color_gradientn(
    colours = c("blue", "green", "red"),
    limits = c(0, 1000),
    oob = scales::squish
  ) +
  ggtitle("UMI") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you
  #want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name = "X",
                     limits = c(NA, NA),
                     expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name = "Y",
                  limits = c(NA, NA),
                  expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim = c(0, 51), ylim = c(51, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    #legend.title = element_text(colour="black", size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()

#Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to
#be close to the maximum number.
pdf(
  file = paste("r_results/Gene_heatmap.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = gene_count)) +
  scale_color_gradientn(
    colours = c("blue", "green", "red"),
    limits = c(0, 1000),
    oob = scales::squish
  ) +
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you
  #want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name = "X",
                     limits = c(NA, NA),
                     expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name = "Y",
                  limits = c(NA, NA),
                  expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim = c(0, 51), ylim = c(51, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    #legend.title = element_text(colour="black", size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()
############################################################################################
######################DEG.R script ##############################
#change filename1 to name of txt file you want to load

data1 <- data_filtered
row.names(data1) = data1[, 1]
data1 = data1[, -1]
data2 <- t(data1)
sample1.name <- "npgdata"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)

#Create Seurat object
ffpe2 <-
  CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
ffpe2 <-
  PercentageFeatureSet(ffpe2, pattern = "^MT-", col.name = "percent.mt")
ffpe2 <-
  SCTransform(ffpe2, vars.to.regress = "percent.mt", verbose = FALSE)
ffpe2 <- RunPCA(ffpe2, verbose = FALSE)
ffpe2 <- RunUMAP(ffpe2, dims = 1:10, verbose = FALSE)
ffpe2 <- FindNeighbors(ffpe2, dims = 1:10, verbose = FALSE)
ffpe2 <- FindClusters(ffpe2, resolution = 0.8, verbose = FALSE)

pdf(
  file = paste("r_results/umap_plot2.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
DimPlot(ffpe2, label = TRUE) + NoLegend()  # the UMAP plot
dev.off()

ffpe2.markers <-
  FindAllMarkers(
    ffpe2,
    only.pos = TRUE,
    min.pct = 0,
    logfc.threshold = 0.01
  )
ffpe2 <- ScaleData(object = ffpe2, features = rownames(ffpe2))
##View(ffpe2.markers)
top10 <-
  ffpe2.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
pdf(
  file = paste("r_results/top10_heatmap.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
DoHeatmap(ffpe2, features = top10$gene) +
  scale_fill_gradientn(colors = c("red", "black", "green"))
dev.off()
write.table(top10, "r_results/top10.txt", sep = "\t")
go <-
  ffpe2.markers %>%
  group_by(cluster) %>%
  top_n(n = 1000, wt = avg_log2FC)
write.table(go, "r_results/go3.txt", sep = "\t")
x = ffpe2.markers$cluster
aa = ffpe2.markers$p_val_adj >= 0.01
ffpe2.markers$pp = as.numeric(aa)
ffpe2.markers$pp <- as.factor(ffpe2.markers$pp)

pdf(
  file = paste("r_results/significant_genes.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(ffpe2.markers, aes(x = cluster, y = avg_log2FC, color = pp)) +
  geom_violin(color = NA, fill = NA) +
  geom_jitter(shape = 16,
              size = 1,
              position = position_jitter(0.2)) +
  ylab("average log2FC") + xlab("Cluster") +
  geom_hline(yintercept = 0, size = 1) +
  scale_color_manual(values = c("red", "black")) +
  geom_text(aes(label = ifelse(
    !(avg_log2FC < 0.25 &
        avg_log2FC > -0.25) &
      pp == 0,
    as.character(gene),
    ''
  )),
  position = position_jitter(width = 0.5, height = 0.3),
  colour = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(colour = "black", size = 15),
    axis.title = element_text(
      colour = "black",
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_blank(),
    #legend.title = element_text(colour="black", size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank()
  )
dev.off()

## not working
####VlnPlot(ffpe2.markers)

deg <- go

top5 <-
  ffpe2.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = p_val)

pdf(
  file = paste("r_results/violin_plot_igf2_prrx1.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)

VlnPlot(ffpe2, features = c("Igf2", "Prrx1"))

dev.off()

pdf(
  file = paste("r_results/dot_plot_top10_genes.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
DotPlot(
  ffpe2,
  features = unique(top10$gene),
  cols = c("blue", "red"),
  dot.scale = 8,
) +
  RotatedAxis()
dev.off()
#################################code end##################
##########################################################
#####################individual_gene_plot code goes here########
######################################
#Create Seurat object
ffpe3 <-
  CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
ffpe3 <-
  NormalizeData(ffpe3,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

genedata <- ffpe3[["RNA"]]@data
genedata <- t(genedata)
gene <- as.data.frame(as.matrix(genedata))
gene$X = row.names(gene)
##View(gene$X)
test <- gene %>%
  separate(X, c("A", "B"),  sep = "x")

#UMI heatmap
pdf(
  file = paste("r_results/Epcam.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = Epcam)) +
  #scale_color_gradientn(colours = c("black", "green")) +
  scale_color_gradientn(colours = c("blue", "green", "red"),
                        oob = scales::squish) +
  ggtitle("Epcam") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name = "X",
                     limits = c(NA, NA),
                     expand = expansion(mult = c(-0.013,-0.013))) +
  scale_y_reverse(name = "Y",
                  limits = c(NA, NA),
                  expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim = c(0, 51), ylim = c(51, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    axis.text = element_text(colour = "black", size = 30),
    axis.title = element_text(
      colour = "black",
      size = 30,
      face = "bold"
    ),
    legend.text = element_text(colour = "black", size = 30),
    legend.title = element_blank(),
    #legend.title = element_text(colour="black", size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(size = 1, fill = NA)
  )
dev.off()

###################code ends here####################################
##########################################################
##############clustering code start#######################

#Create Seurat object
pdf(
  file = paste("r_results/dimplot_seurat_object.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)

DimPlot(ffpe2, label = TRUE) +
  NoLegend()  # the UMAP plot

dev.off()

#Identify the X and Y coordinates
#head(Idents(ffpe2), 5)
#Idents(ffpe2)
ident <- Idents(ffpe2)
df <- data.frame(ident[])
df1 <- data.frame(X = row.names(df), count = df$ident..)
test <- df1 %>%
  separate(X, c("A", "B"),  sep = "x")

#Plot the spatial clusters
pdf(
  file = paste("r_results/spatial_clusters.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) +
  #scale_color_gradientn(colours = c("black", "green")) +
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  annotation_custom(
    g,
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  ) +
  geom_point(shape = 15, size = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name = "X",
                     limits = c(NA, NA),
                     expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name = "Y",
                  limits = c(NA, NA),
                  expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim = c(0, 51), ylim = c(51, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    #legend.title = element_text(colour="black", size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()
####################################clustering code ends here################
####################################################################
#########enrichgo code goes here##################
####################################################

mylist.names = c(paste0("X", seq (0, 8)))
mylist <- vector("list", length(mylist.names))
names(mylist) <- mylist.names
for (i in 0:8) {
  x <- deg$gene[deg$cluster == i]
  gene.df <- clusterProfiler::bitr(
    toupper(x),
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  mylist[[i + 1]] = gene.df$ENTREZID
}
print ("list prepared")

ck <-
  clusterProfiler::compareCluster(geneCluster = mylist,
                                  fun = "enrichGO",
                                  OrgDb = 'org.Hs.eg.db')

print ("compare cluster done")
save.image("compared_cluster.Rdata")
print ("saved image")
pdf(
  file = paste("r_results/dotplot_cp.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)

clusterProfiler::dotplot(ck)

dev.off()
#View(ck)

# ck@compareClusterResult$Description = NULL
fcc <- gofilter(ck, level = 3)
#View(cc)
#cc@compareClusterResult$Description
scc <- simplify(ck)

pdf(
  file = paste("r_results/enrichgo_cp.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)
clusterProfiler::dotplot(fcc)
clusterProfiler::dotplot(scc)

dev.off()

# ck@compareClusterResult$Description
### Do not know the purpose of this code below
# #View(gene.df)
#mapIds(org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL')
#keys()
#################################code ends here###################
#############rungo script starts here####################
###############################################
genetable <-  data.frame(cluster = go$cluster, gene = go$gene)
temp = genetable[genetable$cluster == 0, ]

for (i in c(0:8)) {
  temp = genetable[genetable$cluster == i,]
  xli = temp$gene
result <-
    queryMany(
      xli,
      scopes = "symbol",
      fields = c("uniprot", "ensembl.gene", "reporter"),
      species = "mouse"
    )
  mylist[[i + 1]] = result$`_id`
}

mylist <- lapply(mylist, function(x)
  x[!is.na(x)])

mylist <- lapply(mylist, function(x)
  str_remove(x, "ENSMUSG000000"))

#mylist

res <- compareCluster(mylist, fun = "enrichPathway")

pdf(
  file = paste("r_results/enrich_pw_run_go.pdf", sep = ""),
  width = 8.6,
  height = 8.6
)

dotplot(res)

dev.off()

######################code ends here###############################
######################################################################
### Figures code################
