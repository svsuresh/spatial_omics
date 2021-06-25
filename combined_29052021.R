setRepositories(ind = 1:5)

pkgs = c( "ggplot2", "plyr", "gridExtra", "magrittr", "tidyr", "raster", "OpenImageR", "ggpubr","grid",
          "wesanderson", "dplyr", "Seurat", "SeuratObject", "patchwork", "rhdf5", "Matrix", "sctransform",
          "org.Hs.eg.db",  "clusterProfiler",  "Hmisc", "ReactomePA", "mygene", "stringr","clusterProfiler"
)

#lapply(pkgs, library, character.only=T)

lapply(pkgs, function(x) {
  if (!require(x, character.only = T)) {
    install.packages(x, dependencies = T)
    library(x)
  }
})

args = commandArgs(trailingOnly = TRUE)
out_directory=paste(file.path(args[4]),Sys.Date(), sep="_")

if (!dir.exists(out_directory)){
  print  (paste(out_directory,"doesn't exist. Creating", out_directory, sep = " "))
  dir.create(out_directory)
  print(paste ("all results will be saved to ",file.path(getwd(),out_directory), sep = ""))
}


# Total transcripts and gene counts. R script read in  the coordinates of points lying on top of the tissue.position.txt is generated from matlab script "Pixel_identification.m".

location <-
  read.table(
    args[1],
    sep = ",",
    header = FALSE,
    dec = ".",
    stringsAsFactors = F
  )
x <- as.character(location[1, ])
x = x[-1]

### read image
#imported_raster=OpenImageR::readImage("ventricle.jpg")     #if you want the
#microscope image under the heatmap, then uncomment this line.
imported_raster = OpenImageR::readImage(args[2])

g <-
  rasterGrob(
    imported_raster,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = FALSE
  )

##read expression matrix and generate the Filtered_matrix.tsv, which contains only the useful pixels
my_data <-
  read.table(
    file = args[3],
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
  )
data_filtered <- my_data[my_data$X %in% x, ]

write.table(
  data_filtered,
  file = file.path(out_directory, "Filtered_matrix.tsv"),
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)
# View(data_filtered)

##calculate the total UMI count and Gene count

count = rowSums(data_filtered[, -1])

theme_base= theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25,face = "bold"),
      axis.text = element_text(colour = "black", size = 20),
      axis.title = element_text(colour = "black",size = 25,face = "bold"),
      axis.line.x = element_line(colour = 'black',size = 0.5,linetype = 'solid'),
      axis.line.y = element_line(colour = 'black',size = 0.5,linetype = 'solid')
    )

ggplot_draw = function(input, region, type, color) {
  ggplot(as.data.frame(input), aes(x = input)) +
    geom_histogram(
      aes(y = ..density..),
      binwidth = region / 20,
      colour = "white",
      fill = color
    ) +
    geom_density(
      alpha = .2,
      fill = "red",
      color = "red",
      alpha = 0.8
    ) +
    labs(x = type, y = "Density", title = paste(type, "density", sep = "\t")) +
    xlim(c(0, region)) +
    theme_base
      }

##UMI Count

region <-
  2500  #change the x axis maximum, need to adjust based on different sample

file.path(out_directory, "UMI.pdf")

pdf(
  file = file.path(out_directory, "UMI.pdf"),
  width = 8.6,
  height = 8.6
)
fig1=ggplot_draw(count, region, "UMI", "steelblue")
fig1
dev.off()

##Gene Count

genecounts = sapply(data_filtered[, -1], as.logical) %>%
  rowSums()

#change the x axis maximum, need to adjust based on different sample
region = 1500

pdf(
  file = file.path(out_directory, "Gene.pdf"),
  width = 8.6,
  height = 8.6
)

fig2=ggplot_draw(genecounts, region, "gene", "darkblue")
fig2
dev.off()

# grid.arrange(fig1, fig2)

#UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to
#be close to the maximum number.

ggplot_heatmap = function(coords, data, title) {
  ggplot(coords, aes(
    x = as.numeric(A),
    y = as.numeric(B),
    color = !!data)
  ) +
    ggtitle(title) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    geom_point(shape = 15, size = 3) +
    expand_limits(x = 0, y = 0)  +
    scale_x_continuous(
      name = "X",
      limits = c(NA, NA),
      expand = expansion(mult = c(-0.013, -0.013))
    ) +
    scale_y_reverse(
      name = "Y",
      limits = c(NA, NA),
      expand = expansion(mult = c(-0.013, 0.008))
    ) +
    coord_equal(xlim = c(0, 51), ylim = c(51, 1)) +
    theme_base+
    theme(legend.title=element_blank())
}

test = data_filtered %>%
  select(X) %>%
  separate(X, c("A", "B"),  sep = "x")

pdf(
  file = file.path(out_directory, "UMI_heatmap.pdf"),
  width = 8.6,
  height = 8.6
)

fig3=ggplot_heatmap(test, count,"UMI")+
  scale_color_gradientn(
    colours = c("blue", "green", "red"),
    limits = c(0, 1000),
    oob = scales::squish
  )
fig3

dev.off()

#Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to
#be close to the maximum number.
pdf(
  file = paste("Gene_heatmap.pdf"),
  width = 8.6,
  height = 8.6
)

fig4=ggplot_heatmap(test, genecounts,"Gene") +
  scale_color_gradientn(
    colours = c("blue", "green", "red"),
    limits = c(0, 1000),
    oob = scales::squish
  )

fig4
dev.off()

# DEG.R script
# change filename1 to name of txt file you want to load

data1 <- data_filtered
row.names(data1) = data1[, 1]
data1 = data1[, -1]
#View(data1)
data2 <- t(data1)
sample1.name <- "npgdata"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
#View(matrix1.data)

#Create Seurat object
ffpe1 = CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
ffpe2 = PercentageFeatureSet(ffpe1, pattern = "^MT-", col.name = "percent.mt")
ffpe2 = SCTransform(ffpe2, vars.to.regress = "percent.mt", verbose = FALSE)
ffpe2 =  RunPCA(ffpe2, verbose = FALSE)
ffpe2 =  RunUMAP(ffpe2, dims = 1:10, verbose = FALSE)
ffpe2 =  FindNeighbors(ffpe2, dims = 1:10, verbose = FALSE)
ffpe2 =  FindClusters(ffpe2, resolution = 0.8, verbose = FALSE)

pdf(
  file = file.path(out_directory, "umap_plot2.pdf"),
  width = 8.6,
  height = 8.6
)
fig5=DimPlot(ffpe2, label = TRUE) + NoLegend()  # the UMAP plot
fig5
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
  file = file.path(out_directory, "top10_heatmap.pdf"),
  width = 8.6,
  height = 8.6
)

fig6=DoHeatmap(ffpe2, features = top10$gene) +
  scale_fill_gradientn(colors = c("red", "black", "green"))
fig6
dev.off()

write.table(top10, file.path(out_directory,"top10.txt"), sep = "\t")

go <-
  ffpe2.markers %>%
  group_by(cluster) %>%
  top_n(n = 1000, wt = avg_log2FC)

write.table(go, file.path(out_directory,"top_thousand_per_cluster.txt"), sep = "\t")

x = ffpe2.markers$cluster
aa = ffpe2.markers$p_val_adj >= 0.01
ffpe2.markers$pp = as.numeric(aa)
ffpe2.markers$pp <- as.factor(ffpe2.markers$pp)

pdf(
  file = file.path(out_directory, "significant_genes.pdf"),
  width = 8.6,
  height = 8.6
)

fig7=ggplot(ffpe2.markers, aes(x = cluster, y = avg_log2FC, color = pp)) +
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
  theme_base
fig7

dev.off()

## not working
####VlnPlot(ffpe2.markers)

deg <- go

top5 <-
  ffpe2.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = p_val)

pdf(
  file = file.path(out_directory, "violin_plot_igf2_prrx1.pdf"),
  width = 8.6,
  height = 8.6
)

fig8=VlnPlot(ffpe2, features = c("Igf2", "Prrx1"))

dev.off()

pdf(
  file = file.path(out_directory, "dot_plot_top10_genes.pdf"),
  width = 8.6,
  height = 8.6
)
fig9 = DotPlot(
  ffpe2,
  features = unique(top10$gene),
  cols = c("blue", "red"),
  dot.scale = 8,
) +
  RotatedAxis()
fig9
dev.off()
#################################code end##################
##########################################################
#####################individual_gene_plot code goes here########
######################################
#Create Seurat object
ffpe3 <-
  NormalizeData(ffpe1,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

genedata <- ffpe3[["RNA"]]@data
genedata <- t(genedata)
gene <- as.data.frame(as.matrix(genedata))
gene$X = row.names(gene)
##View(gene$X)
gene = gene %>%
  separate(X, c("A", "B"),  sep = "x")
getwd()
#save.image("st_pipeline_30052021.Rdata")
#UMI heatmap
pdf(
  file = file.path(out_directory, "Epcam.pdf"),
  width = 8.6,
  height = 8.6
)

fig10 = ggplot_heatmap(gene,sym("Epcam"),"Epcam")+
  scale_color_gradientn(
    colours = c("blue", "green", "red"),
    oob = scales::squish
  )

dev.off()

###################code ends here####################################
##########################################################
##############clustering code start#######################

pdf(
  file = file.path(out_directory, "dimplot_seurat_object.pdf"),
  width = 8.6,
  height = 8.6
)

fig11=DimPlot(ffpe2, label = TRUE) +
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
  file = file.path(out_directory, "spatial_clusters.pdf"),
  width = 8.6,
  height = 8.6
)

fig12=ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) +
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
  theme_base
fig12

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

ck = clusterProfiler::compareCluster(geneCluster = mylist,
                                  fun = "enrichGO",
                                  OrgDb = 'org.Hs.eg.db')

#print ("compare cluster done")
#save.image("compared_cluster.Rdata")
#print ("saved image")

# pdf(
#   file = file.path(out_directory, "dotplot_cp.pdf"),
#   width = 8.6,
#   height = 8.6
# )
#
# # clusterProfiler::dotplot(ck)
#
# dev.off()
#View(ck)
# ck@compareClusterResult$Description = NULL
fcc <- gofilter(ck, level = 3)
#View(cc)
#cc@compareClusterResult$Description
scc <- simplify(ck)

 pdf(
   file = file.path(out_directory, "enrichgo_cp.pdf"),
   width = 8.6,
   height = 8.6
 )

 fig13= clusterProfiler::dotplot(fcc)
 fig14=clusterProfiler::dotplot(scc)
fig13
fig14
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
temp = genetable[genetable$cluster == 0,]

for (i in c(0:8)) {
  temp = genetable[genetable$cluster == i, ]
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
   file = file.path(out_directory, "enrich_pw_run_go.pdf"),
   width = 8.6,
   height = 8.6
 )

fig15=dotplot(res)

dev.off()

######################code ends here###############################
### Figures code################
pdf(
   file = file.path(out_directory, "all_figures.pdf"),
   width = 8.6,
   height = 8.6
 )

fig1
fig2
fig3
fig4
fig5
fig6
fig7
fig8
fig9
fig10
fig11
fig12
fig13
fig14
fig15

dev.off()






