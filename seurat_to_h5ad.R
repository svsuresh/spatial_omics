## Both lines below works to load anndata library.
Sys.setenv(RETICULATE_PYTHON = "/home/suresh/.mc3/bin/python3.8")
#import_from_path("anndata", path = "/home/suresh/.mc3/lib/python3.8/site-packages/")

## Load the data
setwd("~/npg_work/DBiT-seq_FFPE/st_results_04052021/")
df1=read.csv("text_outputs/Filtered_matrix.tsv", head=T, sep="\t", row.names=1)
df2=t(df1)
mdf2=Matrix::Matrix(as.matrix(df2), sparse=T)
library(Seurat)
smdf2=CreateSeuratObject(mdf2, min.cells=10, project="test")
smdf2 <- PercentageFeatureSet(smdf2, pattern = "^MT-", col.name = "percent.mt")
smdf2 <- SCTransform(smdf2, vars.to.regress = "percent.mt", verbose = FALSE)
smdf2 <- RunPCA(smdf2, verbose = FALSE)
smdf2 <- RunUMAP(smdf2, dims = 1:10, verbose = FALSE)
smdf2 <- FindNeighbors(smdf2, dims = 1:10, verbose = FALSE)
smdf2 <- FindClusters(smdf2, resolution = 0.8, verbose = FALSE)
ls()
library(sceasy)
library(reticulate)
setwd("~/Desktop/")

convertFormat(smdf2, from="seurat", to="anndata",outFile="test.h5ad")


## Output is produced and but is not loaded in cellxgene.
anndata <- anndata <- import("anndata", convert = FALSE)
test=anndata$AnnData(
    X = t(as.matrix(GetAssayData(object = smdf2))),
    obs = data.frame(smdf2@meta.data),
    obsm  = list(
        "X_pca" = Embeddings(smdf2[["pca"]]),
        "X_umap" = Embeddings(smdf2[["umap"]])
    )
)
anndata$AnnData$write(test, "filename1.h5ad")
dir()
