cat He_scAnalysis_6_to_60.R
library(Seurat)
library(sctransform)

setwd("counts_AJ")

He6hpf <- read.csv("He6hpf_anno.csv", row.names=1)
He9hpf <- read.csv("He9hpf_anno.csv", row.names=1)
He12hpf <- read.csv("He12hpf_anno.csv", row.names=1)
He16hpf <- read.csv("He16hpf_anno.csv", row.names=1)
He20hpf <- read.csv("He20hpf_anno.csv", row.names=1)
He24hpf <- read.csv("He24hpf_anno.csv", row.names=1)
He30hpf <- read.csv("He30hpf_anno.csv", row.names=1)
He36hpf <- read.csv("He36hpf_anno.csv", row.names=1)
He42hpf <- read.csv("He42hpf_anno.csv", row.names=1)
He48hpf <- read.csv("He48hpf_anno.csv", row.names=1)
He54hpf <- read.csv("He54hpf_anno.csv", row.names=1)
He60hpf <- read.csv("He60hpf_anno.csv", row.names=1)





##Convert files to seurat format
He6hpf_seurat <- CreateSeuratObject(counts = He6hpf, project = "6hpf", min.cells = 3, min.features = 200)
He9hpf_seurat <- CreateSeuratObject(counts = He9hpf, project = "9hpf", min.cells = 3, min.features = 200)
He12hpf_seurat <- CreateSeuratObject(counts = He12hpf, project = "12hpf", min.cells = 3, min.features = 200)
He16hpf_seurat <- CreateSeuratObject(counts = He16hpf, project = "16hpf", min.cells = 3, min.features = 200)
He20hpf_seurat <- CreateSeuratObject(counts = He20hpf, project = "20hpf", min.cells = 3, min.features = 200)
He24hpf_seurat <- CreateSeuratObject(counts = He24hpf, project = "24hpf", min.cells = 3, min.features = 200)
He30hpf_seurat <- CreateSeuratObject(counts = He30hpf, project = "30hpf", min.cells = 3, min.features = 200)
He36hpf_seurat <- CreateSeuratObject(counts = He36hpf, project = "36hpf", min.cells = 3, min.features = 200)
He42hpf_seurat <- CreateSeuratObject(counts = He42hpf, project = "42hpf", min.cells = 3, min.features = 200)
He48hpf_seurat <- CreateSeuratObject(counts = He48hpf, project = "48hpf", min.cells = 3, min.features = 200)
He54hpf_seurat <- CreateSeuratObject(counts = He54hpf, project = "54hpf", min.cells = 3, min.features = 200)
He60hpf_seurat <- CreateSeuratObject(counts = He60hpf, project = "60hpf", min.cells = 3, min.features = 200)



He6hpf_seurat$group <- "6hpf"
He9hpf_seurat$group <- "9hpf"
He12hpf_seurat$group <- "12hpf"
He16hpf_seurat$group <- "16hpf"
He20hpf_seurat$group <- "20hpf"
He24hpf_seurat$group <- "24hpf"
He30hpf_seurat$group <- "30hpf"
He36hpf_seurat$group <- "36hpf"
He42hpf_seurat$group <- "42hpf"
He48hpf_seurat$group <- "48hpf"
He54hpf_seurat$group <- "54hpf"
He60hpf_seurat$group <- "60hpf"



##Merging the seurat objects
He_merged <- merge(He6hpf_seurat, y = c(He9hpf_seurat, He12hpf_seurat, He16hpf_seurat, He20hpf_seurat, He24hpf_seurat, He30hpf_seurat), 
                   add.cell.ids = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf","30hpf"), project = "He_merged_series")
He_merged_full <- merge(He6hpf_seurat, y = c(He9hpf_seurat, He12hpf_seurat, He16hpf_seurat, He20hpf_seurat, He24hpf_seurat, He30hpf_seurat,
                                             He36hpf_seurat, He42hpf_seurat, He48hpf_seurat, He54hpf_seurat, He60hpf_seurat)
                        , add.cell.ids = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf","30hpf","36hpf", "42hpf", "48hpf", "54hpf","60hpf"), project = "He_merged_series")

rm(He6hpf, He9hpf, He12hpf, He16hpf, He20hpf, He24hpf, He30hpf,He6hpf_seurat, He9hpf_seurat, He12hpf_seurat, He16hpf_seurat, He20hpf_seurat, He24hpf_seurat, He30hpf_seurat,
   He36hpf_seurat, He42hpf_seurat, He48hpf_seurat, He54hpf_seurat, He60hpf_seurat)

gc()



He_merged  <- PercentageFeatureSet(He_merged, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")
He_merged_full <- PercentageFeatureSet(He_merged_full, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")


He_merged$Stage <- factor(He_merged$group, levels = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf","30hpf"))
He_merged_full$Stage <- factor(He_merged_full$group, levels = c("6hpf", "9hpf", "12hpf", "16hpf", "20hpf", "24hpf","30hpf","36hpf", "42hpf", "48hpf", "54hpf" ,"60hpf"))



Idents(He_merged) <- "Stage"
Idents(He_merged_full) <- "Stage"



pdf(file="violin_plots.pdf")
VlnPlot(He_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3, pt.size = 0)
VlnPlot(He_merged_full, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3, pt.size = 0)
dev.off()


##Subsetting the data for all samples

He_merged <- subset(He_merged, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000 )
dim(He_merged)

He_merged_full <- subset(He_merged_full, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000 )
dim(He_merged_full)



######################### 2000 variable features


He_merged_transformed_1 <- SCTransform(He_merged, method = "glmGamPoi", verbose = T, variable.features.n = 2000 , vars.to.regress = "percent.Rb" )

# Run the standard workflow for dimension reduction
He_merged_transformed_1 <- RunPCA(He_merged_transformed_1, verbose = T, npcs =200)
He_merged_transformed_1 <- RunUMAP(He_merged_transformed_1,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_1 <- FindNeighbors(He_merged_transformed_1, dims = 1:150)
He_merged_transformed_1 <- FindClusters(He_merged_transformed_1, resolution = 0.5)

######################### 3000 variable features


He_merged_transformed_2 <- SCTransform(He_merged, method = "glmGamPoi", verbose = T, variable.features.n = 3000 , vars.to.regress = "percent.Rb" )

# Run the standard workflow for dimension reduction
He_merged_transformed_2 <- RunPCA(He_merged_transformed_2, verbose = T, npcs =200)
He_merged_transformed_2 <- RunUMAP(He_merged_transformed_2,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_2 <- FindNeighbors(He_merged_transformed_2, dims = 1:150)
He_merged_transformed_2 <- FindClusters(He_merged_transformed_2, resolution = 0.5)

######################### 4000 variable features


He_merged_transformed_3 <- SCTransform(He_merged, method = "glmGamPoi", verbose = T, variable.features.n = 4000 , vars.to.regress = "percent.Rb")

# Run the standard workflow for dimension reduction
He_merged_transformed_3 <- RunPCA(He_merged_transformed_3, verbose = T, npcs =200)
He_merged_transformed_3 <- RunUMAP(He_merged_transformed_3,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_3 <- FindNeighbors(He_merged_transformed_3, dims = 1:150)
He_merged_transformed_3 <- FindClusters(He_merged_transformed_3, resolution = 0.5)



save(He_merged_transformed_1, He_merged_transformed_2, He_merged_transformed_3 ,file="He_dataset_single_cell_times_match_half.Rda")



print('saved file 1')

#########################
####     FULL       #####
######################### 2000 variable features


He_merged_transformed_1 <- SCTransform(He_merged_full, method = "glmGamPoi", verbose = T, variable.features.n = 2000 , vars.to.regress = "percent.Rb")


# Run the standard workflow for dimension reduction
He_merged_transformed_1 <- RunPCA(He_merged_transformed_1, verbose = T, npcs =200)
He_merged_transformed_1 <- RunUMAP(He_merged_transformed_1,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_1 <- FindNeighbors(He_merged_transformed_1, dims = 1:150)
He_merged_transformed_1 <- FindClusters(He_merged_transformed_1, resolution = 0.5)

######################### 3000 variable features


He_merged_transformed_2 <- SCTransform(He_merged_full, method = "glmGamPoi", verbose = T, variable.features.n = 3000 ,  vars.to.regress = "percent.Rb")


# Run the standard workflow for dimension reduction
He_merged_transformed_2 <- RunPCA(He_merged_transformed_2, verbose = T, npcs =200)
He_merged_transformed_2 <- RunUMAP(He_merged_transformed_2,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_2 <- FindNeighbors(He_merged_transformed_2, dims = 1:150)
He_merged_transformed_2 <- FindClusters(He_merged_transformed_2, resolution = 0.5)


#########################   4000 variable features

He_merged_transformed_3 <- SCTransform(He_merged_full, method = "glmGamPoi", verbose = T, variable.features.n = 4000  , vars.to.regress = "percent.Rb")


# Run the standard workflow for dimension reduction
He_merged_transformed_3 <- RunPCA(He_merged_transformed_3, verbose = T, npcs =200)
He_merged_transformed_3 <- RunUMAP(He_merged_transformed_3,  dims = 1:150)
# t-SNE and Clustering
He_merged_transformed_3 <- FindNeighbors(He_merged_transformed_3, dims = 1:150)
He_merged_transformed_3 <- FindClusters(He_merged_transformed_3, resolution = 0.5)


save(He_merged_transformed_1, He_merged_transformed_2, He_merged_transformed_3 ,file="He_dataset_single_cell_times_match_full.Rda")

print('saved file 2')
