library(Seurat)
library(sctransform)


He6hpf <- read.csv("counts_AJ/He6hpf_anno.csv", row.names=1)
He9hpf <- read.csv("counts_AJ/He9hpf_anno.csv", row.names=1)
He12hpf <- read.csv("counts_AJ/He12hpf_anno.csv", row.names=1)
He16hpf <- read.csv("counts_AJ/He16hpf_anno.csv", row.names=1)
He20hpf <- read.csv("counts_AJ/He20hpf_anno.csv", row.names=1)
He24hpf <- read.csv("counts_AJ/He24hpf_anno.csv", row.names=1)
He30hpf <- read.csv("counts_AJ/He30hpf_anno.csv", row.names=1)
He36hpf <- read.csv("counts_AJ/He36hpf_anno.csv", row.names=1)
He42hpf <- read.csv("counts_AJ/He42hpf_anno.csv", row.names=1)
He48hpf <- read.csv("counts_AJ/He48hpf_anno.csv", row.names=1)
He54hpf <- read.csv("counts_AJ/He54hpf_anno.csv", row.names=1)
He60hpf <- read.csv("counts_AJ/He60hpf_anno.csv", row.names=1)





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



png(file="violin_plots.png")
VlnPlot(He_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3, pt.size = 0)
VlnPlot(He_merged_full, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3, pt.size = 0)
dev.off()


##Subsetting the data for all samples

He_merged <- subset(He_merged, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000 )
dim(He_merged)

He_merged_full <- subset(He_merged_full, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000 )


dimensions <- seq(160,200, by=5)
vfeatures <- seq(3000,7000, by= 1000)

vfeatures
dimensions

for (x in vfeatures) {
	He_merged <- SCTransform(He_merged, verbose = T, variable.features.n = x,  vars.to.regress = "percent.Rb")
	He_merged <- RunPCA(He_merged,  npcs = 200, verbose = T)
	for (y in dimensions) { 
		print(x)
		print(y)
		He_merged <- RunUMAP(He_merged,  dims = 1:y)
		He_merged <- FindNeighbors(He_merged, dims = 1:y)
		He_merged <- FindClusters(He_merged, resolution = 1)
		Idents(He_merged) <- 'Stage'
		png(file=paste("He_6to30/dimplot_", x, "_vf_", y , "_dim.png", sep = ""))
			print(DimPlot(He_merged,  label=T)) 
		dev.off()
	}
}



for (x in vfeatures) {
        He_merged_full <- SCTransform(He_merged_full, verbose = T, variable.features.n = x,  vars.to.regress = "percent.Rb")
        He_merged_full <- RunPCA(He_merged_full,  npcs = 200, verbose = T)
        for (y in dimensions) {
                print(x)
                print(y)
                He_merged_full <- RunUMAP(He_merged_full,  dims = 1:y)
                He_merged_full <- FindNeighbors(He_merged_full, dims = 1:y)
                He_merged_full <- FindClusters(He_merged_full, resolution = 1)
                Idents(He_merged_full) <- 'Stage'
                png(file=paste("He_6to60/dimplot_", x, "_vf_", y , "_dim.png", sep = ""))
                        print(DimPlot(He_merged_full,  label=T))
                dev.off()
        }
}




print('saved file full')


sessionInfo()
