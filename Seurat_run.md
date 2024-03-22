
```bash
for dir in `cat count_tables/list`; do
cellranger mat2csv ${dir}/outs/filtered_feature_bc_matrix count_tables/${dir}_anno.csv
done
```

```bash
conda activate /data/wraycompute/alejo/conda/Seurat4_3
mkdir He_6to60mix
nano He_scAnalysis_6_to_60.R

```



```R
library(Seurat)
library(sctransform)


He_6_30hpf <- read.csv("count_tables/He_6_30hpf_anno.csv", row.names=1)
He_9_36hpf <- read.csv("count_tables/He_9_36hpf_anno.csv", row.names=1)
He_12_42hpf <- read.csv("count_tables/He_12_42hpf_anno.csv", row.names=1)
He_16_48hpf <- read.csv("count_tables/He_16_48hpf_anno.csv", row.names=1)
He_20_54hpf <- read.csv("count_tables/He_20_54hpf_anno.csv", row.names=1)
He_24_60hpf <- read.csv("count_tables/He_24_60hpf_anno.csv", row.names=1)


##Convert files to seurat format
He_6_30hpf_seurat <- CreateSeuratObject(counts = He_6_30hpf, project = "6_30hpf", min.cells = 3, min.features = 200)
He_9_36hpf_seurat <- CreateSeuratObject(counts = He_9_36hpf, project = "9_36hpf", min.cells = 3, min.features = 200)
He_12_42hpf_seurat <- CreateSeuratObject(counts = He_12_42hpf, project = "12_42hpf", min.cells = 3, min.features = 200)
He_16_48hpf_seurat <- CreateSeuratObject(counts = He_16_48hpf, project = "16_48hpf", min.cells = 3, min.features = 200)
He_20_54hpf_seurat <- CreateSeuratObject(counts = He_20_54hpf, project = "20_54hpf", min.cells = 3, min.features = 200)
He_24_60hpf_seurat <- CreateSeuratObject(counts = He_24_60hpf, project = "24_60hpf", min.cells = 3, min.features = 200)



He_6_30hpf_seurat$group <- "6_30hpf"
He_9_36hpf_seurat$group <- "9_36hpf"
He_12_42hpf_seurat$group <- "12_42hpf"
He_16_48hpf_seurat$group <- "16_48hpf"
He_20_54hpf_seurat$group <- "20_54hpf"
He_24_60hpf_seurat$group <- "24_60hpf"


##Merging the seurat objects
He_merged <- merge(He_6_30hpf_seurat, y = c(He_9_36hpf_seurat, He_12_42hpf_seurat, He_16_48hpf_seurat, He_20_54hpf_seurat, He_24_60hpf_seurat), 
                   add.cell.ids = c("6_30hpf", "9_36hpf", "12_42hpf", "16_48hpf", "20_54hpf", "24_60hpf"), project = "He_merged_series")

rm(He_6_30hpf_seurat, He_9_36hpf_seurat, He_12_42hpf_seurat, He_16_48hpf_seurat, He_20_54hpf_seurat, He_24_60hpf_seurat)

gc()



He_merged  <- PercentageFeatureSet(He_merged, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")


He_merged$Stage <- factor(He_merged$group, levels = c("6_30hpf", "9_36hpf", "12_42hpf", "16_48hpf", "20_54hpf", "24_60hpf"))

Idents(He_merged) <- "Stage"

png(file="violin_plots.png")
VlnPlot(He_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3, pt.size = 0)
dev.off()

##Subsetting the data for all samples

He_merged <- subset(He_merged, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 4000 )
dim(He_merged)

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
		png(file=paste("He_6to60mix/dimplot_", x, "_vf_", y , "_dim.png", sep = ""))
			print(DimPlot(He_merged,  label=T)) 
		dev.off()
	}
}


print('saved file full')


sessionInfo()


```
