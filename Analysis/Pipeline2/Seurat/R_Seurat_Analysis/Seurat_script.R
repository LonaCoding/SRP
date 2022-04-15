library(SeuratObject)
library(Seurat)
raw_counts <- read.table('featurecounts.txt',header=T,stringsAsFactors = FALSE, sep='\t')
counts <- raw_counts[,-(2:6)]
rownames(counts) <- counts[,1]
counts <- as.matrix(counts)
counts <- counts[,-1]
counts <- CreateSeuratObject(counts = counts)

#Normalise data 
norm_counts <- NormalizeData(counts, normalization.method = 'RC', scale.factor = 1e6)


#Highly variable features 
variable_features <- FindVariableFeatures(norm_counts, selection.method = 'vst')

#Data scaling 
genes <- rownames(variable_features)
sn_data <- ScaleData(variable_features, features = genes)

#Linear dimensional reduction through PCA 
pca_data <- RunPCA(sn_data, features = VariableFeatures(sn_data))

#Visualise PCA 
map <- DimHeatmap(pca_data,dims=1, cells = 500, balanced = TRUE)
dimplot <- DimPlot(pca_data,reduction='pca')
dimplot

#Elbow plot to determine dimensions of the data 
ElbowPlot(pca_data,ndims=50)

#JackStrawPlot 
jack_pca <- JackStraw(pca_data, num.replicate = 200, dims = 50)
jack_pca <- ScoreJackStraw(jack_pca, dims=1:50)
JackStrawPlot(jack_pca, dims = 10:30) #13?

#Cluster the cells, chosen dims is 15 
cluster_data <- FindNeighbors(pca_data,dims=1:15)
cluster_data <- FindClusters(cluster_data,resolution=0.5)
head(Idents(cluster_data),5)


# tSNE 
tSNE_data <- RunTSNE(cluster_data, dims=1:15)
DimPlot(tSNE_data, reduction='tsne')

#Defining markers for clusters 
cluster0.markers <- FindMarkers(tSNE_data, ident.1 = 0, ident.2=c(1:5), min.pct=0.25)
cluster1 <- FindMarkers(tSNE_data, ident.1 = 1, ident.2 = c(1,2:5), min.pct = 0.25)
cluster2 <- FindMarkers(tSNE_data, ident.1 = 2, ident.2 = c(0,1,3,4,5), min.pct=0.25 )
cluster3 <- FindMarkers(tSNE_data, ident.1 = 3, ident.2 = c(0,1,2,4,5), min.pct=0.25)
cluster4 <- FindMarkers(tSNE_data, ident.1 = 4, ident.2 = c(0,1,2,3,5), min.pct(0.25))
cluster5 <- FindMarkers(tSNE_data, ident.1 = 5, ident.2 = c(0,1,2,3,4),min.pct(0.25))

#Get markers and order 

