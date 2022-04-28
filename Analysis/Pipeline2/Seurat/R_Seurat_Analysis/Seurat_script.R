## Custom pipeline w/ darmanis data 
# import modules 
library(Seurat)
library(dplyr)
library(sctransform)
library(ggplot2)

######DATA IMPORT / FORMATTING 
raw.data<-read.table('featurecounts.txt',header=T,sep='\t') #table containing genes as rows, cells as columns, and a few other variables as rows (chromosome, etc. )
#replace ensembl IDs with gene IDs 
id.data <- read.table('gene_id_conversion2.txt',header=T,sep='\t',na.strings = c("")) #import ensembl id to gene symbol data 
raw.data<- raw.data[,!names(raw.data) %in% c('Length','Chr','Start','End','Strand','Geneid')] #drop columns that do not correspond to cells 
colnames(raw.data)<- gsub('Aligned.sortedByCoord.out.bam','',colnames(raw.data)) #strip annoying suffix off the end of cell names 
rows<- make.names(id.data$Gene_id,unique=TRUE) #make gene ids unique so they can be used as row names 
rownames(raw.data)<-rows

###### SEURAT ANALYSIS
sro <- CreateSeuratObject(counts=raw.data,project='darmanis',min.cells=23,min.features=800) #create seurat analysis object
#above removes cells w/ < 500 genes and genes seen in 5% of cells
sro

#ADD CELL TYPES FROM META DATA 
celltypes<- read.table('celltype_metadata.csv',header=T,sep=',',row.names=1)
sro<- AddMetaData(sro,celltypes,col.name='biased_type') # add biased cell types to the metadata 
sro@meta.data
sro

#get list of cell names and types to see which are dropped 
cells_preqc<- sro$biased_type

## QC
sro[['percent.mt']] <- PercentageFeatureSet(sro,pattern="^MT.")
sro$percent.mt #HIGHER MITOCHONDRIAL GENE COVERAGE in a lot of cells (! - note in results)

#visualise metdata (for supplementary figures) - violin plots 
VlnPlot(sro,features=c('nFeature_RNA','nCount_RNA','percent.mt'),ncol=3) #COULD COLOUR THESE BY CLUSTER HERE 

# scatter plot of feature count vs molecule count 
FeatureScatter(sro,feature1='nCount_RNA',feature2='nFeature_RNA')
# scatter plot of mt % vs molecule count 
FeatureScatter(sro,feature1='nCount_RNA',feature2='percent.mt')

#filter the data with high mitochondrial % 
sro<- subset(sro,subset=nFeature_RNA > 1000 & nFeature_RNA <10000 & percent.mt < 10 ) 
#no more than 10k genes, no less than 1k, no more than 10% mt taken into final sample 
sro@meta.data
cells_afterqc<- sro$biased_type

#get list of cell types dropped 
drop_cells<-c()

for (x in names(cells_preqc)){
  if (!x %in% names(cells_afterqc)){drop_cells[x] <- cells_preqc[x]}
}
table(drop_cells) #amounts of the different types of cells dropped 

genes_per_cell <- sro$nFeature_RNA
mean(genes_per_cell) #get the mean no. of genes across all 466 cells (for report)
#read counts per cell 
molecules_per_cell <- sro$nCount_RNA
mean(molecules_per_cell) # get mean AND min() no. of molecules 


# plot to show variance of features 
plot1<- VariableFeaturePlot(sro)
plot2<- LabelPoints(plot=plot1,points=top10,repel=T)
plot1
plot2

## data scaling 
#all.genes<- rownames(sro)
#sro<- ScaleData(sro,features=all.genes) #scale expression so mean = 0, variance = 1 for all cells 
## NB: COULD TRY regression out mitochondrial feature influence here

## SCTRANSFORM - for normalisation, variable feature identification (?), and gene scaling 
# (does all 3 - see SCTransform() doc for details)
sro<- SCTransform(sro,vars.to.regress = 'percent.mt',verbose=T)
sro


# plot to show variance of features 
plot1<- VariableFeaturePlot(sro)
plot2<- LabelPoints(plot=plot1,points=top10,repel=T)
plot1
plot2

#look at variable features 
top10<- head(VariableFeatures(sro),10) #find the top 10 most variable features 

# plot to show variance of features 
plot1<- VariableFeaturePlot(sro)
plot2<- LabelPoints(plot=plot1,points=top10,repel=T)
plot1
plot2

## DIMENSIONALITY REDUCTION - PCA + plots
sro<- RunPCA(sro,features=VariableFeatures(sro))

#plot showing gene contributions to PC 
VizDimLoadings(sro,dims= 1:5,reduction='pca')
#PCA scatterplot 
DimPlot(sro,reduction='pca')
#Heatmap of PCA components / gene 
DimHeatmap(sro,dims=5,balanced=T)

#Jackstraw plot - to determine optimal no. of PCs to include 
#sro<- JackStraw(sro)
#sro<- ScoreJackStraw(sro,dims = 1:20)

#JackStrawPlot(sro,dims=1:15) #suggests no less than 4 PCs, use 7+

###CLUSTERING OF CELLS 
sro<- FindNeighbors(sro, dims=1:30) #construct KNN graph based on distances in PCA space, SCTransform means many PCs can be used
sro<- FindClusters(sro,resolution=0.4) #determine cell clusters based on graph. resolution should be lower for smaller data sets (might need to be even lower here)
Idents(sro) #get each cell's cluster label 

#rename idents and seurat clusters, starting labels from 1 
Idents(sro)<- as.factor(as.numeric(as.character(Idents(sro)))+1)
#rename clusters so the indexing begins from 1, not 0 
sro$seurat_clusters <- as.factor(as.numeric(as.character(sro$seurat_clusters)) + 1)

#save number of cells assigned to each cluster 
clusts <- Idents(sro)
table(clusts)
write.csv(clusts,'../cluster_sizes.csv')

biased<- sro$biased_type
table(biased)
# umap plot for clusters 
sro<- RunUMAP(sro,dims=1:30)

biasedumap<- DimPlot(sro,reduction='umap',group.by='biased_type',label=T) + labs(title = 'UMAP w/ biased cell types')#get 2d umap of clusters 
clusterumap<- DimPlot(sro,reduction='umap',group.by='seurat_clusters',label=T) + labs(title= 'UMAP w/ Seurat unbiased clusters')
biasedumap + clusterumap

###FINDING DIFFERENTIALLY EXPRESSED FEATURES OF CLUSTERS 
#get lists of the top markers for all clusters 
cluster1.markers<- FindMarkers(sro,ident.1=1,min.pct = 0.25)
cluster2.markers<- FindMarkers(sro,ident.1=2,min.pct = 0.25)
cluster3.markers<- FindMarkers(sro,ident.1=3,min.pct = 0.25)
cluster4.markers<- FindMarkers(sro,ident.1=4,min.pct=0.25)
cluster5.markers<- FindMarkers(sro,ident.1=5,min.pct=0.25)
cluster6.markers<- FindMarkers(sro,ident.1=6,min.pct=0.25)

#get the top 2 positive markers for each cluster 
sro.markers <- FindAllMarkers(sro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#EXPORT markers 
#do below for all markers or write loop 
clust_export<-cluster6.markers[order(cluster6.markers[,2],decreasing = T),]
clust_export<- head(clust_export,20)
write.csv(clust_export,'../markers/allmarkers_cluster6.csv')

#display the top 2 markers of each cluster using dplyr
sro.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#heatmap of genes by cluster (unbiased groupings)
sro.markers %>%
  group_by(cluster) %>%
  top_n(n=10,wt=avg_log2FC) -> top10

png('../figures/clustered_heatmap.png')
DoHeatmap(sro,features=top10$gene,group.by = 'seurat_clusters') + NoLegend() #use this to determine cell types, annotating manually
DoHeatmap(sro,features=top10$gene,group.by='biased_type') + NoLegend()
dev.off()
##3D PLOTS USING UMAP DATA 
##
library(rgl)
umap3d<- RunUMAP(sro,dims=1:30,n.components=3L) # copy seurat object with 3d umap (not really that efficient but whatever)

#get data for 3D umap 
plot.data<- FetchData(umap3d,vars=c('UMAP_1','UMAP_2','UMAP_3','seurat_clusters','biased_type'))

#3d plot with our unbiased clusters 
#get colours function for consistent colours in plot 
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
} # need this for colours (taken from stackoverflow)
plot1cols<- get_colors(plot.data$seurat_clusters)
plot1cols
plot3d(plot.data$UMAP_1,plot.data$UMAP_2,plot.data$UMAP_3,col=plot1cols, xlab='UMAP1',ylab='UMAP2',zlab='UMAP3',decorate = T,size=4)
leg1col<- c('black','#DF536B','#61D04F','#2297E6','#28E2E5','#CD0BBC') # get colours for plot
bgplot3d({
  plot.new()
  title(main='3D UMAP w/ unsupervised Seurat clusters')
  legend("topright", legend= levels(plot.data$seurat_clusters),col=leg1col, pch = 16, cex=1, inset=c(0.02))
})

#same as above but for biased type clusters 
plot2cols<- get_colors(plot.data$biased_type)
split(plot2cols,names(plot2cols))['OPC'] <- 'purple' #replace black (used twice for some reason) with orange-y colour
plot3d(plot.data$UMAP_1,plot.data$UMAP_2,plot.data$UMAP_3,col=plot2cols,xlab='UMAP1',ylab='UMAP2',zlab='UMAP3',decorate=T,size=4)
leg2col<- unique(plot2cols)
names(leg2col)<- unique(names(plot2cols))
bgplot3d({
  plot.new()
  title(main='3D UMAP w/ biased cell type assignments')
  legend("topright", legend= unique(plot.data$biased_type),col=leg2col, pch = 16, cex=1, inset=c(0.02))
})


#save 3d plots 
snapshot3d(filename='../figures/biased_3dumap.png',fmt='png',webshot = F)

##MARKER PLOTS 
## to show marker representation in each cluster (choose good markers)
#plots showing DE levels of particular genes over clusters 
# NB: do for genes representing clusters in darmanis 
#show violin plots (or box/whisker) of expression of various genes FOR REPORT 
#using best markers identified in our clusters AND unbiased groups of darmanis 
VlnPlot(sro,features= c('CD24','VIP','AQP4','ERMN','PDGFRA','ITGAX'),group.by = 'biased_type') 

#heatmap of gene expression through clusters 
FeaturePlot(sro,features=c(''))
