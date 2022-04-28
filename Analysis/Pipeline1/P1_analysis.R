################# Pipeline 1 Analysis ####################


########## package installation ###########


#install remotes using remotes and devtools to aid installation of packages
install.packages("remotes")
install.packages("devtools")


# install scde dependency 'flexmix' specific version 2.3-13 (release 2015-01-17) to ensure the combination of both packages does not throw errors when generating models; https://github.com/hms-dbmi/scde/issues/40
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

# install scde dependency 'Cairo' version 1.5-15 (release date 16-03-2022); must be installed prior to scde, later release used as early releases threw errors during compiling;
# "ERROR: compilation failed for package ‘Cairo’
install.packages("Cairo")

# boot version 1.3-11 (release date 28-03-2014) is required during reciprocal weighted distance adjustment in scde
download.file(url = "https://cran.r-project.org/src/contrib/Archive/boot/boot_1.3-11.tar.gz", destfile = "boot_1.3-11.tar.gz")
install.packages("boot_1.3-11.tar.gz", repos = NULL, type = 'source', dependencies = TRUE)

## early version of scde (release date 10/02/2016); this actually installs version 1.99.1
##in order to solve problem with this package and a bug with flexmix, this version (and not the previous) had to be used alongside flexmix; https://github.com/hms-dbmi/scde/issues/40
download.file(url = "https://github.com/hms-dbmi/scde/archive/1.99.2.tar.gz", destfile = "1.99.2.tar.gz")
install.packages("1.99.2.tar.gz", repos = NULL, type = 'source', dependencies = TRUE)

# early version of Rtsne (release date 26-05-2015); version 0.10 is first version to have distance matrix handling function.
##Visne package with tsne implementation originally used but this package was stuck behind a paywall, so opted to use Rtsne instead from a similar release date.
download.file(url = "https://cran.r-project.org/src/contrib/Archive/Rtsne/Rtsne_0.10.tar.gz", destfile = "Rtsne_0.10.tar.gz")
install.packages("Rtsne_0.10.tar.gz", repos = NULL, type = 'source', dependencies = TRUE)

# installation of mclust 4.4 (release date 16-09-2014)
install_version("mclust","4.4")

install.packages("scatterplot3d") # Install for graphics; up to date version will not affect analysis results.

# parallel should be installed as a part of R version 4.1.2


#load packages
library(scde)
library(mclust) 
library(Rtsne)
library(zoo)
library(boot)
library(parallel)
library(scatterplot3d)
library(dplyr)

# set working directory and save all outputs to this directory alongside the matrix generated


######### Matrix generation ###########

#combine the counts per million transformed HTSeq outputs for each cell into a single matrix
#CPM_out_nolog10 directory will have to be saved to the working directory

#get all counts filenames starting with "GSM" (or "GSM165") using pattern matching and put into list
inputFileList<-list.files(path="/home/sesm2/SRA_data/CPM_out_nolog10",pattern=NULL, all.files=FALSE, full.names=FALSE)
print(inputFileList)

#all input files has same first column, which acts as index.
#rename 2nd column's header to match source filename
#remove redundant first column
#then, append to bigger dataframe horizontally (by column)
#end result is a single dataframe containing all concatenated data from multiple csv files


NormCounts_df<-data.frame(nrow=56643) #initialize empty 22088x0 dataframe for receiving input files

first=1 #set first input dataframe flag to TRUE
inputFileNameList<-inputFileList #use list of file names gained from previous cell
for (f in inputFileNameList){
  fn<-paste("/home/sesm2/SRA_data/CPM_out_nolog10/", f, sep ="") #reconstruct full relative path and filename with extension already included
  fs<-sub("_cpm3.csv$", "", f) #remove extension from file name
  fi<-paste(fs,"_Index")
  #print(fn) #displays input filenames
  
  df2<-read.csv(fn, header = FALSE, stringsAsFactors = FALSE, sep = ",")
  
  names(df2)[1]<- fi #rename first column's header
  names(df2)[2]<- fs #rename second column's header
  
  if (first==0){ #If current input dataframe is not the first in entire dataframe
    df2[1]<-NULL #drop first column, because it's redundant (deactivate this line to perform all-column append to see if all rows align)
  }
  
  NormCounts_df<-cbind(NormCounts_df,df2) #append new df
  first=0 #no longer the first dataframe, set flag to FALSE
}
NormCounts_df[1]<-NULL #drop first column (not used. final cleanup)
  
NormCounts_df<-rename(NormCounts_df, 'Gene' = 'SRR1974543 _Index') #rename first column to 'Gene'.

write.table(NormCounts_df, file = "P1_Matrix_nolog10.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


######### SCDE ##########

# scde will be used to ascertain pairwise distances as described in;
# Kharchenko PV, Silberstein L, Scadden DT (2014) Bayesian approach to single-cell differential expression analysis. Nat Meth 11:740–742.

### formatting ###

mat<- read.csv('P1_Matrix_nolog10.txt',header=T, sep= "\t") #import the file
rownames(mat) <- mat[,1] #assign rownames as column 1 
mat[,1]<- NULL #remove column 1
t_mat<-t(mat)
#count matrix is now properly formatted 

### build model ###

# generate model, no groups required; linear fit advised in later publications but not prior to 2015 so have used original parameters
o.ifm <- scde.error.models(counts = mat, n.cores =16, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="scde_fit1.RData") # save model to file

# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm, counts = mat) # This must be done prior to reciprochal weighting.

### Adjusted distance meaures ###

#Reciprochal weighting; taken from scde website tutorial https://hms-dbmi.github.io/scde/diffexp.html

# Reciprocal weight used over direct dropout as this method performed better in the original publication referenced and with the GSE data.
# load boot package for the weighted correlation implementation

require(boot)
cell.names <- colnames(mat); names(cell.names) <- cell.names;
k <- 0.95; # this parameter is recommended in both the tutorial and the cited paper
reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = o.ifm[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = o.ifm[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(mat[, nam1], mat[, nam2])+1), w = pnf)
  }))
},mc.cores = 1)), upper = FALSE)


########## Rtsne ##########

#Rtsne will be used to implement tsne dimensionality reduction in place of visne as originally described in;
#Amir E-AD et al. (2013) viSNE enables visualization of high dimensionalsingle-cell data and reveals phenotypic heterogeneity of leukemia. Nature Biotechnology 31:545–552.

set.seed(42) #f or reproducibility, returns same cluster number and cells per cluster everytime
tsne_model_1<- Rtsne(reciprocal.dist, verbose = 1, is_distance = TRUE, dims=3, theta=0, perplexity = 20) # added theta=0.0 for accuracy 
d_tsne_1 = as.data.frame(tsne_model_1$Y) # the reduced dimensions are held in this list
save(d_tsne_1,file="tnse_out_final.RData")

########## Mclust ##########

#Mclust was used to test for all EM algorithm heirarchical clustering initialised Gaussian mixture models (implemented via Bayesian information criterion) and obtain clusters as prevoiusly described;
#Fraley C, Raftery AE (1999) MCLUST: Software for Model-Based Cluster Analysis. J of Classification 16:297–306.

library(mclust)
fit<-Mclust(d_tsne_1, G=1:20) # run default parameters to test all models using BIC  
summary(fit) #check cluster number and distribution of cells in between clusters
plot(fit) #user interactive, must press escape in the terminal to move onto the next section; select 1 for BIC plot.

########## SCDE; differential expression ##########

#SCDE will be used to test for differential expression of each cluster against the remaining cell population with a view to provide cluster markers.
#The package or statistical method used in the original publication was not defined, so scde was used again as a likely candidate due to its previous use for cell distance recalculation.
#This section can handle increased usage of cores (I used n=8) as long as they are dedicated to the job only so will work as a job submission on the cluster.
#Rather than write a function, this section was written out and applied in separate scripts that could run one cluster_id at a time on the HPC, otherwise running time would be too long.

### Cluster1 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=1] <- "RMD"  
cluster_id[cluster_id==1] <- "CL1"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL1).*", "\\1", colnames(mat)), levels = c("RMD", "CL1"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus1_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL1).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL1"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster1_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c1<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)

### Cluster2 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=2] <- "RMD"  
cluster_id[cluster_id==2] <- "CL2"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL2).*", "\\1", colnames(mat)), levels = c("RMD", "CL2"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus2_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL2).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL2"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster2_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c2<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)

### Cluster3 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=3] <- "RMD"  
cluster_id[cluster_id==3] <- "CL3"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL3).*", "\\1", colnames(mat)), levels = c("RMD", "CL3"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus3_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL3).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL3"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster3_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c3<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)


### Cluster4 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=4] <- "RMD"  
cluster_id[cluster_id==4] <- "CL4"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL4).*", "\\1", colnames(mat)), levels = c("RMD", "CL4"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus4_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL4).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL4"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster4_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c4<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)



### Cluster5 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=5] <- "RMD"  
cluster_id[cluster_id==5] <- "CL5"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL5).*", "\\1", colnames(mat)), levels = c("RMD", "CL5"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus5_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL5).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL5"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster5_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c5<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)


### Cluster6 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=6] <- "RMD"  
cluster_id[cluster_id==6] <- "CL6"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL6).*", "\\1", colnames(mat)), levels = c("RMD", "CL6"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus6_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL6).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL6"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster6_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c6<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)


### Cluster7 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=7] <- "RMD"  
cluster_id[cluster_id==7] <- "CL7"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL7).*", "\\1", colnames(mat)), levels = c("RMD", "CL7"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus7_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL7).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL7"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster7_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c7<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)

### Cluster8 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=8] <- "RMD"  
cluster_id[cluster_id==8] <- "CL8"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL8).*", "\\1", colnames(mat)), levels = c("RMD", "CL8"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus8_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL8).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL8"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster8_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c8<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)


### Cluster9 ###

# take the list of cluster allocations from the mclust fit and generate a list that defines if the position of the entry is in the cluster to be tested or the remaining cell population
cluster_id<-fit$'fit.classification'
cluster_id[cluster_id!=9] <- "RMD"  
cluster_id[cluster_id==9] <- "CL9"  
table(cluster_id)

# create a metadata dataframe to combine both cell names and their allocation to either the cluster being tested or the remaining dataset. 
metadata<-data.frame(rbind(cluster_id, colnames(mat)))
metadata<-rbind(metadata, row3 = apply(metadata, 2, paste0, collapse = "-"))
colnames(metadata) <- metadata[3,] #assign row3 as colnames
colnames(mat)<-colnames(metadata) #transfer adapted cell names to the original counts matrix so that factoring can be applied to separate both groups for differential experssion

# define two groups of cells
sg <- factor(gsub("(RMD|CL9).*", "\\1", colnames(mat)), levels = c("RMD", "CL9"))
names(sg) <- colnames(mat)  
table(sg)

# calculate models
o.ifm <- scde.error.models(counts = mat, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
save(o.ifm,file="clus9_oifm.RData")

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = mat, length.out = 466, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(RMD|CL9).*", "\\1", rownames(o.ifm)), levels = c("RMD", "CL9"))
names(groups) <- colnames(mat)  
table(groups)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, mat, o.prior, groups = groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

ediff2<-ediff[order(ediff$Z, decreasing  =  TRUE), ]
# write out a table with all the results, showing most significantly different genes (in both directions) on top. Converted back to p-values
write.csv(ediff2, "cluster9_genes.csv")

#filter for significant genes
ediff2<- filter(ediff2, ediff2$p.adjust<0.05) 

# top 20 upregulated genes
top_20_c9<-head(ediff2[order(ediff2$mle, decreasing  =  TRUE), ], n=20)


############ Figures ############

t_mat<-t(mat)#transpose matrix for use in metadata labeling
fit_class<-as.data.frame(fit$classification) #take cluster classification to be added to metadata
SRR_index <-as.data.frame(rownames(t_mat)) # index of cell cluster metadata must match Rtsne output
cell_clusters <- cbind(SRR_index,fit_class,d_tsne_1) #cell cluster metadata

##### multiqc plots

#number of reads per cell
png(file="/home/sesm2/SRA_data/reads_per_cell.png", width=600, height=350)
STAR_out<-read.csv("STAR_totals_multiqc_star.txt", sep="\t", header=T)
hist(STAR_out$total_reads, main='Number of reads per cell')
dev.off()

#fraction of mapped reads
png(file="/home/sesm2/SRA_data/mapped_frac.png", width=600, height=350)
STAR_out$mapped_frac<-(STAR_out$uniquely_mapped_percent/100)
hist(STAR_out$mapped_frac, main='Fraction of mapped reads')
dev.off()

#gene body coverage
png(file="/home/sesm2/SRA_data/gene_body_cov.png", width=600, height=350)
HTSeq_out<-read.csv("HTSeq_multiqc_general_stats.txt", sep="\t", header=T)
HTSeq_out$gbc<-(HTSeq_out$HTSeq.Count_mqc.generalstats.htseq_count.percent_assigned/100)
hist(HTSeq_out$gbc, main='Gene body coverage', breaks=12)
dev.off()


##### reproduce uncertainty plot
png(file="/home/sesm2/SRA_data/uncertainty_plot.png", width=800, height=450)
box=cbind((as.data.frame(fit$classification)), (as.data.frame(fit$uncertainty)))
grid(nx=x, ny=y) #grid over boxplot
boxplot(fit$uncertainty ~ fit$classification, data = box, xlab= "Cluster", ylab= "Clustering uncertainty", theme_bw())
dev.off()

##### denovo cluster plot
rownames(d_tsne_1)<-rownames(t_mat)
rownames(cell_clusters)<-rownames(t_mat)
lab_cluster<-cell_clusters


png(file="/home/sesm2/SRA_data/denovo_cluster.png", width=600, height=350)
lab_cluster$legend <- lab_cluster$`fit$classification`
lab_cluster$legend = paste0("Cluster ", lab_cluster$legend)
colors <- c("#FB9A99", "#E31A1C","green4","dodgerblue2","#FF7F00", "black", "gold1","skyblue2", "#6A3D9A", "palegreen2")
colors <- colors[as.numeric(lab_cluster$`fit$classification`)]

x <- lab_cluster$V3
y <- lab_cluster$V1
z <- lab_cluster$V2
scatterplot3d(x,y,z, angle = 35, pch = 20, color = colors, grid=TRUE, box=FALSE, xlab="", ylab="", zlab="", mar=c(2,8,0,9)+0.1)
leg<-factor(lab_cluster$legend)
legend("right", legend = levels(leg),
       col = c("#FB9A99", "#E31A1C","green4","dodgerblue2","#FF7F00", "black", "gold1","skyblue2", "#6A3D9A", "palegreen2"), pch = 16, cex=0.85, inset = 0.0000000000001)
dev.o

###### adult vs foetal plot
# check clusters against adult/foetal distribution as in text; import metadata
dev_metadata<- read.csv('dev.csv', sep=',',header=T) #import the file 
cell_clusters_dev <- cbind(cell_clusters,dev_metadata$Dev)

png(file="/home/sesm2/SRA_data/adultvsfoetal.png", width=600, height=350)
lab_cluster<-cbind(lab_cluster, dev_metadata$Dev)
colors <- c("black", "red")
colors <- colors[as.numeric(lab_cluster$`dev_metadata$Dev`)]

x <- lab_cluster$V3
y <- lab_cluster$V1
z <- lab_cluster$V2
scatterplot3d(x,y,z, angle = 35, pch = 20, color = colors, grid=TRUE, box=FALSE, xlab="", ylab="", zlab="", mar=c(2,8,0,9)+0.1)

lab_cluster$legend <- lab_cluster$`fit$classification`
lab_cluster$legend = paste0("Cluster ", lab_cluster$legend)
leg<-factor(dev_metadata$dev_legend)
legend("right", legend = levels(leg),
       col = c("black", "red"), pch = 16, cex=0.85, inset = 0.0000000000001)
dev.off()

##### plot against original paper's cell-type assignment

dev_metadata<- read.csv('dev.csv', sep=',',header=T) #import the file 
cell_clusters_dev <- cbind(cell_clusters,dev_metadata$cell_type)


png(file="/home/sesm2/SRA_data/cell_type.png", width=600, height=350)
lab_cluster<-cbind(lab_cluster, dev_metadata$cell_type)
colors <- c("#FB9A99", "#E31A1C","green4","dodgerblue2","#FF7F00", "black", "gold1","skyblue2", "#6A3D9A", "palegreen2","brown")
colors <- colors[as.numeric(lab_cluster$`dev_metadata$cell_type`)]

x <- lab_cluster$V3
y <- lab_cluster$V1
z <- lab_cluster$V2
scatterplot3d(x,y,z, angle = 35, pch = 20, color = colors, grid=TRUE, box=FALSE, xlab="", ylab="", zlab="", mar=c(2,8,0,9)+0.1)

lab_cluster$legend <- lab_cluster$`dev_metadata$cell_type`
lab_cluster$legend = lab_cluster$`dev_metadata$cell_type`
leg<-factor(dev_metadata$Cell_type)
legend("right", legend = levels(leg),
       col = c("#FB9A99", "#E31A1C","green4","dodgerblue2","#FF7F00", "black", "gold1","skyblue2", "#6A3D9A", "palegreen2","brown"), pch = 16, cex=0.7, inset = 0.00000001)
dev.off()

##### clusters table

# a table consisting of cluster numbers and the distribution of cells taken directly from the mclust output
summary(fit)
cluster_n<-c("31","169","22","52","65","38","15","56","18")
cluster_label<-c(1,2,3,4,5,6,7,8,9)
cluster_tab<-cbind(cluster_label,cluster_n)
colnames(cluster_tab)<-c("Cluster Number","Number of Cells")
write.csv(cluster_tab, "cluster_table.csv")

##### top 20 enriched genes table
c1<-rownames(top_20_c1)
c2<-rownames(top_20_c2)
c3<-rownames(top_20_c3)
c4<-rownames(top_20_c4)
c5<-rownames(top_20_c5)
c6<-rownames(top_20_c6)
c7<-rownames(top_20_c7)
c8<-rownames(top_20_c8)
c9<-rownames(top_20_c9)

top_20_tab<-cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9)
colnames(top_20_tab)<-c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")
write.csv(top_20_tab, "top_20_tab.csv")
