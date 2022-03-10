#required packages:
#scde 3.0 (from bioconductor) 
#CytobankAPI/tsne (for viSNE)
#cyt (for visualizing viSNE) ????
#mclust 4.4
#igraph 0.7.1
#FactoMineR 1.28

#install packages using remotes (with specific version)
install.packages("remotes")  #get and install remotes first
library(remotes) #load remotes before everything else
install_bioc("3.0/scde",upgrade="never") #installing from bioconductor github page
install_version("scde","3.0") #not using the bioconductor install function
install_version("mclust","4.4")
install_version("igraph","0.7.1")
install_version("FactoMineR","1.28")

#install packages without using remotes (no specific version)
install.packages("bioconductor") #??
install.packages("scde") 
install.packages("CytobankAPI") #guaranteed and has ViSNE subpackage, but not sure if it's the one asked by paper (for ViSNE)
install.packages("tsne") #valid package, mentioned in paper but may not be right one because no subpackages (for ViSNE)
install.packages("Rtsne")
install.packages("mclust")
install.packages("igraph")
install.packages("FactoMineR")

#load packages
library(CytobankAPI)
library(bioconductor)
library(scde)
library(tsne)
library(Rtsne)
library(mclust)
library(igraph)
library(FactoMineR)
library(tidyverse) #for tsne
library(knitr) #for outputting and CytobankAPI



#load data


#0. excluded cells that had less than 400000 reads, 
#   reducing the initial dataset of 482 cells to 466


#1. scde: Pairwise distances between cells


#2. viSNE <- tsne: Dimensionality reduction of the distances
#is the source package tsne or Rtsne or cytobankAPI???

##Rtsne method:
#remove missing data
#make sure unique row ID exists

##cytobankAPI method:
#cytobankAPI needs username and password for aunthethication to site???

#endpoint.method -> viSNE , dimensionality_reduction

#copied from https://rdrr.io/cran/CytobankAPI/f/vignettes/cytobank-quickstart.Rmd section: "Making a request"

#Requests to Cytobank API endpoints can only be made after authenticating (see above). 
#The authentication object will be passed as the first parameter of any endpoint request.

# Making a request to the 'experiments' endpoint to create a new experiment
new_experiment <- experiments.new(cyto_session, experiment_name="My New Experiment", purpose="CytobankAPI quickstart vignette")

View(new_experiment)

c1 <- c(1, 0.007, 0, 0.0023, 0, 0.0026)

new_experiment_dataframe <- data.frame(id=c(1),
                                       version=c(42),
                                       purpose=c("CytobankAPI quickstart vignette"),
                                       comments=c(""),
                                       public=c(FALSE),
                                       deleted=c(FALSE),
                                       sources=c(""),
                                       experimentName=c("My New Experiment"),
                                       ...=c("..."),
                                       stringsAsFactors=FALSE)

knitr::kable(new_experiment_dataframe)


#Generating viSNE maps: (see Supplementary Table 2)
#   2a. uniformly subsample between 6,000 and 12,000 cells


#   2b. ViSNE (500 iterations) -> project data to 2D
#       All runs used:
#           1. an identical random seed 
#           2. default t-SNE parameters: 
#               1. perplexity = 30
#               2. momentum = ...
#                   a. 0.5 for initial 250 iterations
#                   b. 0.8 for remaining (250) iterations
#               3. epsilon = 500
#               4. lie factor = ...
#                   a. 4 for initial 100 iterations
#                   b. 1 for remaining (400) iterations 


#   2c. cyt: visualize viSNE maps and figures 
#            (color coding by subtype, by marker expression levels and in plotting expression level densities).
#            non-existent???


#3. mclust: Subsequent clustering 
#           w/ Bayesian information criterion (BIC) for parameterized Gaussian mixture models fitted by EM algorithm initialized by model-based hierarchical clustering


#4. igraph: Minimum spanning trees, community identification and computation of longest paths 


#5. FactoMineR: Principal components analysis (PCA)

