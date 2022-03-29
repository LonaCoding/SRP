#### Convert cluster ENSEMBLE IDs into gene symbols 
library(httr)
library(jsonlite)

url = "https://biotools.fr/human/ensembl_symbol_converter/"

#Get rownames 
E_id0 <- toJSON((row.names(cluster0)))
E_id1 <- toJSON(row.names(cluster1))
E_id2 <- toJSON(row.names(cluster2))
E_id3 <- toJSON(row.names(cluster3))
E_id4 <- toJSON(row.names(cluster4))
E_id5 <- toJSON(row.names(cluster5))

print(attributes(E_id0))
#Send GET requests 
#Cluster 0
body0 <- list(api=1, ids=E_id0)
response_0 <- POST(url, body=body0)
converted0 <- fromJSON(rawToChar(response_0$content),flatten=TRUE)
converted0 <- as.vector(converted0)

#Set rowname 
#nams0 <- make.names(converted0,unique=TRUE) 
#row.names(cluster0) <- nams0

#Cluster 1 
body1 <- list(api=1, ids=E_id1)
response_1 <- POST(url, body=body1)
converted1 <- fromJSON(rawToChar(response_1$content),flatten=TRUE)
converted1 <- as.vector(converted1)

#Cluster 2 
body2 <- list(api=1, ids=E_id2)
response_2 <- POST(url, body=body2)
converted2 <- fromJSON(rawToChar(response_2$content),flatten=TRUE)
converted2 <- as.vector(converted2)

#Cluster 3 
body3 <- list(api=1, ids=E_id3)
response_3 <- POST(url, body=body3)
converted3 <- fromJSON(rawToChar(response_3$content),flatten=TRUE)
converted3 <- as.vector(converted3)

#Cluster 4 
body4 <- list(api=1, ids=E_id4)
response_4 <- POST(url, body=body4)
converted4 <- fromJSON(rawToChar(response_4$content),flatten=TRUE)
converted4 <- as.vector(converted4)

#Cluster 5 
body5 <- list(api=1, ids=E_id5)
response_5 <- POST(url, body=body5)
converted5 <- fromJSON(rawToChar(response_5$content),flatten=TRUE)
converted5 <- as.vector(converted5)
