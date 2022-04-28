# SRP
Steered Research Project

# Directory Overview
*main.py*, *database.py* and *template_inserts.py* as well as **static** and **templates** folders, are for running the website.
1. 'main.py' is the python file that runs the webpage via Flask. 
2. 'database.py' is a file that was created to interact with the database through post-requests and SQL queries from the gene query webpage. 
3. 'template_inserts.py' is a file developed for text processing and directory searching for any required figures and texts required for a particular webpage.
4. The static folder is a folder that contains all the css style sheets, texts and figures required for a particular webpage 
5. The templates folder is a folder that contains all the html templates for a particular webpage. 


**Analysis** folder contains the scripts for processing and analysing data for both the original pipeline and the alternative pipeline.
1. The subfolders: pipeline1 and pipeline2 contain all the scripts required for pre-processing of the Darmanis data as well as clustering for both the original and alternative pipeline.

**database** folder contains database files for use within website's *genesearch* webpage
1. The database folder contains database files GeneQuery1 and GeneQuery2 that contains all the clusters and their associated genes obtained from the original and alternative pipelines. 

