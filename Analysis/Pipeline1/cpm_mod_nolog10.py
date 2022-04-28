
import pandas as pd
import numpy as np
import argparse
import sys



def InputArgumentParser():
    inputParse = argparse.ArgumentParser(description='Input arguments and flags')
    inputParse.add_argument("inputFilename",help="File name of input data.txt file",type=str)
    Input = inputParse.parse_args()
    return Input

def readInputFile(filename): #opens and reads contents of file, outputing its content as list object
    #This assumes the required input files are already located in the same folder as this Python script
    #open files and create file object for it that has not been read yet, with only read permissions (so as to not modify the original file's content)
    fileObject=open(filename, "r") #opens file

    #read file object and put file content inside a list, where each row or line in the file is its own list element/item 
    #(transcribe each line in the file as an individual element that is part of a list)
    fileContent=fileObject.readlines() 

    #close file object after reading from it
    fileObject.close() #close  file object

    return(fileContent)


Input=InputArgumentParser() #capture, read and store arguments in variable

try: #attempt to load file in normal python way as well, as some functions cannot be done solely thorugh Biopython
    Content=readInputFile(Input.inputFilename) 
except:
    sys.exit("ERROR! File cannot be read or does not exist.") #force program to stop if file cannot be opened or does not exist.


Filename=str(Input.inputFilename)
outFilename=Filename.replace(".tsv","_cpm.csv")
fileOutObject=open(outFilename, "w") #opens file
outFilename2=Filename.replace(".tsv","_cpm3.csv")
fileOutObject2=open(outFilename2, "w") #opens file

#Remove unmapped counts from dataframe.
df = pd.read_csv(Filename, sep="\t", names=['Genes', 'Counts'])
df = df[df.Genes != '__no_feature']
df = df[df.Genes != '__ambiguous']
df = df[df.Genes != '__too_low_aQual']
df = df[df.Genes != '__not_aligned']
df = df[df.Genes != '__alignment_not_unique']


#Iterate over genes and list in order.
genes=[]
for i in df['Genes']:
    genes.append(i)
#print(genes)

#Iterate over counts and convert to counts per million (cpm) for values not equal to 0. Store in counts list.
Total_counts = df['Counts'].sum()
counts=[]
for i in df['Counts']:   
    if i ==0:
        i = i
        counts.append(i)
    elif i != 0:
        t=((i/Total_counts)*1000000)
        i=t.round()
        counts.append(i)
#print(counts)

#Reconstruct dataframe (using lists) with original gene names and ammended cpm. Write to file.
df2 = pd.DataFrame(counts,genes)
df2.to_csv(fileOutObject, index=True)

fileOutObject.close()

#Remove top row of file to account for first line error. 
df3 = pd.read_csv(outFilename, sep=",", skiprows=1)
df3.to_csv(fileOutObject2, index=False)
