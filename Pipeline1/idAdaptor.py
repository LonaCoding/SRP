#import re
import argparse
import sys

#script to convert identifier to match printseq requirements.
#usage example: python idAdaptor.py SRR1974543_1_prinseq_good_singletons_5_H4.fastq
#output: originalFilename_adp.fastq


def InputArgumentParser():
    inputParse = argparse.ArgumentParser(description='Input arguments and flags')
    inputParse.add_argument("inputFilename",help="File name of input PDB file",type=str)
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


def replacement(input):
    outsub=input.replace(".1 ","\\1 ")
    out=outsub.replace(".2 ","\\2 ")
    return out


Input=InputArgumentParser() #capture, read and store arguments in variable

try: #attempt to load file in normal python way as well, as some functions cannot be done solely thorugh Biopython
    Content=readInputFile(Input.inputFilename) 
except:
    sys.exit("ERROR! File cannot be read or does not exist.") #force program to stop if file cannot be opened or does not exist.


Filename=str(Input.inputFilename)
outFilename=Filename.replace(".fastq","_adp.fastq")

fileOutObject=open(outFilename, "w") #opens file
c=5 #line counter
for l in Content: #1,4,7
    if c==5:
        out=replacement(l)
        c=1
    else:
        out=l

    fileOutObject.write(out)  
    #fileOutObject.write("\n")

    c=c+1

fileOutObject.close()