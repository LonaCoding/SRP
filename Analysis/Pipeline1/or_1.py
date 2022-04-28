#import re
import argparse
import sys

#script to convert identifier to match printseq requirements.
#usage example: python idAdaptor.py SRR1974543_1_prinseq_good_singletons_5_H4.fastq
#output: originalFilename_adp.fastq


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
outFilename=Filename.replace(".txt","_aseqs.txt")
fileOutObject=open(outFilename, "w") #opens file

found_type = False
t_ype = []

with open(Filename, 'r') as f:
    for line in f:
        if '>>Overrepresented sequences' in line:
            found_type = True
            continue                    

        if found_type:
            if '>>END_MODULE' in line:
                found_type = False               
            elif '#Sequence' not in line:    
                t_line = str.split(line)
                fileOutObject.write(">"+"\n"+(t_line[0])+"\n")

fileOutObject.close()
