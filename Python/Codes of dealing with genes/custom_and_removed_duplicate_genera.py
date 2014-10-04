from Bio import SeqIO
import re
import sys

################################################################################
#
#Here we start with a function to access the data in the fasta file
#and to create a custom fasta that reads like:
#> some_genus_name
#sequence
#
################################################################################
def custom_fasta(fasta_file):
   
    #out_list = open("list"+fasta_file+".txt","w+")
    out_custom = open("custom_fasta.fasta","w+")

    
    for seq_record in SeqIO.parse(fasta_file,"fasta"):

#r is a regular expression. Here I am searching for a white space (\s) followed
#by a capital letter and then lowercase letters
        genus_name = extract(seq_record.description, r'\s[A-Z][a-z]+(.*)16S')
        #print genus_name

#r is a regular expression. Here I am searching for a some number of capital
#letters and underscore [A-Z_] plus some numbers \d+, followed by a period,
#and some numbers \d+. This will give you the accession number
        acc_num = extract(seq_record.description, r'[A-Z_]+\d+.\d+')
        #print acc_num
        genus = str(genus_name)
        
        
        
        out_custom.write(">"+genus.strip()+" "+acc_num+"\n"+str(seq_record.seq)+"\n")

    out_custom.close()

################################################################################
#
#The code that follows takes your custom fasta file and removes duplicates 
#by name. You will end up with a fasta file at the genus level. There is a
#seperate file call genera_list which lists each genera and it's accesion number
#
################################################################################
    sequences={}
    for seq_record in SeqIO.parse("custom_fasta.fasta", "fasta"):
        sequence=extract(seq_record.description,r'[A-Z][a-z]+(.*)16S')
        #print sequence
        acc = extract(seq_record.description, r'[A-Z_]+\d+.\d+')
        if sequence not in sequences:
            sequences[sequence]=str(seq_record.seq)

    output_file=open("singletons_"+fasta_file,"w+")
    out_genera=open("genera_list"+fasta_file+".txt","w+")

    for sequence in sequences:
            output_file.write(">"+sequence+"\n"+sequences[sequence]+"\n")
            out_genera.write(sequence+"\t"+acc+"\n")
    output_file.close()
    out_genera.close()
        
#####
#
#This function sets up the reguar expression routine for pulling out information
#in the fasta file. This will convert the extracted pattern to a string
#
#####
def extract(src, pattern):
     src_str = str(src)
     matchObj = re.search(pattern, src_str)
     if matchObj:
          return matchObj.group()
     else:
          return "NONE"


filename = sys.argv[-1]
custom_fasta(filename)

