###
#Hacking together a script for blasting a fasta file of bacterial 16S which I forgot to capture all the meta data.
#This takes my fuckup of a fasta file and submits it to NCBI and pulls down the top 5 hits for each sequence.
#Hopefully the xml file then goes through a filter that only finds a 100% identity match and then dumps
#out all the information. 
#
#If you are seeing this comment you are looking at the test file which isn't writing anything to a file
#just to the screen.
###



from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq




file_in = open("clear_genus_level.fasta")
fas_rec = SeqIO.parse(file_in,"fasta")
first_seq = fas_rec.next().seq
result = NCBIWWW.qblast("blastn", "nr", first_seq, expect=0.001, hitlist_size=5, megablast=True)


record = NCBIXML.parse(result)
blast_record = record.next()

print "Species, Accession, Length, Sequence"
for alignment in blast_record.alignments:
   
    for hsp in alignment.hsps:
        
        percent_id = (100*hsp.identities)/hsp.align_length
        if hsp.gaps == 0 and percent_id == 100:
           # print "***FOUND***"
           # print "Species:", alignment.title
            title_element = alignment.title.split()
           # print "Species:", title_element[1]+" "+title_element[2]
           # print "Accession:", alignment.accession
           # print "Length:", alignment.length
           # print "gaps:", hsp.gaps
           # print "expect:", hsp.expect
           # print "%id:", percent_id
           # print hsp.sbjct
            print  title_element[1]+" "+title_element[2]+","+" "+alignment.accession+","+" "+str(alignment.length)+","+" "+hsp.sbjct




