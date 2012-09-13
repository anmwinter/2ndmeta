##################
#Region Extract
#Extracts specified regions from a CLUSTAL alignment (.aln) file
#Output a text file with sequence and species
#Output fasta file with the extracted sequences
##################

from Bio import AlignIO
from Bio import SeqIO

print "This program extracts a section of a set of aligned sequences"
print "You will need to know the start and stop positions before you begin"

print ""
print "Your filename should be something like: foo.aln"
file_in = raw_input("Please enter the name of your clustal alignment file:")

align = AlignIO.read(file_in,"clustal")

print "Alignment length %i" % align.get_alignment_length()

print "Primer Region 1"
print align[:,705:722]
SeqIO.write(align[:,705:722], "primer1_region.txt","tab")

out1 = open("primer1.fasta","wb")
AlignIO.write(align[:,705:722], out1,"fasta")
out1.close()

print ""
print "Primer Region 2"
print align[:,1083:1104]
SeqIO.write(align[:,1083:1104], "primer2_region.txt","tab")

out2 = open("primer2.fasta","wb")
AlignIO.write(align[:,1083:1104], out2,"fasta")
out2.close()

##print ""
##print "Hypervariable region 3"
##print align[:,457:480]
##SeqIO.write(align[:,457:480], "hv3_region.txt","tab")
##
##out3 = open("hv3.fasta","wb")
##AlignIO.write(align[:,457:480], out3,"fasta")
##out3.close()
##
##print""
##print "Hypervariable region 4"
##print align[:,837:853]
##SeqIO.write(align[:,837:853], "hv4_region.txt","tab")
##
##out4 = open("hv4.fasta","wb")
##AlignIO.write(align[:,837:853], out4,"fasta")
##out4.close()
##
##print""
##print "Hypervariable region 5"
##print align[:,1030:1076]
##SeqIO.write(align[:,1030:1076], "hv5_region.txt","tab")
##
##out5 = open("hv5.fasta","wb")
##AlignIO.write(align[:,1030:1076], out5,"fasta")
##out5.close()
##
##print""
##print "Hypervariable region 6"
##print align[:,1132:1154]
##SeqIO.write(align[:,1132:1154], "hv6_region.txt","tab")
##
##out6 = open("hv6.fasta","wb")
##AlignIO.write(align[:,1132:1154], out6,"fasta")
##out6.close()
##
##print""
##print "Hypervariable region 7"
##print align[:,1461:1473]
##SeqIO.write(align[:,1461:1473], "hv7_region.txt","tab")
##
##out7 = open("hv7.fasta","wb")
##AlignIO.write(align[:,1461:1473], out7,"fasta")
##out7.close()
