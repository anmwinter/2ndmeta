from Bio import SeqIO
 
def genus_level(fasta_file):

    sequences={}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence=seq_record.id
        if sequence not in sequences:
            sequences[sequence]=str(seq_record.seq)

    output_file=open("clear_"+fasta_file,"w+")
    out_genera=open("genera_list"+fasta_file,"w+")

    for sequence in sequences:
            output_file.write(">"+sequence+"\n"+sequences[sequence]+"\n")
            out_genera.write(sequence+"\n")
    output_file.close()
    out_genera.close()

genus_level("actino_refseq_cem.fasta")
#Call the f(x) like this:
#genus_level("my_fasta.fasta")
