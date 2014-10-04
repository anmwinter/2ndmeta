from Bio import SeqIO

# barcode_reads = (rec for rec in \
#                 SeqIO.parse("SRR342214.fastq", "fastq") \
#                 if rec.seq.startswith("GTCAGAGA"))
# count = SeqIO.write(barcode_reads, "AB_KS.fastq", "fastq")
# print "Saved %i reads" % count


from Bio import SeqIO
trimmed_primer_reads = (rec[8:] for rec in \
                       SeqIO.parse("SRR342214.fastq", "fastq") \
                       if rec.seq.startswith("AACCTACG"))
count = SeqIO.write(trimmed_primer_reads, "UT_KSa_wo_barcode.fastq", "fastq")
print "Saved %i reads" % count

SeqIO.convert("UT_KSa_wo_barcode.fastq", "fastq", "UT_KSa.qual", "qual")
SeqIO.convert("UT_KSa_wo_barcode.fastq", "fastq", "UT_KSa.fasta", "fasta")