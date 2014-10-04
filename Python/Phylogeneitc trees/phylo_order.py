import csv
import collections

#from Bio import Phylo
#tree = Phylo.read("complete.tree", "newick")
#for leaf in tree.get_terminals(): print leaf.name



index = collections.defaultdict(list)

out_writer = csv.writer(open('finsihed.csv','wb'))

file1= open( "data.csv", "rb" )
rdr = csv.DictReader( file1 )
for row in rdr:
    index[row['Genus']].append( row )
file1.close()

file2= open( "ordered_from_tree.csv", "rb" )
rdr= csv.DictReader( file2 )
for row in rdr:
    data = row, index[row['Genus']]
    out_writer.writerow(data)
file2.close()



inp = file("finished.csv","r")
outp = open("final_data.csv","w")

for line in inp:
   line = line.replace('"','')
   line = line.replace('{','')
   line = line.replace('}','')
   line = line.replace(']','')
   line = line.replace('[','')
   line = line.replace('\'','')
   outp.write(line)

inp.close()
outp.close()
