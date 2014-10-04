from Bio import Phylo


tree = Phylo.read("genus_level.tree","newick")
#for node in tree.get_nonterminals():print node.name
for node in tree.find_clades():
    if node.name:
        print node.name
        
