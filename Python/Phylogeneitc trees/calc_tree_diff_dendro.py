#!/usr/bin/python                                                                          import dendropy
from dendropy import tree_source_iter 

best_tree = dendropy.Tree.get_from_path("best_tree.tre", schema="nexus")

output_f = open('treedist_output', 'w')
counter = 1

for tree in tree_source_iter(
        stream=open('all_trees', 'rU'),
        schema='nexus',
        taxon_set=best_tree.taxon_set):

    distance = best_tree.symmetric_difference(tree)
    output_f.write("1 " + str(counter) + " " + str(distance) + "\n")
    counter += 1

output_f.close()
