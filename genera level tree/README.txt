Towards a Genera Level Bacterial Tree for Overlaying Data
Journy of a Crypto-geologist in a Phylogentic World

FASTA FILES
clear_genus_level.fasta

Contains 2228 (with four Archeal outgroups) genera each with one 16S representative sequence. The core of the tree is built with 16S sequences mined from the NCBI RefSeq database. I gathered the remaining sequences from NCBI(1) and greengenes(2). 

genera_list_genus.txt

Contains a list of all the genera names that were included in this study.

ALIGNMENT FILES
There are five alignment files to go with each tree that was built. 
muscle_align.fasta 		COMMANDS -maxiters 2
infernal_align.fasta 		COMMANDS default Pyro Pipeline
clustalo_align.fasta 		COMMANDS --full
nastmask_align.fasta 		COMMANDS
infernal_mask_align.fasta 	COMMANDS 



TREES
There five trees: muscle.tree, infernal.tree, clustalo.tree, nast_mask.tree, infernal_mask.tree. Each tree is from the same initial fasta file but aligned with different software as describe above.

Worst delta-Log as report by FastTree
muscle		65.593
infernal		
clustalo	19.747		
nast
infernal mask

PDFS



References:
(1) Geer LY, Marchler-Bauer A, Geer RC, Han L, He J, He S, Liu C, Shi W, Bryant SH. The NCBI BioSystems database. Nucleic Acids Res.

(2) Jeffrey J Werner, Omry Koren, Philip Hugenholtz, Todd Z DeSantis, William A Walters, J Gregory Caporaso, Largus T Angenent, Rob Knight and Ruth E Ley,Impact of training sets on classification of high-throughput bacterial 16S rRNA gene surveys. 2011, ISME J. 


