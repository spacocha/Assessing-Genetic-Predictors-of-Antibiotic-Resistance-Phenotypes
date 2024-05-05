#!/bin/bash
# written by vappiah: https://github.com/vappiah/create-phylogenetic-trees/blob/main/generate_phylo.sh

MSA_phylip=$1
bname=$(basename $MSA_phylip)
apps/phyml -i $MSA_phylip -m JC69 -o tlr

mv ${bname}_phyml_stats.txt phyml_stats.txt 
mv ${bname}_phyml_tree.txt phyml_tree.txt
