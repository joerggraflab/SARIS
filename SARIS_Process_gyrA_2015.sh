#!/bin/sh

#  Process_gyrA.sh
#
#  A shell script for the processing of gyrase A amplicon data generated on the MiSeq.
#  Uses QIIME for aplicon clustering and taxonomy identification.
#  A custom perl script is used for determining presence of the gyrA mutation (S83I).
#  The input data must have been pre-processed correctly prior to running this script.
#
#  Created by Michael C. Nelson on 2013-09-16.
#  Version: 3
#  Last revised: 2015-04-6
#  Copyright (c) 2013 Michael C. Nelson/University of Connecticut. All rights reserved.

summary_table(){
    # Make a new file header by tacking on info abt resistance
    head -n 2 $1 | tail -n 1 >h1
    echo Resistance$'\t'Codon>h2
    paste h1 h2 > header
    # Sort the OTU data and paste on the resistance data
    tail -n +3 $1 | sort -k 1 -t $'\t' > data_sort
    cut -d $'\t' -f 2,3 $2 > res
    paste data_sort res > table
    cat header table > $3
    rm h1 h2 header data_sort res table
}

### Initialize timestamp variables
DATE=`date +%Y-%m-%d`
TIME=`date +%H:%M`
TM=`date +%Y%m%d-%H%M`
LOG=Process_gyrA_log_$TM.txt
### Create log file
echo "Executing Qiime_Process_gyrA.sh on $DATE at $TIME " | tee $LOG

# Set up working environment
clear
source /macqiime/configs/bash_profile.txt
echo ''

# Set script variables
#Paths to reference files
REF_SEQ=gyrA_ref.fasta
REF_TAX=gyrA_ref_taxonomy.txt
REF_ALN=gyrA_ref_aln.fasta
#Input file
INFILE=$1
INMD5=`md5 $INFILE`
FILENAME=${INFILE##*/}
BASENAME=${FILENAME%.*}
echo "Using $FILENAME as the input for processing." | tee -a $LOG
echo $INMD5 | tee -a $LOG

# Cluster the sequences based on 100% identity.
#echo '' | tee -a $LOG
#echo '' | tee -a $LOG
#echo 'Step 1: Clustering sequences into OTUs based on 99.7% identity' | tee -a $LOG
#pick_otus.py -i $INFILE -o otus/ -s 1.00 -m usearch61

# Pick representative sequences
#echo '' | tee -a $LOG
#echo 'Step 2: Picking representative sequences (Rep Seqs) from the full dataset' | tee -a $LOG
#echo 'Out is Rep_997.fasta' | tee -a $LOG
#pick_rep_set.py -i otus/"$BASENAME"_otus.txt -f $INFILE -o Rep_Seqs.fasta -m most_abundant

# Assign taxonomy to the Rep_Seqs using RDP classifier trained on gyrA reference sequences
#echo '' | tee -a $LOG
#echo 'Step 3: Assigning taxonomy to the Rep Seqs using BLAST against a set of gyrA reference sequences' | tee -a $LOG
#assign_taxonomy.py -i Rep_Seqs.fasta -o taxonomy_assignment/ -t $REF_TAX -r $REF_SEQ -m blast

# Align the sequences to make a tree to examine representatives that seem interesting
#echo '' | tee -a $LOG
#echo 'Step 4: Aligning the Rep Seqs  file against the aligned gyrA reference sequences' | tee -a $LOG
#align_seqs.py -i Rep_Seqs.fasta -o alignment/ -t $REF_ALN | tee -a $LOG

# Make the phylogenetic tree including the reference sequences
#echo '' | tee -a $LOG
#echo 'Step 5: Makign a phylogenetic tree of the Rep Seqs combined with the Reference Seqs' | tee -a $LOG
#cat alignment/Rep_Seqs_aligned.fasta $REF_ALN > Ref_plus_Rep.fasta
#make_phylogeny.py -i Ref_plus_Rep.fasta -o Ref_plus_Rep_tree.tree

# Create an OTU table for use in examining abundances
#echo '' | tee -a $LOG
#echo 'Step 6: Making an OTU table for the full dataset, output is raw_table.biom and raw_table.txt' | tee -a $LOG
#make_otu_table.py -i otus/"$BASENAME"_otus.txt -o raw_table.biom -t taxonomy_assignment/Rep_Seqs_tax_assignments.txt

# Filter out low abundance sequences as noise
#echo '' | tee -a $LOG
#echo 'Step 7: Filtering out singleton OTUs, as these most likely represent just sequencing noise' | tee -a $LOG
#echo 'The filtered table is table.biom' | tee -a $LOG
#filter_otus_from_otu_table.py -i raw_table.biom -o table.biom -n 2
#biom convert -i table.biom -o table.txt -b --header-key taxonomy
#biom summarize_table -i table.biom -o table_stats.txt

# Filter the aligned representative sequences and tree to just the reduced set
#echo '' | tee -a $LOG
#echo 'Step 8: Filtering out the low abundance sequences from the aligned Rep Seqs file, output is Rep_Seqs_aligned_filtered.fasta' | tee -a $LOG
#filter_fasta.py -f alignment/Rep_Seqs_aligned.fasta -o alignment/Rep_Seqs_aligned_filtered.fasta -b table.biom

# Determine Cipro sensitivity/resistance profiles based on the sequence for amino acid position 84 and output an overall summary table
#echo '' | tee -a $LOG
#echo 'Step 9: Checking for Cipro resistance of the Rep_Seqs_aligned_filtered.fasta sequences' | tee -a $LOG
#perl gyrA_Resistance.pl alignment/Rep_Seqs_aligned_filtered.fasta
#echo 'Creating an overall summary table including OTU abundances, taxonomy, and Cipro resistance: Results_Summary.txt'
#summary_table table.txt alignment/Rep_Seqs_aligned_filtered_resistance.txt Results_Summary.txt

# Calculate alpha diversity measures
#echo '' | tee -a $LOG
#echo 'Step 10: Calculating alpha diversity measures from the filtered table' | tee -a $LOG
#alpha_diversity.py -i table.biom -o alpha.txt -m osd,simpson,shannon

# Rarefy the OTU table, recalculate alpha diversity, calculate beta diversity measures, and output an overall results summary table.
#DEPTH=`grep 'Min:' table_stats.txt | sed 's/ Min: //' | sed 's/\.0//'`
DEPTH=5000
echo '' | tee -a $LOG
echo "Step 11: Rarefying the filtered table to $DEPTH sequences" | tee -a $LOG
single_rarefaction.py -i table.biom -o table_even.biom -d $DEPTH
biom convert -i table_even.biom -o table_even.txt -b --header-key taxonomy
echo 'Step 12: Creating an overall summary table' | tee -a $LOG
filter_fasta.py -f alignment/Rep_Seqs_aligned_filtered.fasta -o alignment/Rep_Seqs_aligned_filtered_even.fasta -b table_even.biom
perl gyrA_Resistance.pl alignment/Rep_Seqs_aligned_filtered_even.fasta
summary_table table_even.txt alignment/Rep_Seqs_aligned_filtered_even_resistance.txt Rarefied_Results_Summary.txt
echo "Recalculating alpha diversity measures" | tee -a $LOG
alpha_diversity.py -i table_even.biom -o alpha_even.txt -m osd,simpson,shannon
echo "Calculating beta diversity measures and creating 3D PCoA plots using Bray-Curtis and weigthed UniFrac distances" | tee -a $LOG
beta_diversity.py -i table_even.biom -o beta/ -m bray_curtis,weighted_unifrac -t Ref_plus_Rep_tree.tree
principal_coordinates.py -i beta/bray_curtis_table_even.txt -o beta/bray_curtis_PCoA.txt
principal_coordinates.py -i beta/weighted_unifrac_table_even.txt -o beta/weighted_unifrac_PCoA.txt
make_3d_plots.py -i beta/bray_curtis_PCoA.txt -o beta/BC/ -m gyrA_Map.txt
make_3d_plots.py -i beta/weighted_unifrac_PCoA.txt -o beta/wUF/ -m gyrA_Map.txt
mkdir beta/3D
find ./ -name '*.kin' -exec mv {} beta_div/3D/ ';' ; rm -r beta/BC beta/wUF

