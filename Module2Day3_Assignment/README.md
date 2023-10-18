#  CPBS Module2 Day3 Assignment 

## Overview
The python script generates two populations of subnetworks and performs permutation test:<br>
- Generate 5000 random subnetworks with disease associated genes (one gene from each loci)
- Generate 5000 random non-informative subnetworks 
- Run non-parametric test (Permutation test)
- Compute empirical p-value

## Input Arguments

There are two input files:<br> 
- -i or --input: Disease genes file in Gene Map Table (GMT) format from OMIM <br>
- -n or --network: Protein protein interactions network, for example, from STRING database. Column1 is gene1, column2 is gene2, column3 is edge<br>

Arguments:
- nSub: Number of subnetworks. Default is 5000 <br>
- nBins: Number of bins for edge counts. Default is 500 <br>
- nPerm: Number of permutations. Default is 10000 <br>

## Output files

- FA_subnetworks.txt
- Null_subnetworks.txt
- FA_subnetwork_distribution.png 
- Null_subnetwork_Distribution.png 
- Permutation_statistic.png

## Usage

`python Subnetworks_Significance.py Disease_gene_file.txt Protein_Protein_Interaction.txt --nSub --nBins --nPerm`

Example:

`python Subnetworks.py Input.gmt.txt STRING 1.txt --nSub 5000 --nBins 500 --nPerm 10000`

## Additional Documentation

https://docs.google.com/document/d/1LAkTVxXKhVa0joW9RsNr0VPWyctVbuGfIAhhJcwRDg0/edit?usp=sharing
