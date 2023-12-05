## Import modules and packages
import sys  
import math
import matplotlib.pyplot as plt
import networkx as nx
from netgraph import Graph
import argparse
import time
import random
import mplcursors
import community
from matplotlib.patches import Patch
from GA_Functions import data_to_diction,string_to_diction,first_parent_population,parent_mean_edgeweight,generate_random_excluding,ga_mutation,ga_mating
from GA_Functions import gene_frequency,bin_dictionary,match_fa_nonfa_by_bucket,generate_nonfa_from_fa,generate_random_loci_NonFA
from GA_Functions import gene_score,gene_scores_avg

# Run time start
start = time.time()

string_path = "./STRING.txt"
input_path = "./Input.gmt.txt"
num_iterations=5000
min_edges=0
num_bins=100
loci_diction,fa_set,fa_diction,fa_string_set,string_set,all_genes_set = data_to_diction(input_path,string_path)
string_gene_set,string_diction = string_to_diction(string_path)

parent_subnetwork,subnetwork_diction,subnetwork_diction_edge_counts=first_parent_population(loci_diction,fa_diction,num_iterations=num_iterations,min_edges=min_edges)
first_population=parent_subnetwork
normalized_weights_parent,subnetwork_diction_updated,parent_subnetwork,subnetwork_normalized_weighted=ga_mutation(loci_diction,fa_diction,parent_subnetwork,num_iterations=num_iterations)
parent_subnet_edgeweight,parent_mean_edgeweights = parent_mean_edgeweight(subnetwork_diction,num_iterations=num_iterations)
mate_population,mate_parent_subnetwork,mate_weighted,mean_delta,mean_delta_perc=ga_mating(fa_diction,normalized_weights_parent,parent_mean_edgeweights,parent_subnetwork,num_iterations=num_iterations,min_edges=min_edges)

parent_ = parent_subnetwork
mate_population_opt=mate_population
gen=1
FA_stats=["Generation","Parent_Mean_edgeweight","Mate_Mean_edgeweight","Mean_difference","Percentage_difference"]
while abs(mean_delta_perc) > 0.5:

      normalized_weights_parent_opt,subnetwork_diction_updated_opt,parent_subnetwork_opt,subnetwork_normalized_weighted_opt=ga_mutation(loci_diction,fa_diction,parent_subnetwork=parent_,num_iterations=num_iterations)
      parent_subnet_edgeweight,parent_mean_edgeweights_opt = parent_mean_edgeweight(mate_population_opt,num_iterations=num_iterations)
      mate_population_opt,mate_parent_subnetwork_opt,mate_weighted_opt,mean_delta,mean_delta_perc=ga_mating(fa_diction,normalized_weights_parent_opt,parent_mean_edgeweights_opt,parent_subnetwork=parent_,num_iterations=num_iterations,min_edges=min_edges)
      print(gen,parent_mean_edgeweights_opt,sum(mate_weighted_opt.values())/num_iterations,mean_delta,mean_delta_perc)
      FA_stats.append([gen,parent_mean_edgeweights_opt,sum(mate_weighted_opt.values())/num_iterations,mean_delta,mean_delta_perc])
      parent_ = mate_parent_subnetwork_opt
      fa_pop_opt=mate_population_opt
      final_sub=mate_parent_subnetwork_opt
      gen+=1

with open('FA_subnetwork_Opt.txt', 'w') as file:
    for k, v in fa_pop_opt.items():
        for sub in v:
            file.write("{}\t{}\n".format(k,'\t'.join(sub)))

scores=[float(vals) for vals in subnetwork_normalized_weighted.values()]
plt.hist(scores, color='steelblue', edgecolor='black', linewidth=1.2,bins=10)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Normalized density score')
plt.ylabel('Frequency')
plt.title('Distribution of normalized density in the Initial FA population (Mutation)')
plt.savefig('Firstpop_Mutation_Normalized_weights.png')
plt.close()

scores=[float(vals) for vals in subnetwork_diction_edge_counts.values()]
plt.hist(scores, color='steelblue', edgecolor='black', linewidth=1.2,bins=10)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('average density score')
plt.ylabel('Frequency')
plt.title('Distribution of average density score in the Initial FA population')
plt.savefig('Firstpop_avg_density.png')
plt.close()

scores=[float(vals) for vals in mate_weighted.values()]
plt.hist(scores, color='steelblue', edgecolor='black', linewidth=1.2,bins=10)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('average density score')
plt.ylabel('Frequency')
plt.title('Distribution of average density score in the First mating population')
plt.savefig('First_matingpop_avg_density.png')
plt.close()


gene_freq=gene_frequency(string_diction,all_genes_set)
bin_diction,filtered_bin_diction,empty_keys=bin_dictionary(gene_freq,fa_set,num_bins=num_bins,filter_bin_data=True)

nonFA_opt_pops={}
nonFA_opt_pops_mean={}
for i in range(1000):
    non_fa_loci=[]
    nonFA_parent,nonFA_parent_size,nonFA_parent_genes=generate_nonfa_from_fa(optimized_FA_subnetwork=parent_subnetwork,random_seed=i,bin_diction=bin_diction,string_diction=string_diction)
    non_fa_loci=generate_random_loci_NonFA(nonFA_parent_genes,all_genes_set)
    nonFA_normalized_weights_parent,nonFA_subnetwork_diction_updated,nonFA_parent_subnetwork,nonFA_subnetwork_normalized_weighted=ga_mutation(non_fa_loci,string_diction,parent_subnetwork=nonFA_parent_genes,num_iterations=num_iterations)
    NonFA_parent_mean_edgeweights=parent_mean_edgeweight(nonFA_parent,num_iterations=num_iterations)
    nonFA_mate_population,nonFA_mate_parent_subnetwork,nonFA_mate_weighted,nonFA_mean_delta,nonFA_mean_delta_perc=ga_mating(string_diction,nonFA_normalized_weights_parent,NonFA_parent_mean_edgeweights,parent_subnetwork=nonFA_parent_genes,num_iterations=num_iterations,min_edges=min_edges)
    parent_ = nonFA_parent_genes
    nonFA_opt_pops["population_"+str(i+1)]=[]
    gen=1
    nonFA_parent_opt=nonFA_parent
    while abs(nonFA_mean_delta_perc) > 0.5:
        nonFA_normalized_weights_parent_opt,nonFA_subnetwork_diction_updated_opt,nonFA_parent_subnetwork_opt,nonFA_subnetwork_normalized_weighted_opt=ga_mutation(non_fa_loci,string_diction,parent_subnetwork=parent_,num_iterations=num_iterations)
        NonFA_parent_subnet_edgeweight,NonFA_parent_mean_edgeweights_opt=parent_mean_edgeweight(nonFA_parent_opt,num_iterations=num_iterations)
        nonFA_mate_population_opt,nonFA_mate_parent_subnetwork_opt,nonFA_mate_weighted_opt,nonFA_mean_delta,nonFA_mean_delta_perc=ga_mating(string_diction,nonFA_normalized_weights_parent_opt,NonFA_parent_mean_edgeweights_opt,parent_subnetwork=parent_,num_iterations=num_iterations,min_edges=min_edges)
        print(gen,NonFA_parent_mean_edgeweights_opt,sum(nonFA_mate_weighted_opt.values())/num_iterations,nonFA_mean_delta,nonFA_mean_delta_perc)
        parent_ = mate_parent_subnetwork_opt
        non_fa_loci=generate_random_loci_NonFA(parent_,all_genes_set)
        nonFA_parent_opt=nonFA_mate_population_opt
        nonFA_opt_pops["population_"+str(i+1)].append([gen,NonFA_parent_mean_edgeweights_opt,sum(nonFA_mate_weighted_opt.values())/num_iterations,nonFA_mean_delta,nonFA_mean_delta_perc])
        gen+=1
    
    nonFA_opt_pops_mean["population_"+str(i+1)]=sum(nonFA_mate_weighted_opt.values())/num_iterations

for sub,vals in parent_subnet_edgeweight.items():
    subnet_pval={}
    for pop, mean in nonFA_opt_pops_mean.items():
        sum=0
        if mean>=vals:
            sum+=1
    subnet_pval[sub]=sum/1000
    
gene_scores_loci = {loci: {inner_key: 0 for inner_key in loci_genes} for loci, loci_genes in loci_diction.items()}
for subnetwork,locus_genes in parent_.items():
    iteration_fa_genes_diction_2={'Locus_'+str(ind):loc_gene for ind,loc_gene in enumerate(locus_genes)}
    gene_score_func = gene_score(iteration_fa_genes_diction_2,fa_diction,loci_diction)

    for loci_num,inner_diction in gene_score_func.items():
          for loci_gene,score in inner_diction.items():
              gene_scores_loci[loci_num][loci_gene]=gene_scores_loci[loci_num].get(loci_gene,0)+score

    gene_scores_final=gene_scores_avg(loci_diction,gene_scores_loci,fa_string_set,num_iterations=num_iterations) 

loci_data={loci: [] for loci in loci_diction.keys()}
for locus_num,genes_scores in gene_scores_final.items():
    for genes,scores in genes_scores.items():
        loci_data[locus_num].append(scores)

with open('Day3_Output.gmt.txt', 'w') as file:
    for locus, genes in gene_scores_final.items():
        file.write(f"{locus} ")
        for gene, value in genes.items():
            file.write(f"{gene} {value} ")
        file.write('\n')

hist_data=[]
for scores in loci_data.values():
    for score in scores:
        if score != 'NA':
            hist_data.append(float(score))

plt.hist(hist_data,color='steelblue', edgecolor='black', linewidth=1.2,bins=10)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Gene score')
plt.ylabel('Frequency')
plt.title('Distribution of Gene Scores')
plt.savefig('GeneScores_Distribution.png')
plt.close()

boxplot_data = [float(scores) for loci_data, scores in loci_data.items() if 'NA' not in scores]
 #Generate boxplot
plt.boxplot(boxplot_data, labels=loci_data.keys())
plt.xticks(rotation=90)
plt.xlabel('Loci Number')
plt.ylabel('Scores')
plt.title('Scores per Loci')
plt.savefig('GeneScores_Distribution.png')
plt.close()

# run time
end = (time.time())
print(end-start)

