#---------------------------------------------------------------------------------#
#- Generate 5000 random subnetworks with one gene from each loci                 -#
#- Generate 5000 random subnetworks of non-informative genes (Null Distribution) -#
#- Run Permutation test for statistical significance                             -#
#- Compute Empirical p-value                                                     -#
#---------------------------------------------------------------------------------#             

import matplotlib.pyplot as plt
import sys
import argparse
import time

start = time.time()

# Function to convert the string file to a nested dictionary and set of all the genes from the string file
def string_to_diction(string_path):
    string_diction={}
    string_gene_set = set()
    with open(string_path,'r') as string_file:
        for edge in string_file:
            genes = sorted(edge.strip().split('\t'),reverse=True)  #Remove duplicates
            string_gene_set.update(genes[:2])
            if genes[0] in string_diction:
                if ((genes[1] not in string_diction[genes[0]].keys()) and genes[1] not in string_diction[genes[0]].keys()):
                    string_diction[genes[0]][genes[1]]=genes[2]
            elif genes[0] not in string_diction:
                string_diction[genes[0]] = {genes[1]:genes[2]}
    return string_gene_set,string_diction    #string_gene_set has all the unique gene names from the string_diction

# Function to convert input file into loci dictionary, set of loci, filtered FA dictionary from string dictionary
def input_to_diction(input_path,string_path):
    with open(input_path,'r') as input:  
        loci_diction = {'Locus_'+str(i):locus.strip('\n').split('\t')[2:] for i,locus in enumerate(input)} #loci dictionary with loci number as key and genes as values

    fa_set = set(val for values in loci_diction.values() for val in values) #fa disease associated genes from the loci dictionary

    fa_diction = {}
    fa_string_set = set()
    filtered_loci_diction = {}
    with open(string_path,'r') as string_file:
        for edge in string_file:
            genes=sorted(edge.strip().split('\t'),reverse=True) #remove duplicates
            if ((genes[0] in fa_set) and (genes[1] in fa_set)):
                fa_string_set.update(genes[:2]) # fa_string_set has 330 FA genes that are in the string file
                if genes[0] in fa_diction.keys():
                    if((genes[1] not in fa_diction[genes[0]].keys()) and (genes[2] not in fa_diction[genes[0]].keys())):
                        fa_diction[genes[0]].update({genes[1]:genes[2]})
                elif genes[0] not in fa_diction.keys():
                    fa_diction[genes[0]] = {genes[1]:genes[2]} # fa_diction is a subset of string_diction with FA genes connecting to FA genes (has 1024 connections)
     
    filtered_loci_diction = {k:[element for element in v if element in fa_string_set] for k,v in loci_diction.items()} # Filtering the loci dictionary to remove FA genes with 0 interactions resulting in 330 genes

    return loci_diction,filtered_loci_diction,fa_set,fa_string_set, fa_diction

# Number of edges for each gene in the string file
def gene_frequency(string_diction,string_gene_set):
    gene_freq = {gene:0 for gene in string_gene_set}

    for outer_key,inner_diction in string_diction.items():
        for inner_key, inner_value in inner_diction.items():
            if outer_key in string_gene_set:
                gene_freq[outer_key]+=1
            if inner_key in string_gene_set:
                gene_freq[inner_key]+=1 # gene_freq is a dictionary with gene name as the key and number of edges as the value
    return gene_freq

# segregate the genes into bins
# creating a nested dictionary with bin range as the outer key and fa_genes/non-fa genes as the inner key and list of genes as the inner value
def bin_dictionary(gene_freq,num_bins=500,filter_bin_data=True):
    min_freq = min(gene_freq.values())
    max_freq = max(gene_freq.values())
    bin_width = (max_freq - min_freq)//num_bins

    bin_diction = {}
    for i in range(num_bins):
        lower_bound = int(min_freq + i*bin_width) #lower  bound of each bin
        upper_bound = int(lower_bound + bin_width) #upper bound of each bin
        bin_diction[str(lower_bound)+'-'+str(upper_bound)]={}
        bin_diction[str(lower_bound)+'-'+str(upper_bound)]["fa_genes"]=[]
        bin_diction[str(lower_bound)+'-'+str(upper_bound)]["non_fa_genes"]=[]
        for gene,freq in gene_freq.items():
            if freq >= lower_bound and freq < upper_bound: 
                if gene in fa_set:
                    bin_diction[str(lower_bound)+'-'+str(upper_bound)]["fa_genes"].append(gene)  #assigning FA genes to bins
                else:
                    bin_diction[str(lower_bound)+'-'+str(upper_bound)]["non_fa_genes"].append(gene) # assigning non-FA genes to bins

    # filter the bin_diction to remove non-fa genes that do not have any fa genes in a specific bin
    if filter_bin_data == True:
        empty_keys = []
        filtered_bin_diction = {}
        for outer_key, inner_dict in bin_diction.items():
            for inner_key, inner_value in inner_dict.items():
                if not inner_value and inner_key == 'fa_genes':
                   empty_keys.append((outer_key))  # empty keys are the bins that do not have fa genes 
                   break

        filtered_bin_diction = {k:v for k,v in bin_diction.items() if k not in empty_keys} # remove non-fa genes that do not have fa genes in the associated bin
        return bin_diction,filtered_bin_diction,empty_keys

    else:
        return bin_diction
    
# Generate two populations of 5000 random subnetworks 

# For every FA gene, match the bin category
def match_fa_nonfa_by_bucket(nested_dict,disease_gene):
    for out_g, inner_dict in nested_dict.items():
        bin_cat = out_g 
        if isinstance(inner_dict,dict): # Recursion
            result = match_fa_nonfa_by_bucket(inner_dict,disease_gene)
            if result is not None:
                return out_g
        elif disease_gene in inner_dict:
            return bin_cat                 # returns bin category or key of the dictionary that has the fa gene being searched
        
# 5000 subnetworks for FA genes
def subnetworks(filtered_loci_dictionary,num_iterations=5000,print_subnetwork=False):
    seed1 = 123  #Random seed
    seed2 = 345  #Random seed
    subnetwork_diction = {}
    null_subnetwork_diction = {}
    fa_subnetwork={}
    null_subnetwork={}
    for i in range(num_iterations):
        iteration_fa_genes=set()
        iteration_nonfa_genes=set()
        for key,value in filtered_loci_dictionary.items():
            random_index=(seed1*seed2)%len(value) #random index for FA gene will be between 0 to len(value)-1
            iteration_fa_genes.add(value[random_index]) # iteration_fa_genes is a set of 12 FA genes for a specific subnetwork
        
            # Match the bin which has that specific FA gene
            bin=match_fa_nonfa_by_bucket(bin_diction,value[random_index])
            non_fa_values = bin_diction[bin]['non_fa_genes']
            random_index_nonfa=(seed1*seed2)%len(non_fa_values) #random index for non FA gene will be between 0 and len(non_fa_values)-1 
            iteration_nonfa_genes.add(non_fa_values[random_index_nonfa]) # iteration_nonfa_genes is a set of 12 genes for a specific subnetwork
            seed1+=1
            seed2+=1
        
        # subnetworks for iteration_fa_genes list
        temp_list = []
        for outer_g, inner_d in fa_dictionary.items():
            for inner_g, edge_value in inner_d.items():
                if outer_g in iteration_fa_genes and inner_g in iteration_fa_genes:
                    temp_list.append([outer_g,inner_g,edge_value])
            subnetwork_diction['subnetwork_'+str(i+1)]=len(temp_list)
            if print_subnetwork:
                fa_subnetwork['subnetwork_'+str(i+1)]=temp_list
        
        # subnetworks for iteration_nonfa_genes list
        temp_list_2 = []
        for outer_g_2, inner_d_2 in string_dictionary.items():
            for inner_g_2, edge_value_2 in inner_d_2.items():
                if outer_g_2 in iteration_nonfa_genes and inner_g_2 in iteration_nonfa_genes:
                    temp_list_2.append([outer_g_2,inner_g_2,edge_value_2])
            null_subnetwork_diction['null_subnetwork_'+str(i+1)]=len(temp_list_2)
            if print_subnetwork:
                null_subnetwork['null_subnetwork_'+str(i+1)]=temp_list_2

    if print_subnetwork:
        return subnetwork_diction,fa_subnetwork, null_subnetwork_diction, null_subnetwork
    else:      
        return subnetwork_diction, null_subnetwork_diction
    

def permutation_test(subnetwork_diction,null_subnetwork_diction,num_permutations=10000):
    
    # observed mean difference
    # Calculate mean difference of edges in 5000 FA subnetworks
    mean_fa_sub = sum(subnetwork_diction.values())/len(subnetwork_diction.values())
    mean_nonfa_sub = sum(null_subnetwork_diction.values())/len(null_subnetwork_diction.values())
    mean_diff_observed = abs(mean_nonfa_sub - mean_fa_sub)
    all_networks = {**subnetwork_diction,**null_subnetwork_diction}
    extreme=0
    perm_stat_list = []

    seed1 = 123
    seed2 = 345
    
    # calculate permutation statistic (mean)
    for i in range(num_permutations):
        random_index_set = set()
        all_network_keys = list(subnetwork_diction.keys())+list(null_subnetwork_diction.keys()) #merge labels
        swap_labels = [net for net in all_network_keys] 
        for index in range(len(all_network_keys)):
            random_index=(seed1*seed2)%len(all_network_keys) #random index for swapping
            seed1+=1
            seed2+=1
            random_index_set.add(random_index)

            # swapping labels without replacement 
            # swap each index with a random_index
            swap_labels[index], swap_labels[random_index] = swap_labels[random_index], swap_labels[index]
        
        seed1+=1
        seed2+=1

        #mean difference/permutation statistic
        perm_mean1 = sum(all_networks[subnet] for subnet in swap_labels[:len(subnetwork_diction.keys())])/len(subnetwork_diction.keys())
        perm_mean2 = sum(all_networks[subnet] for subnet in swap_labels[len(subnetwork_diction.keys()):])/len(subnetwork_diction.keys())
        perm_stat = abs(perm_mean1-perm_mean2)
        perm_stat_list.append(perm_stat) #append permutation statistics
    
        # Number of times the permutation statistic is as extreme and more extreme than the observed statistic
        if perm_stat >= mean_diff_observed:
           extreme+=1

    # empirical p value is the probability that the permutation statistic is as extreme or more extreme than the observed statistic
    empirical_p_value = extreme/num_permutations

    return perm_stat_list, extreme, empirical_p_value

# parsing command line arguments
parser=argparse.ArgumentParser()
parser.add_argument('-dl','--diseaseloci',type=str,help='Disease loci file with list of genes')
parser.add_argument('-n','--network',type=str,help='Network file in tab delimited format')
parser.add_argument('--nBins',type=int,help='Number of bins to divide the data into based on the number of edges')
parser.add_argument('--nSub',type=int,help='Number of random subnetworks to generate')
parser.add_argument('--nPerm',type=int,help='Number of iterations for the permutation test')

args = parser.parse_args()

# calling the functions
string_gene_set,string_dictionary = string_to_diction(args.network)
loci_dictionary,filtered_loci_dictionary,fa_set,fa_string_set,fa_dictionary = input_to_diction(args.diseaseloci,args.network)
gene_freq = gene_frequency(string_dictionary,string_gene_set)
bin_diction,filtered_bin_diction,empty_keys=bin_dictionary(gene_freq,num_bins=args.nBins)
subnetwork_diction, fa_subnetwork, null_subnetwork_diction,null_subnetwork = subnetworks(filtered_loci_dictionary,num_iterations=args.nSub,print_subnetwork=True)
perm_stat_list,extreme,empirical_p_value = permutation_test(subnetwork_diction,null_subnetwork_diction,num_permutations=args.nPerm)

mean_fa_sub = sum(subnetwork_diction.values())/len(subnetwork_diction.values())
mean_nonfa_sub = sum(null_subnetwork_diction.values())/len(null_subnetwork_diction.values())
mean_diff_observed = abs(mean_nonfa_sub - mean_fa_sub)

# Distribution of gene frequencies
plt.hist(list(gene_freq.values()), bins=100, color='steelblue', edgecolor='black', linewidth=1.2)
# Add a grid to the plot
plt.grid(axis='y', alpha=0.75)
# Add labels and title
plt.xlabel('Number of edges')
plt.ylabel('Frequency')
plt.title('Number of edges in each gene')
plt.close()

# Distribution of 5000 FA subnetwork edge counts
plt.hist(list(subnetwork_diction.values()),color='steelblue', edgecolor='black', linewidth=1.2)
# Add a grid to the plot
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Edges in each subnetwork')
plt.ylabel('Frequency')
plt.title('Distribution of FA gene subnetworks')
plt.savefig('FA_subnetwork_distribution.png')
plt.close()

# Distribution of 5000 non-FA subnetwork edge counts
plt.hist(list(null_subnetwork_diction.values()), color='steelblue', edgecolor='black', linewidth=1.2)
# Add a grid to the plot
plt.grid(axis='y', alpha=0.75)
# Add labels and title
plt.xlabel('Edges in each subnetwork')
plt.ylabel('Frequency')
plt.title('Distribution of FA gene subnetworks (Null case)')
plt.savefig('Null_subnetwork_Distribution.png')
plt.close()

# Distribution of 10000 permutation statistic 
plt.hist(perm_stat_list, color='steelblue', edgecolor='black', linewidth=1.2)
plt.axvline(x = mean_diff_observed, color = 'r', linestyle = 'dashed')
# Add a grid to the plot
plt.grid(axis='y', alpha=0.75)
# Add labels and title
plt.xlabel('Permutation statistic')
plt.ylabel('Frequency')
plt.title('permutation statistic vs observed statistic')
plt.savefig('Permutation_statistic.png')
plt.close()

# save the 5000 random FA subnetworks
with open('FA_subnetwork.txt', 'w') as file:
    for k, v in fa_subnetwork.items():
        for sub in v:
            file.write("{}\t{}\n".format(k,'\t'.join(sub)))

# save the 5000 random non-FA subnetworks
with open('Null_subnetworks.txt', 'w') as file:
    for k, v in null_subnetwork.items():
        for sub in v:
            file.write("{}\t{}\n".format(k,'\t'.join(sub)))

print("Empirical p value is ",empirical_p_value)

# run time
end = (time.time())
print(end-start)