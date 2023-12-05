#------------------------------------------------------------------------------
# Section 1: Data Wrangling; Functions: data_to_diction; string_to_diction
# Section 2a: Disease associated gene subnetwork; Functions: first_parent_population
# Section 2b: GA optimization (Mutation); Functions: parent_mean_edgeweight, generate_random_excluding, ga_mutation,
# Section 2c: GA optimization (Mating); Functions: ga_mating
# Section 3a: Null case

import random
#-------------------------------------------------------------------------------------------------------------------#
# Function to convert input files into objects for the later steps                                                  #
# loci_diction: dictionary of loci (Key) and genes (value)                                                          #
# fa_set: set with all the fa genes from loci_diction                                                               #
# fa_diction: Subset of the STRING.txt file but with FA gene connected to another FA gene                           #
# Inputs: input path for Input.gmt.txt and string path for STRING.txt                                               #
# Outputs: loci_diction (key as loci ID and values as genes),                                                       #
# Outputs: fa_set (set of 584 genes that are in the Input file),                                                    #
# Outputs: fa_diction (dictionary of FA gene to FA gene connections total 1024),                                    #
# Outputs: fa_string_set (set of FA genes that are in the above 1024 list. Total 329 genes),                        #
# Outputs: string_set (set of FA genes that exist in STRING file, they may be connected to non-FA genes. Total 508) #
#-------------------------------------------------------------------------------------------------------------------#
## Function converts input files into dictionaries
def data_to_diction(input_path,string_path):
    with open(input_path, 'r') as input:
        # a) loci 1 to 12 in dictionary
        loci_diction = {'Locus_'+str(i):locus.strip('\n').split('\t')[2:] for i,locus in enumerate(input)}

        # b) set with all fa genes
        fa_set = set(val for values in loci_diction.values() for val in values)

        fa_diction = {}
        fa_string_set = set()
        string_set = []
        all_genes_set = set()
        
        with open(string_path,'r') as string_file:
            for edge in string_file:
                genes=sorted(edge.strip().split('\t'),reverse=True)
                if((genes[0] in fa_set) and (genes[1] in fa_set)):
                    # c) set with fa genes that are connected to other fa genes in string file
                    fa_string_set.update(genes[:2])
                    if genes[0] in fa_diction.keys():
                        if((genes[1] not in fa_diction[genes[0]].keys()) and (genes[2] not in fa_diction[genes[0]].keys())):
                            # d) dictionary with fa gene to fa gene connections
                            fa_diction[genes[0]].update({genes[1]:genes[2]})
                    elif genes[0] not in fa_diction.keys():
                        fa_diction[genes[0]] = {genes[1]:genes[2]}

                elif((genes[0] in fa_set) or (genes[1] in fa_set)):
                    if(genes[0] in fa_set):
                        # e) all fa gene connections from string file (i.e., fa to non-fa and fa to fa)
                        string_set.append(genes[0])
                    if(genes[1] in fa_set):
                        string_set.append(genes[1]) 
                        
                # f) all genes from the string file and the loci file
                all_genes_set.add(genes[0])
                all_genes_set.add(genes[1])

        all_genes_set|=fa_set
                
    string_set = set(string_set)

    return loci_diction,fa_set,fa_diction,fa_string_set,string_set,all_genes_set


# Function to convert the string file to a nested dictionary and set of all the genes from the string file
def string_to_diction(string_path):
    string_diction={}
    string_gene_set = set()
    with open(string_path,'r') as string_file:
        for edge in string_file:
            genes = sorted(edge.strip().split('\t'),reverse=True)  #Remove duplicates
            # a) all genes from the string file (does not have the fa genes without any connections)
            string_gene_set.update(genes[:2])
            # b) string file in the form of a dictionary
            if genes[0] in string_diction:
                if ((genes[1] not in string_diction[genes[0]].keys()) and genes[1] not in string_diction[genes[0]].keys()):
                    string_diction[genes[0]][genes[1]]=genes[2]
            elif genes[0] not in string_diction:
                string_diction[genes[0]] = {genes[1]:genes[2]}
    return string_gene_set,string_diction    #string_gene_set has all the unique gene names from the string_diction


## Step1: Generate 5000 random subnetworks i.e., the parent population
def first_parent_population(loci_diction,fa_diction,num_iterations=5000,min_edges=0):
    
    subnetwork_diction = {}
    subnetwork_diction_edge_counts={}
    first_parent_subnetwork={}
    random_numbers={}
    #iteration_fa_genes_diction = {}
    
    i=0
    random.seed(123)

    while len(subnetwork_diction.keys()) < num_iterations:

        iteration_fa_genes=[]
        random_index_set=[]
    
        for key,value in loci_diction.items():
            random_index=random.randint(0,len(value)-1)
            iteration_fa_genes.append(value[random_index])
            random_index_set.append(random_index)

        # if the subnetwork already doesn't exist then add it to the dictionary 
        if not any(iteration_fa_genes == values for values in first_parent_subnetwork.values()):
            first_parent_subnetwork['subnetwork_'+str(i+1)] = iteration_fa_genes 
    
        temp_list = []
        for outer_g, inner_d in fa_diction.items():
            for inner_g, edge_value in inner_d.items():  
                if outer_g in first_parent_subnetwork['subnetwork_'+str(i+1)] and inner_g in first_parent_subnetwork['subnetwork_'+str(i+1)]:
                    temp_list.append([outer_g,inner_g,edge_value])

                if (len(temp_list) >= min_edges):
                    subnetwork_diction_edge_counts['subnetwork_'+str(i+1)]=len(temp_list)
                    subnetwork_diction['subnetwork_'+str(i+1)]=temp_list
        
        i+=1
        
    return first_parent_subnetwork,subnetwork_diction,subnetwork_diction_edge_counts

def parent_mean_edgeweight(subnetwork_diction,num_iterations=5000):
    #parent_subnet_edgeweight={}
    mean_parent_edgeweight=0
    for subnet_index, edge_list in subnetwork_diction.items():
        #parent_subnet_edgeweight[subnet_index]=[]
        for edge_num,edges in enumerate(edge_list):
            #parent_subnet_edgeweight[subnet_index].append(edges[2])
            mean_parent_edgeweight+=float(edges[2])

    return float(mean_parent_edgeweight/num_iterations)      
#parent_subnet_edgeweight    

def generate_random_excluding(range_start, range_end, exclude_number):
    valid_numbers = [num for num in range(range_start, range_end + 1) if num != exclude_number]
    return random.choice(valid_numbers)

## Step 2: Mutation with 5% probability
## Given a parent population, the following function returns a mutated population and the normalized weights of this subnetwork

def ga_mutation(loci_diction,fa_diction,parent_subnetwork,num_iterations=5000):

    subnetwork_diction_updated={}
    #subnetwork_mutated_weighted={}
    subnetwork_normalized_weighted={}
    
    random.seed(1234)

    for subnetwork_index,subnet_genes_list in parent_subnetwork.items():
        for locus,mut_genes in enumerate(subnet_genes_list):
        
            prob_5 = random.randint(1,20)
            if  prob_5 == 1 :
                mutate_gene=loci_diction['Locus_'+str(locus)].index(mut_genes)
                rand_replace=generate_random_excluding(0,len(loci_diction['Locus_'+str(locus)])-1,mutate_gene)
                parent_subnetwork[subnetwork_index][locus]=loci_diction['Locus_'+str(locus)][rand_replace]

    for subnetwork_index,genes_list in parent_subnetwork.items():
        temp_list_2 = []
        
        for outer_gene,inner_diction in fa_diction.items():
            if outer_gene in genes_list: 
                for inner_gene,edge_value in inner_diction.items():
                    #if all(current_fa_gene in genes_list for current_fa_gene in (outer_gene,inner_gene)):
                    if inner_gene in genes_list:
                        temp_list_2.append([outer_gene,inner_gene,edge_value])

        subnetwork_diction_updated[subnetwork_index]=temp_list_2
        #subnetwork_mutated_weighted[subnetwork_index]=(sum([float(edge_weights[2]) for edge_weights in temp_list_2]))
        subnetwork_normalized_weighted[subnetwork_index]=(sum([float(edge_weights[2]) for edge_weights in temp_list_2]))**3
        
    sum_norm_weights_parent=sum(subnetwork_normalized_weighted.values())
    
    normalized_weights_parent={subnet:edge_weighted/sum_norm_weights_parent for subnet,edge_weighted in subnetwork_normalized_weighted.items()}

    return normalized_weights_parent,subnetwork_diction_updated,parent_subnetwork,subnetwork_normalized_weighted#,subnetwork_mutated_weighted,


## Step 3: Mating 
## Given a mutated parent population of subnetworks, this function enriches for dense subnetworks
def ga_mating(fa_diction,normalized_weights_parent,parent_mean_edgeweights,parent_subnetwork,num_iterations=5000,min_edges=0):
    random.seed(123)

    mate_population = {}
    mate_weighted={}
    iteration=0
    mate_parent_subnetwork={}
    
    while len(mate_population.keys()) < num_iterations:
        subnet_pair=random.choices(list(normalized_weights_parent.keys()),
                              weights=list(normalized_weights_parent.values()),
                              k=2)
        sub_parent={}
        mate_parent_subnetwork['subnetwork_'+str(iteration+1)]=[]
    
        for r in range(12):
            pair_random=random.randint(1,2)
            mate_parent_subnetwork['subnetwork_'+str(iteration+1)].append(parent_subnetwork[subnet_pair[pair_random-1]][r])

        temp_list_3=[]
        for outer_gene,inner_diction in fa_diction.items():
            if outer_gene in mate_parent_subnetwork['subnetwork_'+str(iteration+1)]:
                for inner_gene,edge_value in inner_diction.items(): 
                    if inner_gene in mate_parent_subnetwork['subnetwork_'+str(iteration+1)]:
                        temp_list_3.append([outer_gene,inner_gene,edge_value])

        #if (len(temp_list_3) >= min_edges):
        mate_population['subnetwork_'+str(iteration+1)]=temp_list_3
        mate_weighted['subnetwork_'+str(iteration+1)]=(sum([float(edge_weight[2]) for edge_weight in temp_list_3]))
            #mate_normalized_weighted['subnetwork_'+str(iteration+1)]=(sum([float(edge_weights[2]) for edge_weights in temp_list_3]))**3
            
        iteration+=1 
    
    tot_mate_weighted=sum(mate_weighted.values())
    #mate_normalized_weighted = {key: value / tot_mate_normalized_weighted  for key, value in mate_normalized_weighted.items()}

    mean_current_weights=(sum(mate_weighted.values()))/num_iterations
    mean_parent_weights=parent_mean_edgeweights
    mean_delta = mean_current_weights - mean_parent_weights
    mean_delta_perc = (mean_delta*100)/mean_parent_weights
    #mean_delta_perc = (mean_delta)/(sum(list(normalized_weights_parent.values()))/num_iterations)

    return mate_population,mate_parent_subnetwork,mate_weighted,mean_delta,mean_delta_perc

# Number of edges for each gene in the string file
def gene_frequency(string_diction,all_genes_set):
    gene_freq = {gene:0 for gene in all_genes_set}

    for outer_key,inner_diction in string_diction.items():
        for inner_key, inner_value in inner_diction.items():
            if outer_key in all_genes_set:
                gene_freq[outer_key]+=1
            if inner_key in all_genes_set:
                gene_freq[inner_key]+=1 # gene_freq is a dictionary with gene name as the key and number of edges as the value
    return gene_freq

# segregate the genes into bins
# creating a nested dictionary with bin range as the outer key and fa_genes/non-fa genes as the inner key and list of genes as the inner value
def bin_dictionary(gene_freq,fa_set,num_bins=100,filter_bin_data=True):
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

def generate_nonfa_from_fa(optimized_FA_subnetwork,random_seed,bin_diction,string_diction):
    nonFA_parent={}
    nonFA_parent_size={}
    iteration=0
    random.seed(random_seed)
    all_nonFA_null=set()
    nonFA_parent_genes={}
    
    for subnet_index,list_genes in optimized_FA_subnetwork.items():
        iteration_nonfa_genes=[]
        nonFA_parent_genes['subnetwork_'+str(iteration+1)]=[]
        for genes in list_genes:
            bin=match_fa_nonfa_by_bucket(bin_diction,genes)
            non_fa_values = bin_diction[bin]['non_fa_genes']

            random_index_nonfa=random.randint(0,len(non_fa_values)-1)#random index for non FA gene will be between 0 and len(non_fa_values)-1 
            iteration_nonfa_genes.append(non_fa_values[random_index_nonfa])
            all_nonFA_null.add(non_fa_values[random_index_nonfa])
        nonFA_parent_genes['subnetwork_'+str(iteration+1)]=iteration_nonfa_genes

        temp_list_4=[]

        for outer_gene, inner_diction in string_diction.items():
            if outer_gene in iteration_nonfa_genes:
                for inner_gene,inner_edge in inner_diction.items():
                    if inner_gene in iteration_nonfa_genes:
                        temp_list_4.append([outer_gene, inner_gene, inner_edge])
    
        iteration+=1
        nonFA_parent['subnetwork_'+str(iteration+1)]=temp_list_4
        nonFA_parent_size['subnetwork_'+str(iteration+1)]=len(temp_list_4)

    return nonFA_parent,nonFA_parent_size,nonFA_parent_genes

## Create random loci for non-FA genes
def generate_random_loci_NonFA(nonFA_parent_genes,all_genes_set):

    non_fa_loci = {'Locus_{}'.format(l): [] for l in range(12)}

    all_nonfa_subnet_genes=set()

    for nonfa_sub,non_fa in nonFA_parent_genes.items():
        for nonfa_locus,non_fa_g in enumerate(non_fa):
            all_nonfa_subnet_genes.add(non_fa_g)
            nonfa_index=non_fa.index(non_fa_g)
            non_fa_loci['Locus_{}'.format(nonfa_index)].append(non_fa_g)
            
    for i in all_genes_set:
        if i not in all_nonfa_subnet_genes:
            rand_nonfa_loci=random.randint(0,11)
            non_fa_loci['Locus_'+str(rand_nonfa_loci)].append(i)

    return non_fa_loci


def gene_score(iteration_fa_genes_diction_2,fa_diction,loci_diction):
    gene_score = {}
    for key,value in iteration_fa_genes_diction_2.items():
    
        chosen_gene = iteration_fa_genes_diction_2[key]
        temp_list_4=[]
        for outer_g, inner_d in fa_diction.items():
            for inner_g, edge_value in inner_d.items():
                if outer_g != chosen_gene and inner_g != chosen_gene:
                    if outer_g in iteration_fa_genes_diction_2.values() and inner_g in iteration_fa_genes_diction_2.values():
                        temp_list_4.append([outer_g,inner_g,edge_value])
                        
        for loci_fa_gene in loci_diction[key]:
            iteration_fa_genes_diction_2[key] = loci_fa_gene
            temp_list_5=[]
            gene_score_diction = {}

            for outer_g, inner_d in fa_diction.items():
                for inner_g, edge_value in inner_d.items():
                    if outer_g in iteration_fa_genes_diction_2.values() and inner_g in iteration_fa_genes_diction_2.values():
                        temp_list_5.append([outer_g,inner_g,edge_value])

            if key not in gene_score.keys():
                gene_score[key]={loci_fa_gene:(len(temp_list_5)-len(temp_list_4))}
            else:
                gene_score[key].update({loci_fa_gene:(len(temp_list_5)-len(temp_list_4))})
 
    return(gene_score)

