#--------------------------------------------------------------------------#
#- Subnetwork 1                                                           -#
#- FA gene nodes connected on different loci                              -#
#- Does not take the FA genes connected within the same loci into account -#
#--------------------------------------------------------------------------#
#- Subnetwork 2                                                           -#
#- FA gene nodes connected on different and same loci                     -#
#- Takes into account the FA genes connected within the same loci         -#
#--------------------------------------------------------------------------#
#- Subnetwork 3                                                           -#
#- FA gene node interactions regardless of loci.                          -#
#- Non-FA gene nodes connecting two or more FA genes                      -#
#- Does not take the FA genes connected within the same loci into account -#
#--------------------------------------------------------------------------#

import time
import sys

start = time.time()

# Convert loci file into a dictionary of lists with loci number as keys and list of genes as values
def loci_dictionary(loci_path):
    loci_diction = {}
    with open(loci_path,'r') as loci:
       for i,locus in enumerate(loci):
        geneList = []
        loci_diction[i]=locus.strip('\n').split('for ')[1].split('\t')  # Data Wrangling to fetch the gene names in each loci
    return loci_diction

# Convert loci file into a list 
def loci_list(loci_path):
    with open(loci_path,'r') as input_file:
        loci=input_file.readlines()
        loci=[locus.strip().split('for ')[1].split('\t') for locus in loci]
        fa_list=[gene for genes in loci for gene in genes]
    return(fa_list)

# Convert the tab-delimited String db file into a nested dictionary
# Removes duplicate interactions (i.e., one of the two rows "geneA geneB edge" and "geneB geneA edge" is removed)
def string_dictionary(string_path):
    string_diction={}
    with open(string_path, 'r') as string_file:
       for edge in string_file:
           genes = sorted(edge.strip().split('\t'),reverse=True)  # sort column1 and column2 of each row
           if genes[0] in string_diction.keys():
              if (genes[1] not in string_diction[genes[0]].keys()) and (genes[2] not in string_diction[genes[0]].get(genes[1], {})):
                string_diction[genes[0]].update({genes[1]:genes[2]})
           elif genes[0] not in string_diction.keys():
               string_diction[genes[0]] = {genes[1]:genes[2]} 
    return(string_diction)

#--------------------------------------------------------------------------#
#- Subnetwork 1                                                           -#

# Retrieve subnetwork dictionary (Only includes between loci gene interactions)
# Loop through locus key and check if the FA gene is from the same loci
def subnetwork_1(loci_geneList,string_geneList):
    subnetwork_diction={}
    for locus_key,locus_value in loci_geneList.items():
        include_loci = [k for k in loci_geneList.keys() if k!= locus_key]  
        compare_list=[]
        for i in include_loci:
            compare_list+=loci_geneList[i]
        for node in locus_value:
            if node in string_geneList.keys():
               subnetwork_diction[node] = {}
               for comp in compare_list:
                   if comp in string_geneList[node].keys():
                      subnetwork_diction[node][comp]=string_geneList[node][comp]
    return(subnetwork_diction)

#--------------------------------------------------------------------------#
#- Subnetwork 2                                                           -#

# Retrieve subnetwork dictionary (Includes between and within loci gene interactions)
def subnetwork_2(loci_geneList,string_geneList):
    subnetwork_diction={}
    for locus_key,locus_value in loci_geneList.items():
        include_loci = [k for k in loci_geneList.keys()]
        compare_list=[]
        for i in include_loci:
            compare_list+=loci_geneList[i]
        for node in locus_value:
            if node in string_geneList.keys():
               subnetwork_diction[node] = {}
               for comp in compare_list:
                   if comp in string_geneList[node].keys():
                      subnetwork_diction[node][comp]=string_geneList[node][comp]
    return(subnetwork_diction)

#--------------------------------------------------------------------------#
#- Subnetwork 3                                                           -#

# Fetch non-fa genes that are connected to fa genes and append to a list
# count the number of fa genes that each non-fa gene is connected to 
# keep the genes that are connected to at least two fa genes
# This would be a network path connecting fa genes
# string_subset also gives subnetwork 2 (alternate code to the above function)
def non_fa_list(fa_list,string_diction):
    string_subset={}
    not_fa=[]
    for gene1, sub_string in string_diction.items():
        for gene2 in sub_string.keys():
            if gene1 in fa_list and gene2 in fa_list:
                if gene1 not in string_subset.keys():
                    string_subset[gene1] = {gene2:sub_string[gene2]} 
                else:
                    string_subset[gene1].update({gene2:sub_string[gene2]}) # also gives subnetwork2
            elif gene1 not in fa_list and gene2 in fa_list:
                if gene1 not in not_fa:
                    not_fa.append(gene1)
            elif gene1 in fa_list and gene2 not in fa_list:
                if gene1 not in not_fa:
                    not_fa.append(gene2)

    not_fa_nodes_count={}  
    for node in not_fa:
        if node in not_fa_nodes_count.keys():
            not_fa_nodes_count[node]+=1        # counting the number of fa genes a non-fa gene is connected to
        else:
            not_fa_nodes_count[node]=1

    not_fa_nodelist = [k for k,v in not_fa_nodes_count.items() if v >= 2] # keep non-fa genes that are connected to atleast 2 fa genes
    final_nodelist = fa_list+not_fa_nodelist    # append to the fa list
    return string_subset,final_nodelist

def subnetwork_3(string_diction,final_nodelist):
    string_nonfa={}
    for gene1, sub_string in string_diction.items():
        for gene2 in sub_string.keys():
            if gene1 in final_nodelist and gene2 in final_nodelist:
                if gene1 not in string_nonfa.keys():
                    string_nonfa[gene1] = {gene2:sub_string[gene2]} 
                else:
                    string_nonfa[gene1].update({gene2:sub_string[gene2]})
    return(string_nonfa)

# Dictionary to a text file
def diction_to_text(nested_dict):
    result = []
    for key1, inner_diction in nested_dict.items():
        for key2, value in inner_diction.items():
            result.append((key1, key2, value))
    return result

#loci_path = 'Input.gmt.txt'
loci_path=sys.argv[1]
loci_geneList = loci_dictionary(loci_path)
fa_list = loci_list(loci_path)

#string_path='STRING 1.txt'
string_path=sys.argv[2]
string_genediction = string_dictionary(string_path)

subnetwork_diction_between = subnetwork_1(loci_geneList,string_genediction)
sub_data_between = diction_to_text(subnetwork_diction_between)

subnetwork_diction = subnetwork_2(loci_geneList,string_genediction)
sub_data = diction_to_text(subnetwork_diction)

string_subset,final_nodelist=non_fa_list(fa_list,string_genediction)
string_nonfa = subnetwork_3(string_genediction,final_nodelist)
subnetwork3 = diction_to_text(string_nonfa)

with open('subnetwork1.txt', 'w') as file:
    for row in sub_data_between:
        file.write('\t'.join(row) + '\n')

with open('subnetwork2.txt', 'w') as file:
    for row in sub_data:
        file.write('\t'.join(row) + '\n')

with open('subnetwork3.txt', 'w') as file:
    for row in subnetwork3:
        file.write('\t'.join(row) + '\n')

# run time
end = (time.time())
print(end-start)