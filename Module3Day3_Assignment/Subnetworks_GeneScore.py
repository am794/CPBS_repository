#---------------------------------------------------------------------------------#
#- Generate 5000 random subnetworks with one representative gene from each loci  -#
#- Calculate average gene scores for candidate and representative genes          -#
#- Visualizing the subnetworks                                                   -#
#---------------------------------------------------------------------------------#             

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

# Run time start
start = time.time()

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
def data_to_diction(input_path,string_path):
    with open(input_path, 'r') as input:
        loci_diction = {'Locus_'+str(i):locus.strip('\n').split('\t')[2:] for i,locus in enumerate(input)}

        fa_set = set(val for values in loci_diction.values() for val in values)

        fa_diction = {}
        fa_string_set = set()
        string_set = []

        with open(string_path,'r') as string_file:
            for edge in string_file:
                genes=sorted(edge.strip().split('\t'),reverse=True)
                if((genes[0] in fa_set) and (genes[1] in fa_set)):
                    fa_string_set.update(genes[:2])
                    if genes[0] in fa_diction.keys():
                        if((genes[1] not in fa_diction[genes[0]].keys()) and (genes[2] not in fa_diction[genes[0]].keys())):
                            fa_diction[genes[0]].update({genes[1]:genes[2]})
                    elif genes[0] not in fa_diction.keys():
                        fa_diction[genes[0]] = {genes[1]:genes[2]}

                elif((genes[0] in fa_set) or (genes[1] in fa_set)):
                    if(genes[0] in fa_set):
                        string_set.append(genes[0])
                    if(genes[1] in fa_set):
                        string_set.append(genes[1])

    string_set = set(string_set)

    return loci_diction,fa_set,fa_diction,fa_string_set,string_set


#----------------------------------------------------------------#
# Calculating gene scores                                        #
# Generate empty case score                                      #
# Generate candidate gene score                                  #
# Subtract empty case score from candidate score                 #
# Returns the sum of these scores for each subnetwork            #
# Input parameters: All representative genes for each iteration  #
# Output: Raw candidate gene scores                              #
#----------------------------------------------------------------#
def gene_score(iteration_fa_genes_diction_2):
    gene_score = {}

    for key,value in iteration_fa_genes_diction_2.items():

        # Empty case
        chosen_gene = iteration_fa_genes_diction_2[key]
        temp_list_2=[]
        for outer_g, inner_d in fa_diction.items():
            for inner_g, edge_value in inner_d.items():
                if outer_g != chosen_gene and inner_g != chosen_gene:
                    if outer_g in iteration_fa_genes_diction_2.values() and inner_g in iteration_fa_genes_diction_2.values():
                        temp_list_2.append([outer_g,inner_g,edge_value])

        #raw candidate gene scores per iteration
        for loci_fa_gene in loci_diction[key]:
            iteration_fa_genes_diction_2[key] = loci_fa_gene
            temp_list=[]
            gene_score_diction = {}

            for outer_g, inner_d in fa_diction.items():
                for inner_g, edge_value in inner_d.items():
                    if outer_g in iteration_fa_genes_diction_2.values() and inner_g in iteration_fa_genes_diction_2.values():
                        temp_list.append([outer_g,inner_g,edge_value])

            # In each iteration, subtract the empty case score from the candidate scores
            if key not in gene_score.keys():
                gene_score[key]={loci_fa_gene:(len(temp_list)-len(temp_list_2))}
            else:
                gene_score[key].update({loci_fa_gene:(len(temp_list)-len(temp_list_2))})
 
    return(gene_score)

#------------------------------------------------------------------------------------------------------------#
# Calculate gene averge, given the sum of all scores                                                         #
# If the FA gene exists in the STRING file, then compute the average                                         #
# Else give an NA to that specific gene                                                                      #
# Input parameters: loci_diction, gene_scores_loci (sum of raw candidate scores),                            #
# Input parameters: compare_set (if gene doesn't exist in this set, give NA) Use string_set or fa_string_set #
# Input parameters: num_subnet (number of subnetworks, default is 5000)                                      #
# Output: Final average gene scores                                                                          #
#------------------------------------------------------------------------------------------------------------#
def gene_scores_avg(loci_diction,gene_scores_loci,compare_set,num_subnet=5000):
    gene_scores_final = {loci: {inner_key: 0 for inner_key in loci_genes} for loci, loci_genes in loci_diction.items()}

    for locus,inner_d in gene_scores_loci.items():
        for inner_g,score_sum in inner_d.items():
            if inner_g not in compare_set:
                gene_scores_final[locus][inner_g]='NA'
            else:
                gene_scores_final[locus][inner_g]="{:.4f}".format(float(score_sum)/num_subnet) #averge 

    return gene_scores_final

#---------------------------------------------------------------------------------------------------#
# Generate 5000 random subnetworks with one representative gene from each loci                      #
# Create a dictionary of all the genes from each iteration                                          #
# Create a dictionary of all the subnetworks                                                        #
# Call gene_score and compute the sum of gene scores for candidate genes                            #
# Input Parameters: loci_diction, num_subnet (default = 5000)                                       #
# Input parameters: min_edges (minimum number of edges that a subnetwork should have. Default is 0) #
# Ouputs: fa_subnetwork (all subnetworks),                                                          #
# Outputs: gene_scores_loci (sum of candidate gene scores across the iterations.                    #
# Input this to gene_scores_avg functio)                                                            #
# Outputs: iteration_fa_genes_diction_2 (dictionary with candidate genes across iterations)         #
#---------------------------------------------------------------------------------------------------#
def subnetwork_genescore(loci_diction,num_subnet=5000,min_edges=0):
    
    subnetwork_diction = {}
    fa_subnetwork={}
    iteration_fa_genes_diction = {}
    iteration_fa_genes_diction_2 = {}
    gene_score_func = {}
    #fa_subnetwork_list=[]
    gene_scores_loci = {loci: {inner_key: 0 for inner_key in loci_genes} for loci, loci_genes in loci_diction.items()}
    random.seed(123)
    i=0

    # Create 5000 subnetworks
    while len(fa_subnetwork.keys()) < num_subnet:
         iteration_fa_genes=set()
    
         for key,value in loci_diction.items():  #iterate through each loci and genes
             random_index = random.randint(0, len(value)-1) # generate a random number between 0 and the last index of that specific loci
             iteration_fa_genes.add(value[random_index]) #Make a list of representative genes for each loci

         # if the subnetwork already doesn't exist then add it to the dictionary 
         if not any(iteration_fa_genes == values for values in iteration_fa_genes_diction.values()):
            iteration_fa_genes_diction['iteration_'+str(i+1)] = iteration_fa_genes      

         temp_list = []
         for outer_g, inner_d in fa_diction.items():
             for inner_g, edge_value in inner_d.items():  
                 if outer_g in iteration_fa_genes_diction['iteration_'+str(i+1)] and inner_g in iteration_fa_genes_diction['iteration_'+str(i+1)]:
                    temp_list.append([outer_g,inner_g,edge_value])
                    #fa_subnetwork_list.append([outer_g,inner_g,edge_value])
              
                 if (len(temp_list) > min_edges):
                    subnetwork_diction['subnetwork_'+str(i+1)]=len(temp_list)
                    fa_subnetwork['subnetwork_'+str(i+1)]=temp_list
                    iteration_fa_genes_diction_2={'Locus_'+str(ind):loc_gene for ind,loc_gene in enumerate(iteration_fa_genes)}  

         gene_score_func = gene_score(iteration_fa_genes_diction_2)
         for loci_num,inner_diction in gene_score_func.items():
             for loci_gene,score in inner_diction.items():
                 gene_scores_loci[loci_num][loci_gene]=gene_scores_loci[loci_num].get(loci_gene,0)+score
         i+=1

    return fa_subnetwork,gene_scores_loci,iteration_fa_genes_diction_2


#------------------------------------------------------------------------------#
# Visualize First n subnetworks                                                #
# Input parameters: gene_scores_final, fa_subnetwork, loci_diction             #
# Input params: file_path is the path for saving the network                   #
# Input params: num_sub_viz (First N subnetworks to visualize)                 #
#------------------------------------------------------------------------------#
def fa_subnet_viz(gene_scores_final,fa_subnetwork,loci_diction,file_path="./network.png",num_subs_viz=5000):
    nodes_subs=set()
    sub_fa_network = []
    for i, fa in enumerate(fa_subnetwork.values()):
        if i==num_subs_viz:
            break
        else:
            for fa_subs in fa:
                nodes_subs.add(fa_subs[0])
                nodes_subs.add(fa_subs[1])
                sub_fa_network.append(fa_subs)

    gene_scores_nodes={}

    for l,g in loci_diction.items():
        for gns in g:
            if gns in nodes_subs:
                gene_scores_nodes[gns]=l

    node_size={}
    for inner_scores in gene_scores_final.values():
        for g,s in inner_scores.items():
            if g in nodes_subs:
                node_size[g] = math.log10(float(s)+0.01)


    #G.clear()
    G=nx.Graph()

    for nd in nodes_subs:
        G.add_node(nd)

    for edge_ in sub_fa_network:
        j = ([edge_elem.strip().split(',') for edge_elem in edge_])
        k=j[0]+j[1]
        G.add_edge(*k)

    unique_nodes = set(gene_scores_nodes.values())
    random_colors = {node: (random.random(), random.random(), random.random()) for node in unique_nodes}
    node_colors = {node: random_colors[loci_id] 
                  for node, loci_id in gene_scores_nodes.items()}

    fig, ax = plt.subplots(figsize=(15, 15))
    pos = nx.spring_layout(G)
    Graph(G,
          node_color=node_colors, # indicates the community each belongs to  
          node_edge_width=0.2,     # no black border around nodes 
          node_size=node_size,
          edge_width=0.1,        # use thin edges, as they carry no information in this visualisation
          edge_alpha=0.6,        # low edge alpha values accentuates bundles as they appear darker than single edges
          node_alpha=0.6,
          node_layout='community', 
          edge_layout='bundled',
          node_layout_kwargs=dict(node_to_community=gene_scores_nodes),
          ax=ax,
          pos=pos,
          with_labels=False
    )

    legend_elements = [Patch(color=color, label=loci_id) for loci_id, color in random_colors.items()]

    ax.legend(handles=legend_elements, scatterpoints=1, frameon=False, labelspacing=1, title='LOCUS ID')

    fig.savefig(file_path)
    plt.close(fig)

    return file_path


# parsing command line arguments
parser=argparse.ArgumentParser()
parser.add_argument('-dl','--diseaseloci',type=str,help='Disease loci file with list of genes')
parser.add_argument('-n','--network',type=str,help='Network file in tab delimited format')
parser.add_argument('--nSub',type=int,help='Number of random subnetworks to generate. Default is 5000')
parser.add_argument('--minEdges',type=int,help='Minimum number of edges for a subnetwork to be considered. Default is 0')
parser.add_argument('--nSubViz',type=int,help='First N subnetworks to visualize. Default is 200')
#parser.add_argument('--out',type=int,help='Path for saving the network visualization')

args = parser.parse_args()

loci_diction,fa_set,fa_diction,fa_string_set,string_set = data_to_diction(input_path=args.diseaseloci,string_path=args.network)
fa_subnetwork,gene_scores_loci,iteration_fa_genes_diction_2= subnetwork_genescore(loci_diction,num_subnet=args.nSub,min_edges=args.minEdges)
gene_scores_final=gene_scores_avg(loci_diction,gene_scores_loci,compare_set=string_set,num_subnet=args.nSub)
file_path=fa_subnet_viz(gene_scores_final,fa_subnetwork,loci_diction,file_path="./network.png",num_subs_viz=args.nSubViz)

# save the 5000 random FA subnetworks
with open('FA_subnetwork.txt', 'w') as file:
    for k, v in fa_subnetwork.items():
        for sub in v:
            file.write("{}\t{}\n".format(k,'\t'.join(sub)))

# run time
end = (time.time())
print(end-start)