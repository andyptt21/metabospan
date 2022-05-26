## Create dictionary
import pdb
import json
import numpy as np
import pandas as pd

with open("/Users/pattac/Desktop/KEGG/compound/compound") as file:
    file_contents = file.read()

with open("/Users/pattac/Desktop/KEGG/parsed_hsa_only_map_numbers.txt") as file:
    human_pathway_numbers = file.read()

human_pathway_numbers = human_pathway_numbers.split("\n")
human_pathway_ids = []
for i in human_pathway_numbers:
    human_pathway_ids.append("map"+i)

file_split = file_contents.split("///")

def splitter():
    my_dict = {}
    for i in file_split:
        lines = i.split("\n")
        str1 = "ENTRY"
        compound = [i for i in lines if str1 in i]
        if isinstance(compound, list):
            if len(compound) == 1:
                compound_name = compound[0].split()[1]
                str2 = "map"
                pathways = [i for i in lines if str2 in i]
                pathways_split_list = list()
                for temp in pathways:
                    pathways_split = temp.split()
                    pathways_split = [temp for temp in pathways_split if str2 in temp]
                    if (pathways_split[0] in human_pathway_ids):
                        pathways_split_list.append(pathways_split)
                my_dict[compound_name] = pathways_split_list
            else:
                continue
        else:
            print(compound)
            break
                
    return my_dict

kegg_pathways = splitter()
kegg_pathways = {k: v for k, v in kegg_pathways.items() if v}

## Build adjacency matrix

def getList(dict):
    return dict.keys()

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

compounds = getList(kegg_pathways)

adjacency_matrix = np.array([0]*len(compounds)**2).reshape(len(compounds),len(compounds))
adjacency_matrix = pd.DataFrame(adjacency_matrix,index = compounds, columns = compounds)

for i in list(compounds):
    for j in list(compounds):
        cpd1 = [''.join(ele) for ele in kegg_pathways[i]]
        cpd2 = [''.join(ele) for ele in kegg_pathways[j]]
        overlap = jaccard(tuple(cpd1),tuple(cpd2))
        adjacency_matrix.loc[i,j] = overlap

adjacency_matrix.to_csv("/Users/pattac/Documents/metabospan/data-raw/kegg_network/kegg_overlap.csv")

