#!/bin/bash
#install ete3

import os
import sys
from ete3 import Tree, TreeStyle, TextFace,NodeStyle,faces, COLOR_SCHEMES


WD=os.getcwd()
FILENAMES=sys.argv[1]
t=sys.argv[2]

#Species Tree
sptree_file = WD+"/4_coalescent_trees/"+FILENAMES+t+"/"+FILENAMES+t+"_coalescent.tre"

#Phyparts Output Files
histfile = WD+"/5_phypartstopiecharts/"+FILENAMES+t+"/"+FILENAMES+t+".hist"
keyfile = WD+"/5_phypartstopiecharts/"+FILENAMES+t+"/"+FILENAMES+t+".node.key"
concon_treefile = WD+"/5_phypartstopiecharts/"+FILENAMES+t+".concon.tre"
total_genes = 176

ptree = Tree(sptree_file)
sptree.convert_to_ultrametric()

def bootstrap_fn(mynode):
    if not mynode.is_leaf():
        F = TextFace(str(mynode.support))
        faces.add_face_to_node(F,mynode,0,"branch-top")

        

bs = TreeStyle()
bs.layout_fn = bootstrap_fn
bs.mode="r"
bs.show_leaf_name = True        


sptree.render("%%inline",tree_style=bs)