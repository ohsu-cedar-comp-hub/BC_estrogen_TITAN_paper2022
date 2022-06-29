#!/usr/bin/env python

import pandas as pd
import numpy as np
import os, glob
import pickle

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.utils import load_motifs

import seaborn as sns

DATABASE_FOLDER = "<location_of_databases>"
RESOURCES_FOLDER = "<location_of_resources>"

DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38*.mc9nr.feather")

SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "MCF7_E2_expmat.csv")

ex_matrix = pd.read_csv(SC_EXP_FNAME, sep = ',', header = 0, index_col = 0)

print("trying to load databases")

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

#print("everything loaded. Running grnboost2")

#tf_names = load_tf_names("/home/groups/CEDAR/doe/projects/LDA/SCENIC_prac/TF.tsv")

#if __name__ == '__main__':

#       adjacencies = grnboost2(expression_data=ex_matrix, verbose = True, client_or_address=client)

#       adjacencies.to_csv("adjacent_output.tsv", index=False, sep = '\t')

print("calculating modules")

adjacencies = pd.read_csv("adj.tsv", sep = "\t")

modules=list(modules_from_adjacencies(adjacencies, ex_matrix))

with open("module_output.pickle", 'wb') as f:
    pickle.dump(modules, f)

#print("prune step")

#df = prune2df(dbs, modules, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

#df.to_csv("motif_SCENIC_output.csv")

print("regulon step")

df = load_motifs("prune_MCF7_E2_output.csv")

genesig = df['Enrichment']['TargetGenes']

goodlist = []
for i in range(len(genesig)):
    if len(genesig[i]) > 0:
        goodlist.append(i)

df = df.iloc[goodlist]

regulons = df2regulons(df)

with open("regulon_T4D_E2_output.pickle", 'wb') as f:
    pickle.dump(regulons, f)

with open("regulon_MCF7_E2_output.pickle", 'rb') as f:
    regulons = pickle.load(f)

auc_mtx = aucell(ex_matrix, regulons, num_workers = 10)

auc_mtx.to_csv("auc_MCF7_E2_mtx.csv")

cmap = sns.clustermap(auc_mtx, figsize = (12,12))

cmap.savefig("regulon_MCF7_E2_clustermap.png")
