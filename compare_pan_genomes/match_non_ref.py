"""
Match two sets of proteins based on BLAST results,
using the maximum weight bipartite matching method.
Parameters:
1. set1 proteins fasta
2. set2 proteins fasta
3. set1 vs. set2 blast6 result 
   (must use -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen")
4. set2 vs. set1 blast6 result
   (must use -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen")
5. Weight parameter - e.g. bitscore or pident
6. Min weight - minnimal weight value allowed in a match
* can normalize weights according to Emms, David M., and Steven Kelly, 2015

Output:
TSV file with the matches pairs and their respective weight
"""

from __future__ import print_function, division
from networkx import bipartite, matching
from networkx.algorithms.bipartite.matching import minimum_weight_full_matching
import argparse
import pandas as pd
import numpy as np

### FUNCTIONS
def assign_ids(fasta, start=0):
  i = start
  d = {}
  with open(fasta) as f:
    for line in f:
      if line.startswith('>'):
        d[i] = line.strip()[1:]
        i += 1
  return d
blast6_headers = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
def parse_blast6(blast6, weight_param):
  df = pd.read_csv(blast6, sep='\t', names=blast6_headers)
  # keep only best hit per pair
  df_top = df.loc[df.groupby(["qseqid", "sseqid"])[weight_param].idxmax()]
  return df_top

def normalize_score(blast6_df, weight_param, bin_size=500):
  # calculate Lqh
  blast6_df['length_product'] = blast6_df["qlen"] * blast6_df["slen"]
  # assign hits to equal-sized bins
  n_bins = int(blast6_df.shape[0]/bin_size)
  blast6_df['length_bin'] = pd.qcut(blast6_df['length_product'].rank(method='first'),q=n_bins, labels=False)
  # add normalized score column
  blast6_df['norm_'+weight_param] = 0
  
  bin_dfs = []
  for lbin in blast6_df['length_bin'].unique():
    bin_df = blast6_df.loc[blast6_df['length_bin'] == lbin].copy()
    # fetch the top 5% hits
    n_top_hits = int(0.05 * bin_df.shape[0])
    bin_df_top = bin_df.nlargest(n_top_hits, weight_param)
    # fit a linear model in log-log space
    bin_df_top['log_length_product'] = np.log10(bin_df_top['length_product'])
    bin_df_top['log_'+weight_param] = np.log10(bin_df_top[weight_param])
    a,b = np.polyfit(bin_df_top['log_length_product'], bin_df_top['log_'+weight_param], 1)
    # normalize scores according to model
    bin_df['norm_'+weight_param] = bin_df[weight_param]/(10**b * bin_df['length_product'] ** a)
    # also normalize by q/s coverage
    bin_df['qcov'] = (bin_df['qend']-bin_df['qstart'])/bin_df['qlen']
    bin_df['scov'] = (bin_df['send']-bin_df['sstart'])/bin_df['slen']
    bin_df['norm_'+weight_param] = bin_df['norm_'+weight_param] * bin_df['qcov'] * bin_df['scov']
    
    bin_dfs.append(bin_df)

  return pd.concat(bin_dfs)

def calculate_biderectional_weight(set1_vs_set2_weights, set2_vs_set1_weights, p1, p2):
  w1 = 0
  w2 = 0
  if (p1,p2) in set1_vs_set2_weights:
    w1 +=  set1_vs_set2_weights[(p1,p2)]
  if (p2,p1) in set2_vs_set1_weights:
    w2 += set2_vs_set1_weights[(p2,p1)]
  return -(w1+w2)/2

def df_to_dict(blast6_df, weight_param):
  """
  {(qseqid, sseqid): weight,...}
  """
  return {(row['qseqid'],row['sseqid']): row[weight_param] for ind,row in blast6_df.iterrows()}

### MAIN
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('set1_fasta', help="Proteins set 1 fasta")
  parser.add_argument('set2_fasta', help="Proteins set 2 fasta")
  parser.add_argument('set1_vs_set2_blast', help="Blast6 result of set1 vs set2")
  parser.add_argument('set2_vs_set1_blast', help="Blast6 result of set2 vs set1")
  parser.add_argument('out_tsv', help="Path to output file")
  parser.add_argument('--weight_param', default='bitscore', help="Blast column to use as edge weight" )
  parser.add_argument('--normalize_weight', default='False', action='store_true', help="Apply a normalization procedure to the weights to eliminate length bias" )
  parser.add_argument('--filter', default=None, help="A boolean expression to use for filtering blast hits. The expression will be evaluated by pd.df.query(). E.g: pident > 90 & qcov > 0.9 & scov > 0.9. Default: no filter")
  parser.add_argument('--set1_name', default="set1", help="Name of set1")
  parser.add_argument('--set2_name', default="set2", help="Name of set2")
  args = parser.parse_args()
  
  # go over fasta files and assign an id to each protein
  print("Parsing fasta files...")
  set1_proteins = assign_ids(args.set1_fasta)
  set1_size = len(set1_proteins)
  set2_proteins = assign_ids(args.set2_fasta, start=set1_size)
  set2_size = len(set2_proteins)
  
  # create bipartite graph
  bg = bipartite.complete_bipartite_graph(set1_size,set2_size)
  
  # parse blast6 files
  print("Parsing BLAST files...")
  set1_vs_set2_weights = parse_blast6(args.set1_vs_set2_blast, args.weight_param)
  set2_vs_set1_weights = parse_blast6(args.set2_vs_set1_blast, args.weight_param)

  # normalize weights
  if args.normalize_weight:
    print("Normalizing weights...")
    set1_vs_set2_weights = normalize_score(set1_vs_set2_weights, args.weight_param)
    set2_vs_set1_weights = normalize_score(set2_vs_set1_weights, args.weight_param)

  # filter blast hits
  if args.filter:
    print("Filtering BLAST hits...")
    set1_vs_set2_weights.query(args.filter, inplace=True)
    set2_vs_set1_weights.query(args.filter, inplace=True)
  
  # go over all BG edges and assign weights
  print("Assigning edge weights...")
  if args.normalize_weight:
    weight_param_name = 'norm_'+args.weight_param
  else:
    weight_param_name = args.weight_param
  set1_vs_set2_weights = df_to_dict(set1_vs_set2_weights, weight_param_name)
  set2_vs_set1_weights = df_to_dict(set2_vs_set1_weights, weight_param_name)
  for edge in bg.edges:
    q = set1_proteins[edge[0]]
    s = set2_proteins[edge[1]]
    bidirectional_weight = calculate_biderectional_weight(set1_vs_set2_weights, set2_vs_set1_weights, q, s)
    bg.edges[edge]['weight'] = bidirectional_weight
  
  # find max weight matching
  print("Finding max weight matching...")
  mw = minimum_weight_full_matching(bg)
  
  # print out matches
  print("Writing out results...")
  with open(args.out_tsv, 'w') as fo:
    print("%s\t%s\t%s weight" %(args.set1_name, args.set2_name, args.weight_param), file=fo)
    for match in mw:
      source, target = match, mw[match]
      if target < source:
        # each match appears in both directions in the output dict
        continue
      match_weight = -1 * bg[source][target]['weight']
      if match_weight > 0:
        q_name = set1_proteins[source]
        s_name = set2_proteins[target]
        print("%s\t%s\t%s" %(q_name, s_name, match_weight), file=fo)
