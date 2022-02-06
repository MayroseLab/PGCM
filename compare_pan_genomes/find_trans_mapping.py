from __future__ import print_function, division
import sys
import os
import pandas as pd

paf_dir = sys.argv[1]
out_tsv = sys.argv[2]

def format_mapping(row):
  if row["Target_sequence_name"] == '*':
    return ''
  return "%s:%s-%s;%s" %(row["Target_sequence_name"],row["Target_start"],row["Target_end"],round(row["Number_of_residue_matches"]/row["Query_sequence_length"],4))

paf_colnames = ["Query_sequence_name", "Query_sequence_length", "Query_start", "Query_end", "strand", "Target_sequence_name", "Target_sequence_length", "Target_start", "Target_end", "Number_of_residue_matches", "Alignment_block_length", "Mapping_quality"]
mapping_df = pd.DataFrame()
for paf in os.listdir(paf_dir):
  if not paf.endswith('.paf'):
    continue
  sample = os.path.splitext(paf)[0]
  paf_path = os.path.join(paf_dir, paf)
  paf_df = pd.read_csv(paf_path, sep='\t', usecols=range(12), names=paf_colnames)
  paf_df['sample'] = sample
  # calculate fraction of matches
  paf_df['frac_matches'] = paf_df["Number_of_residue_matches"]/paf_df["Query_sequence_length"]
  # only keep best hit (based on fraction of matches)
  paf_df = paf_df.loc[paf_df.groupby('Query_sequence_name')['frac_matches'].idxmax()]
  # format mapping
  paf_df['Mapping'] = paf_df.apply(format_mapping, axis=1)
  # Only keep needed columns
  paf_df = paf_df[["Query_sequence_name","sample","Mapping"]]
  # add
  mapping_df = pd.concat([mapping_df, paf_df])

# pivot
mapping_df = mapping_df.pivot_table(columns="sample",values="Mapping", index="Query_sequence_name", aggfunc='first')
# write
mapping_df.to_csv(out_tsv, sep='\t')
