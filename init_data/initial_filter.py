import numpy as np
import pandas as pd
import pickle as pl
import pdb

raw = pd.read_csv("raw_big_data/big_data.ehr_stop_at_baseline.txt.gz", sep="\t", header=0)
colnames = raw.columns.to_numpy()
raw = raw.to_numpy()

count_na = []
for i in range(raw.shape[1]):
  count_na.append(sum(np.isnan(raw[:,i])))


count_low = []
for i in range(raw.shape[1]):
  vals = raw[np.logical_not(np.isnan(raw[:,i])),i]
  u_len = len(np.unique(vals))
  if u_len == 2:
    count_bin = []
    u_bin = np.unique(vals)
    for j in range(2):
      count_bin.append(sum(vals == u_bin[j]))
    count_low.append(min(count_bin))
  elif u_len < 2:
    count_low.append(0)
  else:
    count_low.append(100000)



#pdb.set_trace()
#count_na = pl.load( open( "count_na.p", "rb" ) )
#count_low = pl.load( open( "count_low.p", "rb" ) )

raw = raw[:,np.logical_and(np.array(count_low) > 300, np.array(count_na) < 100000)]
colnames = colnames[np.logical_and(np.array(count_low) > 300, np.array(count_na) < 100000)]

cor_matrix = pd.DataFrame(raw).corr(method = "pearson")

upper_tri = cor_matrix.where(np.triu(np.ones(cor_matrix.shape),k=1).astype(np.bool))
to_drop = [column for column in upper_tri.columns if any(upper_tri[column] > 0.95)]

raw = np.delete(raw, to_drop, 1)
colnames = np.delete(colnames, to_drop)

raw = pd.DataFrame(raw)
raw.columns = colnames

raw.to_csv("clean_big_data/init_clean_big_data.txt", sep = "\t")
print("done")
