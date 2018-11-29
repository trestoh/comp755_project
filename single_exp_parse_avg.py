import string
import sys
import numpy as np

with open(sys.argv[1]) as f:
        content = f.readlines()

true_pos_pcts = []
false_pos_pcts = []

for line in content:
    vals = line.split('\t')
    pct_true_pos = float(vals[3]) / float(vals[1])
    pct_false_pos = float(vals[4]) / float(vals[2])
    true_pos_pcts.append(pct_true_pos)
    false_pos_pcts.append(pct_false_pos)


print (np.mean(true_pos_pcts))
print (np.std(true_pos_pcts))
print (np.mean(false_pos_pcts))
print (np.std(false_pos_pcts))