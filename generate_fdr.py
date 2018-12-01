import numpy as np
import sys

filename = sys.argv[1]

with open(filename) as f:
    content = f.readlines()

data = []

for line in content:    
    lines = line.split("\t")
    dat = [int(lines[0]), float(lines[1])]
    data.append(dat)

data = np.array(data)

data = data[np.argsort(data[:,1])]

print (data)

#with open("test_fdr_sort.txt", "a") as myfile:

true_pos = sum(int(dat[0]) == 1 for dat in data)
false_pos = sum(int(dat[0]) == -1 for dat in data)
pos_found = true_pos

false_discovery = (2 * false_pos / (pos_found + false_pos)) * 100.0
true_discovery = (pos_found/true_pos) * 100.0

with open("roc_test_5.txt", "a") as myfile:
    for i in range(0, len(data)):
        if int(data[i,0]) == -1:
            false_pos -= 1
        else:
            pos_found -= 1
        
        #false_discovery = (2 * false_pos / (pos_found + false_pos)) 
        #true_discovery = (pos_found/true_pos) 
        
        expected_true = pos_found - false_pos
        false_found = 2 * false_pos

        #print (false_discovery, true_discovery)
        #myfile.write("%f,%f\n" % (false_discovery, true_discovery))
        #myfile.write("%f,%f\n" % (false_discovery, pos_found))
        myfile.write("%f,%f,%f\n" % (false_found, expected_true, data[i,1]))



print (pos_found, false_pos)
print (false_discovery, true_discovery)

