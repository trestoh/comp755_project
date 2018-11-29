import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import urllib
import os.path
import sys
from sklearn import linear_model
from data_import import percolator_import, summary_stats

try:
    import cPickle as pickle
    kwargs = {}
except:
    import _pickle as pickle
    kwargs = {'encoding':'bytes'}
    
import gzip

data, true_hits, peptide_number = percolator_import(sys.argv[1])

data2, true_hits_2, peptide_number_2 = percolator_import(sys.argv[2])

means, sdevs = summary_stats(data)

pre_trans_data = data

pre_trans_data2 = data2

#save scans and labels
scan_nums = data[:,0]
labels = data[:,1]

scan_nums_2 = data2[:,0]
labels_2 = data2[:,1]
data = data [:, [i for i in range(2, 20)]]
data2 = data2 [:, [i for i in range(2, 20)]]

means, sdevs = summary_stats(data)

for i in range (0, len(data[0,:])):
    data[:,i] = (data[:,i] - means[i]) / sdevs[i]

for i in range (0, len(data2[0,:])):
    data2[:,i] = (data2[:,i] - means[i]) / sdevs[i]

data[:,[7,12,13]] = 0
data2[:,[7,12,13]] = 0

classes = [-1,1]

sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True)

#train on first N scans of data1
count = 0
N_scans = 2000
for scan in scan_nums:
    if int(scan) > N_scans:
        break
    #init_set.append(data[count,:])
    count += 1 

#labels.reshape(-1,1)

#subset for training
init_set = data[[i for i in range(0,count)], :]
init_labels = labels[[i for i in range(0,count)]]

#print (init_set)
#print (init_labels)

#fit on subset
sgd_clf.fit(init_set, init_labels)

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

#print(pre_trans_data[count])

#print(pre_trans_data[count])


#evaluate results
current_scan = int(pre_trans_data[count][0])
last_scan = 0

true_true_hits = 0
false_true_hits = 0
expected_true_hits = 0

while (count < len(data)):
    current_scan = int(pre_trans_data[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        if true_hits.get(current_scan) is not None:
            expected_true_hits += 1
        pred = sgd_clf.predict(data[count, :].reshape(1,-1))
        #print(int(pre_trans_data[count][1]), pred[0])
        if int(pre_trans_data[count][1]) == int(pred[0]):
            if int(pred[0]) == 1:
                if true_hits.get(current_scan) == peptide_number[int(pre_trans_data[count][20])]:
                    true_true_hits += 1
                else:
                    false_true_hits += 1
                #print (peptide_number[int(pre_trans_data[count][20])])
                true_pos += 1
            else:
                true_neg += 1
        else:
            if int(pred[0]) == 1:
                false_pos += 1
            else:
                false_neg += 1
    count += 1

print ("True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (true_pos, true_neg, false_pos, false_neg) )
print ("Final true pos: %d, final false neg: %d, expected true hits: %d" % (true_true_hits, false_true_hits, expected_true_hits))

print ("Full training on set 1 then running on set 2")

#
#
#   Run below for training partially
#
#

curr_scan = int(scan_nums[0])
cand_per_scan = 0
include = []

num_cands = 6

for i in range(0, len(data)):
    if int(scan_nums[i]) != curr_scan:
        cand_per_scan = 0
    curr_scan = int(scan_nums[i])
    if cand_per_scan < num_cands:
        include.append(i)

init_set = data[include, :]
init_labels = labels[include]

#print (init_set)
#print (init_labels)

#stable classifier
#sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True, n_iter=1000, shuffle=True)

sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True, n_iter=1000, shuffle=True, penalty='none')

#first is training on all, second is only training on top by some xcorr
sgd_clf.fit(data, labels)
#sgd_clf.fit(init_set, init_labels)

count = 0

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

current_scan = int(pre_trans_data2[count][0])
last_scan = 0

true_true_hits = 0
false_true_hits = 0
expected_true_hits = 0

while (count < len(data2)):
    current_scan = int(pre_trans_data2[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        if true_hits_2.get(current_scan) is not None:
            expected_true_hits += 1
        pred = sgd_clf.predict(data2[count, :].reshape(1,-1))
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        if int(pre_trans_data2[count][1]) == int(pred[0]):
            if int(pred[0]) == 1:
                if true_hits_2.get(current_scan) == peptide_number_2[int(pre_trans_data2[count][20])]:
                    true_true_hits += 1
                else:
                    false_true_hits += 1
                #print (peptide_number[int(pre_trans_data[count][20])])
                true_pos += 1
            else:
                true_neg += 1
        else:
            if int(pred[0]) == 1:
                false_pos += 1
            else:
                false_neg += 1
    count += 1


static_true_pos = 0
static_true_neg = 0
static_false_pos = 0 
static_false_neg = 0

count = 0

current_scan = int(pre_trans_data2[count][0])
last_scan = 0

while (count < len(data2)):
    current_scan = int(pre_trans_data2[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        #if true_hits_2.get(current_scan) is not None:
        #    expected_true_hits += 1
        xcorr = pre_trans_data2[count][5]
        charge2 = pre_trans_data2[count][10]
        if int(charge2) == 1:
            thresh = 1.3
        else:
            thresh = 1.5
        #print(int(pre_trans_data[count][1]), pred[0])

        if xcorr > thresh:
            if true_hits_2.get(current_scan) == peptide_number_2[int(pre_trans_data2[count][20])]:
                static_true_pos += 1
            else: 
                static_false_pos += 1
        else:
            if true_hits_2.get(current_scan) == peptide_number_2[int(pre_trans_data2[count][20])]:
                static_false_neg += 1
            else:
                static_true_neg += 1

    count += 1


print ("True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (true_pos, true_neg, false_pos, false_neg) )
print ("Final true pos: %d, final false positive: %d, expected true hits: %d" % (true_true_hits, false_true_hits, expected_true_hits))


print ("Static: True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (static_true_pos, static_true_neg, static_false_pos, static_false_neg) )

#below is pred method with a score
#pred = sgd_clf.decision_function(data2[count, :].reshape(1,-1))

#below we are redoing true pos etc. using decision function
count = 0

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

current_scan = int(pre_trans_data2[count][0])
last_scan = int(pre_trans_data2[count][0])

true_true_hits = 0
false_true_hits = 0
expected_true_hits = 0

thresh = 1.0

#
# Below looks at top overall match
#
#

'''
best_score = -90.0
best = pre_trans_data2[count]

while (count < len(data2)):
    current = pre_trans_data2[count]
    current_scan = int(current[0])
    if current_scan == last_scan:
        pred = sgd_clf.decision_function(data2[count, :].reshape(1,-1))
        if pred > best_score:
            best_score = pred
            best = current
    else:
        last_scan = current_scan
        if true_hits_2.get(int(best[0])) is not None:
            expected_true_hits += 1
        pred = best_score
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        if pred > thresh:
            if int(best[1]) == 1:
                true_pos += 1
                if true_hits_2.get(int(best[0])) == peptide_number_2[int(best[20])]:
                    true_true_hits += 1
                else:
                    false_true_hits += 1
            else:
                false_pos += 1
        else:
            if int(best[1]) == 1:
                false_neg += 1
            else:
                true_neg += 1
        best = current
        best_score = sgd_clf.decision_function(data2[count, :].reshape(1,-1))

    count += 1
'''

#
#
# Below only looks at top Xcorr
#
#

while (count < len(data2)):
    current_scan = int(pre_trans_data2[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        if true_hits_2.get(current_scan) is not None:
            expected_true_hits += 1
        pred = sgd_clf.decision_function(data2[count, :].reshape(1,-1))
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        if pred > thresh:
            if int(pre_trans_data2[count][1]) == 1:
                true_pos += 1
                if true_hits_2.get(current_scan) == peptide_number_2[int(pre_trans_data2[count][20])]:
                    true_true_hits += 1
                else:
                    false_true_hits += 1
            else:
                false_pos += 1
        else:
            if int(pre_trans_data2[count][1]) == 1:
                false_neg += 1
            else:
                true_neg += 1

    count += 1


print ("True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (true_pos, true_neg, false_pos, false_neg) )
print ("Final true pos: %d, final false positive: %d, expected true hits: %d" % (true_true_hits, false_true_hits, expected_true_hits))


with open("previous_exp_train_1.0.txt", "a") as myfile:
    myfile.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (sys.argv[1], sys.argv[2], static_true_pos, static_false_pos, true_true_hits, false_true_hits, expected_true_hits) )

