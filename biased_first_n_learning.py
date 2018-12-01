import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import urllib
import os.path
import sys
import random
from datetime import datetime
from sklearn import linear_model
from data_import import percolator_import, summary_stats

data, true_hits, peptide_number = percolator_import(sys.argv[1])
pre_trans_data = data


scan_nums = data[:,0]
orig_labels = data[:,1]
labels = orig_labels

#strip scans and peps out 
data = data [:, [i for i in range(2, 20)]]

#
# Creat First N Scan Subset
#
count = 0
N_scans = 50000
for scan in scan_nums:
    if int(scan) > N_scans:
        break
    #init_set.append(data[count,:])
    count += 1 

end_train = count
#labels.reshape(-1,1)

#subset for training
init_set = data[[i for i in range(0,count)], :]


init_orig_labels = labels[[i for i in range(0,count)]]
init_labels = init_orig_labels

print(init_labels)

decoy_indexes = np.where(init_labels == -1)[0]
print (decoy_indexes[0])

np.random.shuffle(decoy_indexes)
print (decoy_indexes[0])

init_scans = scan_nums[[i for i in range(0,count)]]
pre_trans_init = pre_trans_data[[i for i in range(0,count)], :]

#
# Compute summary data over First N subset
#
means, sdevs = summary_stats(init_set)

#
# Scale data
#
for i in range (0, len(data[0,:])):
    data[:,i] = (data[:,i] - means[i]) / sdevs[i]

#
# Features are Always Zero, set them as such
#
data[:,[7,12,13]] = 0

#
# Reset Init Set to get scaling and zeroed columns
#
init_set = data[[i for i in range(0,count)], :]


#
# Build initial positives from training set
#
xcorr_col = []
xcorr_lim = 0.5

for i in range(0, len(init_set)):
    xcorr = pre_trans_init[i][5]
    label = int(init_labels[i])
    xcorr_col.append((i, label, xcorr))

pass_thresh = []

false_discovery = .3

while True:
    pass_thresh = []
    for (number, label, xcorr) in xcorr_col:
        if xcorr > xcorr_lim:
            pass_thresh.append(label)
    true_pos = sum(label == 1 for label in pass_thresh)
    false_pos = sum(label == -1 for label in pass_thresh)
    if 2 * false_pos / (true_pos + false_pos) < false_discovery:
        break
    xcorr_lim += .01

print ("XCorr lim = %f" % xcorr_lim)

training_set = []

for i in range(0, len(init_set)):
    xcorr = pre_trans_init[i][5]
    if int(init_orig_labels[i]) == 1 and xcorr > xcorr_lim:
        init_labels[i] = 1
        training_set.append(i)
    else:
        init_labels[i] = -1

random.seed(datetime.now())

print ("Positive examples in first %d scans: %d" % (N_scans, len(training_set)))

pos_len = len(training_set)
count = 0
iterations = 0

np.random.shuffle(decoy_indexes)
for index in range(0,pos_len):
    training_set.append(decoy_indexes[index])

#while True:
#    if int(init_orig_labels[count]) == -1 and init_labels[count] == -1:
#        flip = random.randint(0, 1)
#        if flip == 1:
#            training_set.append(count)
#        if len(training_set) >= 2 * pos_len:
#            break
#    if count == len(init_set) - 1:
#        count = 0
#        iterations += 1
#        if iterations >= 100:
#            break
#    else:
#        count += 1

print ("Total examples in first %d scans: %d" % (N_scans, len(training_set)))


classes = [-1,1]
sgd_clf = linear_model.SGDClassifier(loss="hinge", max_iter=1000, tol=1e-3, shuffle=True)

init_train_set = init_set[training_set, :]
init_train_labels = init_labels[training_set]

sgd_clf.fit(init_train_set, init_train_labels)

print ("1 Iteration Done")

for j in range (0, 9):
    training_set = []

    thresh = 0.0

    predictions = []
    pass_thresh = []

    count = 0
    current_scan = int(pre_trans_init[count][0])
    last_scan = 0

    false_discovery = 0.1


    #
    # Consider everything, not just top
    #
    while (count < len(init_set)):
        current_scan = int(pre_trans_init[count][0])
        #uncomment and tab for less consideration
        if current_scan != last_scan:
            last_scan = current_scan
            pred = sgd_clf.decision_function(init_set[count, :].reshape(1,-1))
            #print (pred[0])
            #print(int(pre_trans_data[count][1]), pred[0])
            predictions.append((int(pre_trans_init[count][1]), pred[0]))

        count += 1

    for (label, prediction) in predictions:
        if prediction > thresh:
            pass_thresh.append(label)
    true_pos = sum(label == 1 for label in pass_thresh)
    false_pos = sum(label == -1 for label in pass_thresh)

    if 2 * false_pos / (true_pos + false_pos) >= false_discovery:
        increase = True
    else:
        increase = False

    while True:
        pass_thresh = []
        for (label, prediction) in predictions:
            if prediction > thresh:
                pass_thresh.append(label)
        true_pos = sum(label == 1 for label in pass_thresh)
        false_pos = sum(label == -1 for label in pass_thresh)
        if increase:
            if 2 * false_pos / (true_pos + false_pos) <= false_discovery:
                break
            thresh += .01
        else:
            if 2 * false_pos / (true_pos + false_pos) >= false_discovery:
                break
            thresh -= .01
    
    print("Thresh = %f" % thresh)

    for i in range(0, len(init_set)):
        pred = sgd_clf.decision_function(data[i, :].reshape(1,-1))
        if pred > thresh  and init_orig_labels[i] == 1:
            init_labels[i] = 1
            training_set.append(i)
        else:
            init_labels[i] = -1
    
    random.seed(datetime.now())

    pos_len = len(training_set)
    print ("Number of positives in train set: %d" % pos_len)
    
    np.random.shuffle(decoy_indexes)
    for index in range(0,pos_len):
        training_set.append(decoy_indexes[index])

    sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True, max_iter=1000, tol=1e-3, shuffle=True)

    init_train_set = init_set[training_set, :]
    init_train_labels = init_labels[training_set]

    sgd_clf.fit(init_train_set, init_train_labels)
    print ("%d Iteration Done" % (j + 2))

static_true_pos = 0
static_true_neg = 0
static_false_pos = 0 
static_false_neg = 0

static_true_true = 0
static_false_true = 0

count = end_train

current_scan = int(pre_trans_data[count][0])
last_scan = 0

#print (true_hits)


while (count < len(data)):
    current_scan = int(pre_trans_data[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        #if true_hits_2.get(current_scan) is not None:
        #    expected_true_hits += 1
        xcorr = pre_trans_data[count][5]
        charge2 = pre_trans_data[count][10]
        if int(charge2) == 1:
            thresh = 1.3
        else:
            thresh = 1.5
        #print(int(pre_trans_data[count][1]), pred[0])

        if xcorr > thresh:
            #print ("Above Thresh")
            #print (true_hits.get(current_scan))
            if int(pre_trans_data[count][1]) == 1:
                static_true_pos += 1
                if true_hits.get(current_scan) == peptide_number[int(pre_trans_data[count][20])]:
                    static_true_true += 1
                else:
                    static_false_true += 1
            else:
                static_false_pos += 1
        else:
            if int(pre_trans_data[count][1]) == 1:
                static_false_neg += 1
            else:
                static_true_neg += 1
    count += 1

print ("Static: True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (static_true_pos, static_true_neg, static_false_pos, static_false_neg) )
print ("Static: Final true pos: %d, Final false pos: %d" % (static_true_true, static_false_true) )

count = end_train

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0


true_true_hits = 0
false_true_hits = 0
expected_true_hits = 0

thresh = 0.0

predictions = []
pass_thresh = []

count = 0
current_scan = int(pre_trans_init[count][0])
last_scan = 0

false_discovery = 0.1

while (count < len(init_set)):
    current_scan = int(pre_trans_init[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        pred = sgd_clf.decision_function(init_set[count, :].reshape(1,-1))
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        predictions.append((int(pre_trans_init[count][1]), pred[0]))

    count += 1

for (label, prediction) in predictions:
    if prediction > thresh:
        pass_thresh.append(label)
true_pos = sum(label == 1 for label in pass_thresh)
false_pos = sum(label == -1 for label in pass_thresh)

if 2 * false_pos / (true_pos + false_pos) >= false_discovery:
    increase = True
else:
    increase = False

#print (2 * false_pos / (true_pos + false_pos))

while True:
    pass_thresh = []
    for (label, prediction) in predictions:
        if prediction > thresh:
            pass_thresh.append(label)
    true_pos = sum(label == 1 for label in pass_thresh)
    false_pos = sum(label == -1 for label in pass_thresh)
    if increase:
        #print (2 * false_pos / (true_pos + false_pos))
        if 2 * false_pos / (true_pos + false_pos) <= false_discovery:
            break
        thresh += .01
    else:
        #print (2 * false_pos / (true_pos + false_pos))
        if 2 * false_pos / (true_pos + false_pos) >= false_discovery:
            break
        thresh -= .01
    #print (2 * false_pos / (true_pos + false_pos))
    #print (thresh)

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

#thresh = -3.0

print("Thresh = %f" % thresh)
count = end_train

fdr_file = open("first_50k_30_10_fdr.txt", "a")

while (count < len(data)):
    current_scan = int(pre_trans_data[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        pred = sgd_clf.decision_function(data[count, :].reshape(1,-1))

        fdr_file.write("%d\t%f\n" % (int(pre_trans_data[count][1]), pred) )

        if true_hits.get(current_scan) is not None:
            expected_true_hits += 1

        if pred > thresh:
            if int(pre_trans_data[count][1]) == 1:
                true_pos += 1
                if true_hits.get(current_scan) == peptide_number[int(pre_trans_data[count][20])]:
                    true_true_hits += 1
                else:
                    false_true_hits += 1
            else:
                false_pos += 1
        else:
            if int(pre_trans_data[count][1]) == 1:
                false_neg += 1
            else:
                true_neg += 1
    count += 1


print ("True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (true_pos, true_neg, false_pos, false_neg) )
print ("Final true pos: %d, final false positive: %d, expected true hits: %d" % (true_true_hits, false_true_hits, expected_true_hits))

with open("restart_biased_first_30k_20fdr_in_10_out.txt", "a") as myfile:
    myfile.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (sys.argv[1], static_true_pos, static_false_pos, true_true_hits, false_true_hits, expected_true_hits) )

fdr_file.close()