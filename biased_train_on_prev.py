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
data2, true_hits_2, peptide_number_2 = percolator_import(sys.argv[2])
pre_trans_data = data
pre_trans_data2 = data2

scan_nums = data[:,0]
orig_labels = data[:,1]
labels = orig_labels

scan_nums_2 = data2[:,0]
orig_labels_2 = data2[:,1]
labels_2 = orig_labels_2

decoy_indexes = np.where(labels == -1)[0]

np.random.shuffle(decoy_indexes)


#strip scans and peps out 
data = data [:, [i for i in range(2, 20)]]
data2 = data2 [:, [i for i in range(2, 20)]]

means, sdevs = summary_stats(data)

for i in range (0, len(data[0,:])):
    data[:,i] = (data[:,i] - means[i]) / sdevs[i]

for i in range (0, len(data2[0,:])):
    data2[:,i] = (data2[:,i] - means[i]) / sdevs[i]

data[:,[7,12,13]] = 0
data2[:,[7,12,13]] = 0

xcorr_col = []
xcorr_lim = 0.5

for i in range(0, len(data)):
    xcorr = pre_trans_data[i][5]
    label = int(orig_labels[i])
    xcorr_col.append((i, label, xcorr))

pass_thresh = []

false_discovery = .30

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

for i in range(0, len(data)):
    xcorr = pre_trans_data[i][5]
    if int(orig_labels[i]) == 1 and xcorr > xcorr_lim:
        labels[i] = 1
        training_set.append(i)
    else:
        labels[i] = -1

random.seed(datetime.now())

pos_len = len(training_set)
print ("Number of positives in train set: %d" % pos_len)
count = 0
iterations = 0
np.random.shuffle(decoy_indexes)
for index in range(0,pos_len):
    training_set.append(decoy_indexes[index])

classes = [-1,1]

#from initial data run changed this to L2 penalty rather than none
sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True, max_iter=1000, tol=1e-3, shuffle=True)

init_set = data[training_set, :]
init_labels = labels[training_set]

sgd_clf.fit(init_set, init_labels)

print ("1 Iteration Done")

for j in range (0, 9):
    #previously I was not wiping the training set between iterations
    training_set = []
    for i in range(0, len(data)):
        pred = sgd_clf.predict(data[i, :].reshape(1,-1))
        if pred[0] == 1  and orig_labels[i] == 1:
            labels[i] = 1
            training_set.append(i)
        else:
            labels[i] = -1
    
    random.seed(datetime.now())

    pos_len = len(training_set)
    print ("Number of positives in train set: %d" % pos_len)
    
    np.random.shuffle(decoy_indexes)
    for index in range(0,pos_len):
        training_set.append(decoy_indexes[index])

    sgd_clf = linear_model.SGDClassifier(loss="hinge", warm_start=True, max_iter=1000, tol=1e-3, shuffle=True)

    init_set = data[training_set, :]
    init_labels = labels[training_set]

    sgd_clf.fit(init_set, init_labels)
    print ("%d Iteration Done" % j)

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

print ("Static: True pos: %d, True neg: %d, False pos: %d, False neg: %d" % (static_true_pos, static_true_neg, static_false_pos, static_false_neg) )

count = 0

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

current_scan = int(pre_trans_data[count][0])
last_scan = int(pre_trans_data[count][0])

true_true_hits = 0
false_true_hits = 0
expected_true_hits = 0

thresh = 0.0

predictions = []


#
#
# This rate is not being calculated properly!!!
# Need to redo to get threshold better!
# Maybe not?
#
#
'''
while (count < len(data)):
    current_scan = int(pre_trans_data[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        pred = sgd_clf.decision_function(data[count, :].reshape(1,-1))
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        predictions.append((int(pre_trans_data[count][1]), pred[0]))

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
        '''

thresh = 0.0

predictions = []
pass_thresh = []

count = 0
current_scan = int(pre_trans_data[count][0])
last_scan = 0

false_discovery = 0.1

while (count < len(data)):
    current_scan = int(pre_trans_data[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        pred = sgd_clf.decision_function(data[count, :].reshape(1,-1))
        #print (pred[0])
        #print(int(pre_trans_data[count][1]), pred[0])
        predictions.append((int(pre_trans_data[count][1]), pred[0]))

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


#
# Just set thresh to 0.0
#
#thresh = 0.0

true_pos = 0
true_neg = 0
false_pos = 0 
false_neg = 0

thresh = -3.0

print("Thresh = %f" % thresh)
count = 0
while (count < len(data2)):
    current_scan = int(pre_trans_data2[count][0])
    if current_scan != last_scan:
        last_scan = current_scan
        pred = sgd_clf.decision_function(data2[count, :].reshape(1,-1))

        if (true_hits_2.get(current_scan)) is not None:
            expected_true_hits += 1

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

with open("biased_previous_exp_train_30_fdr_.txt", "a") as myfile:
    myfile.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (sys.argv[1], sys.argv[2], static_true_pos, static_false_pos, true_true_hits, false_true_hits, expected_true_hits) )
