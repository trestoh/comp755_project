This read_me is to give an idea of what the various scripts do and what the various files are for, how to invoke them and how to replicate results.

I find script writing much more natural than the jupyter notebook style so this is how I've decided to do things. You will find a LaTeX formatted formal writeup in this repository as well.

The descriptions are in alphabetical order, and though this isn't the order I wrote the scripts, it is the order I wrote the read me.

The most important thing to note about all the scripts is that they assume python3. Also note that all scripts in the results are "biased" in spite of the somewhat confusing naming.

All invocations are bash/shell invocations.

Notes on output
---------------
Most console ouptut are developer notes. The main console output of importance is at the end.

We have the same outputs for both the static threshold, and the classfier

The differentiation between "true positives" and "final true positives" relate to Target Matches, those that are considered matches by the algorithm we are trying to approximate, and Targets, which are just real peptide sequences as opposed to the jumbled up ones.

Then we have two seperate breakdowns for the two.

First is True/False Pos/Neg
These are the counts of Targets identified as Targets (true pos), Decoys identified as Decoys (true neg), Decoys identified as Targets (false pos), and Targets identified as Decoys (false neg).

Then we have a breakdown of true pos into Target Matches identified as Target Matches (final true pos) and non-matching Targets identified as Target Matches (final false pos), as a decoy would never actually be considered a match since we know it is definitely a decoy. 

The expected true hits are the number of Target Matches over the test set from the classifier we're trying to approximate.

Finally, there are also two sets of file outputs that occur: labels and scores that allow us to make our modified ROC curves and comparisons of static and classifier final true pos/neg. This second set is not really used for analysis and is more or less deprecated at this point.


MA1750_decoy.pin/MA1750_target.pin
----------------------------------
Example input data. The other set used in the report was too large for github so it is not here. All input is formatted as [experiment]_[specific_file].[ext] and the scripts rely on everything after [experiment] being the same for all potential input.

MA1750_percolator_psms.txt
--------------------------
These are the results we are trying to approximate and where the "true matches" come from. The relevant columns are scan number, q score (the false discovery rate at which the match is "discovered", they report 1% so those are the targets we consider "true matches"), and the peptide identified.

biased_first_n_learning.py
--------------------------
This script trains on the first N scans of the given experiment, N being statically specified and changed in the script itself. Each scan has around 5 targets and 5 decoys. It then tests that classifier on the rest of the scans without subsequent training. The specifics of training are more or less the same for all these scripts, and specified in more detail in the writeup. The biased in the name is in reference to the fact that we get our initial labels and training set using a single variable as the only selection criteria.

Invocation: python3 biased_first_n_learning.py [experiment]

biased_train_on_prev.py
-----------------------
This script trains on one experiment and then tests on a second. I am writing this read me before the writeup so it may or may not be included in the report. As of this writing we do not control for false discovery training, but use a threshold of 0.0 until all training iterations are finished, then pick a controlled threshold. Also note that in testing this script I used experiments that are not in this repo, I can give them to you if you like.

Invocation: python3 biased_train_on_prev.py [training_expeirment] [test_experiment]

data_import.py
--------------
Helper functions used by all scripts. Contains functions for importing and sorting input data (to simulate acquisition over time) and the "ideal results" in the form of the output from the algorithm we're trying to approximate. Also contains a function to aid in standardizing data that could as easily be a one off.

Invocation: N/A

first_n.sh
----------
Unix shell script to call biased_first_n_learning.py on many experimental data sets. Note almost all of these are not in the repo.

Deprecated

Invocation: ./first_n.sh

first_n_total_retrain.py
------------------------
A refinement on first_n_with_online.py, full details in the writeup. Trains repeatedly on the first N scans, then retrains once every batch size number of scans using all (N + i * batch_size) scans. The script tests using whatever the current model and threshold are, which are redone every batch size scans. 

Invocation: python3 first_n_total_retrain.py [experiment]

first_n_with_online.py
----------------------
Adds simple online learning to biased_first_n_learning.py. This script does the same initial learning as biased_first_n but then does a partial fit using any positively identified example targets in each batch size along with an equivalent number of decly examples from the same batch. Uses a static threshold for reasons explained in the writeup. 

Invocation: python3 first_n_with_online.py [experiment]

generate_fdr.py
---------------
Used to generate data for modified ROC curves using the label and score data put out by other scripts.

Invocation: email me if you really need this, I don't think it is useful to anyone but me personally

main.py
-------
Deprecated. Early attempts are here, perhaps interesting if you want to see what's changed or early, silly bugs.

mass_comp.sh
------------
Unix shell script to call biased_train_on_prev.py on pairs of similar experimental data sets. Note almost all of these are not in the repo.

Deprecated

Invocation: ./mass_comp.sh

single_exp_parse_avg.py
-----------------------
Deprecated. Was used to analyze results across lots of experiments in a way that we decided was not useful.
