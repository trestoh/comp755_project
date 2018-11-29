This read_me is to give an idea of what the various scripts do and what the various files are for, how to invoke them and how to replicate results.

I find script writing much more natural than the jupyter notebook style so this is how I've decided to do things. You will find a LaTeX formatted formal writeup in this repository as well.

The most important thing to note about all the scripts is that they assume python3.

MA1750_decoy.pin/MA1750_target.pin
----------------------------------
Example input data. The other set used in the report was too large for github so it is not here. All input is formatted as [experiment]_[specific_file].[ext] and the scripts rely on everything after [experiment] being the same for all potential input.

MA1750_percolator_psms.txt
--------------------------
These are the results we are trying to approximate and where the "true matches" come from. The relevant columns are scan number, q score (the false discovery rate at which the match is "discovered", they report 1% so those are the targets we consider "true matches"), and the peptide identified.

biased_first_n_learning.py
--------------------------

biased_train_on_prev.py
-----------------------

data_import.py
--------------

first_n.sh
----------

first_n_total_retrain.py
------------------------

first_n_with_online.py
----------------------

generate_fdr.py
---------------

main.py
-------
Deprecated. Early attempts are here, perhaps interesting if you want to see what's changed or early, silly bugs.

mass_comp.sh
------------

single_exp_parse_avg.py
-----------------------
Deprecated
