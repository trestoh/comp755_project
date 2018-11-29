import numpy as np

def percolator_import(file):
    filename = "/Users/trentstohrer/Desktop/UNC CS/COMP755/final_project/" + file + "_target.pin"
    decoy_file = "/Users/trentstohrer/Desktop/UNC CS/COMP755/final_project/" + file + "_decoy.pin"
    true_hit_file = "/Users/trentstohrer/Desktop/UNC CS/COMP755/final_project/" + file + "_percolator_psms.txt"

    with open(filename) as f:
        content = f.readlines()

    with open(decoy_file) as f:
        content += f.readlines()


    with open(true_hit_file) as f:
        true_hit_text = f.readlines()

    count = 0

    data = []

    peptides = {}
    peptide_number = {}

    pept_count = 0

    for line in content:
        
        lines = line.split("\t")

        if lines[0] != "SpecId":

            stripped_peps = lines[23][2:(len(lines[23])-2)]

            if not peptides.get(stripped_peps):
                peptides[stripped_peps] = pept_count
                peptide_number[pept_count] = stripped_peps
                pep = pept_count
                pept_count += 1

            else:
                pep = peptides[stripped_peps]

            dat = [int(lines[2]), int(lines[1]), float(lines[5]), float(lines[6]), float(lines[7]), float(lines[8]), float(lines[9]), float(lines[10]), 
                int(lines[11]), int(lines[12]), int(lines[13]), int(lines[14]), int(lines[15]), int(lines[16]), int(lines[17]),
                int(lines[18]), int(lines[19]), float(lines[20]), float(lines[21]), float(lines[22]), pep]
            data.append(dat)
        count += 1

    data = np.array(data)

    data = data[np.argsort(data[:,5])][::-1]
    data = data[np.argsort(data[:,0], kind='mergesort')]

    true_hits = {}

    for line in true_hit_text:
        
        lines = line.split("\t")
        if lines[0] != "file_idx":
            #print ("q = %s" % lines[7])
            if float(lines[7]) > .01:
                break
            true_hits[int(lines[1])] = lines[10]
            #print (true_hits[lines[1]])

    return (data, true_hits, peptide_number)

def summary_stats (data):

    means = []
    sdevs = []

    for i in range (0, len(data[0,:])):
        means.append(np.mean(data[:,i]))
        sdevs.append(np.std(data[:,i]))



    print (means)
    print (sdevs)
    return (means, sdevs)