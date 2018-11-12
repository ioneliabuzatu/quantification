import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator as o
import sys
import os
import numpy as np
from glob import glob


fig = plt.figure()


def barplot(blastn): #, reads_number):
    '''
    Create a barchart for data across different categories with
    multiple conditions for each category.

    @param ax: The plotting axes from matplotlib.
    @param dpoints: The data set as an (n, 3) numpy array
    '''

    ax = fig.add_subplot(111)
    #abs_blastn = os.path.abspath(blastn)
    #blastn = [y for x in os.walk(abs_blastn) for y in glob(os.path.join(x[0], '*.blastn'))]

    table_blast = {}

    all_species = []
    for file in blastn:
        with open(file, "r") as fh:
            for line in fh:
                specie = line.split("\t")
                if len(specie) > 1:
                    name = " ".join(specie[1].split(" ")[:2])
                    all_species.append(name.rstrip())
    all_species_uniq = list(set(all_species))

    for file in blastn:
        barcode = file.split("/")[-1].split(".")[0]
        table_blast[barcode] = []
        with open(file, "r") as fh:
            for line in fh:
                specie = line.split("\t")
                if len(specie) > 1:
                    name = " ".join(specie[1].split(" ")[:2])
                    table_blast[barcode] = table_blast[barcode] + [name.rstrip()]

    #complete = []
    complete_species = []
    for key in table_blast:
        uniq_key = list(set(table_blast[key]))
        species_not_present = list(set(all_species_uniq) - set(uniq_key))
        np_number = []
        for ns in species_not_present:
            complete_species.append([key, ns, 0])
        aa = [[x, table_blast[key].count(x)] for x in set(table_blast[key])]
        for line in aa:
            line = [key] + line
            complete_species.append(line)
    dpoints = np.array(complete_species)
    conditions = [(c, np.mean(dpoints[dpoints[:,0] == c][:,2].astype(float))) for c in np.unique(dpoints[:,0])]
    categories = [(c, np.mean(dpoints[dpoints[:,1] == c][:,2].astype(float))) for c in np.unique(dpoints[:,1])]
    conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
    categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]
    dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))
    space = 0.3
    n = len(conditions)
    width = (1 - space) / (len(conditions))
    for i, cond in enumerate(conditions):
        indeces = range(1, len(categories)+1)
        vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.int)
        pos = [j - (1 - space) / 2. + i * width for j in indeces]
        ax.bar(pos, vals, width=width, label=cond, align="center", color=cm.Accent((float(i)+1) / n))

    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=90)

    # Add the axis labels
    plt.yscale('log')
    plt.autoscale(enable=True, axis='y')
    ax.set_ylabel("RMSD")
    ax.set_xlabel("Structure")

    # Add a legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')

if __name__ == '__main__':
    barplot(*sys.argv[1:])
    plt.show()