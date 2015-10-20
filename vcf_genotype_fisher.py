# coding: utf-8
from pathlib import Path
import sys
import argparse
import yaml
import pandas as pd
import linecache
from progressbar import *
from datetime import datetime

# argument parse
parser = argparse.ArgumentParser(description="Input the list of control and patient vcf file in yaml type.")
parser.add_argument('-l','--list', type=str, help="the yaml file must put with all the vcf files in the same folder.", required=True, metavar=".yaml")
args = parser.parse_args()
list_path = args.list
p = Path(list_path).parent


# read the yaml file listed the file names of controls and patients
vcf_list = []
with open(list_path) as f:
    file_yaml = yaml.load(f)
    Controls = file_yaml['Controls']
    for con in Controls:
        vcf_list.append(Path(con))
    Patients = file_yaml['Patients']
    for pat in Patients:
        vcf_list.append(Path(pat))

    Con_count = len(Controls)
    Pat_count = len(Patients)
    print("There are {} control samples and {} patient samples".format(Con_count, Pat_count),flush=True)


# identified the union snp location.
print("[", str(datetime.now()), "]", "starting building union snp location table.")
location = set()
for n, vcf_single in enumerate(vcf_list):
    with vcf_single.open() as vcf:
        for num, line in enumerate(vcf):
            if line.startswith("#") != True:
                skipline_num = num
                break

    with vcf_single.open() as vcf:
        for num,  line in enumerate(vcf):
            if num >= skipline_num:
                block = line.split('\t')
                chromosome = block[0]
                position = block[1]
                ref = block[3]
                PK = chromosome + ":" + position + ":" + ref
                location.add(PK)

Location = []
for i in range(0,len(location)):
    Location.append(location.pop())

Location = sorted(Location)


## transform the GT in vcf file into Genotype.
#from IPython.html.widgets import FloatProgress
#from IPython.display import display
#from time import sleep
starttime = datetime.now()
Genotype_all = pd.DataFrame()


for n, vcf_single in enumerate(vcf_list):
    print("[" + str(datetime.now()) + "]" + "processing " + str(n+1) +"/" + str(len(vcf_list)) + " sample: " + vcf_single.name, flush=True)
    with vcf_single.open() as vcf:
        PK_sample = []

        # define the header length of .vcf, prepare to skip it.
        for num, line in enumerate(vcf):
            if line.startswith("#") != True:
                skipline_num = num
                break


    # build the snp table for each sample first
    with vcf_single.open() as vcf:
        PK_sample = dict()
        Sample = []
        line_count = 0

        for num, line in enumerate(vcf):
            if num >= skipline_num:
                block = line.split('\t')
                chromosome = block[0]
                position = block[1]
                ref = block[3]
                PK = chromosome + ":" + position + ":" +ref
                PK_sample[PK] = line_count
                line_count += 1
                Sample.append(line)


        #f = FloatProgress(min = 0, max = len(Location)-1)
        #display(f)
        Genotype = ["None"]*len(Location)
        for i in range(0,len(Location)):
            #f.value = i
            row = Location[i]
            if row in PK_sample.keys():
                index = PK_sample[row]
                block = Sample[index].split('\t')
                ref = block[3]
                atl = block[4].split(',')
                if atl != '<CNV>':
                    table = [ref] + [a for a in atl]
                    GT = block[9].split(':')[0]
                    genotype_row = []
                    for geno_num in GT.split('/'):
                        genotype_row.append(table[int(geno_num)])

                    Genotype[i] = ''.join(genotype_row)
            else:
                ref = row.split(':')[2]
                Genotype[i] = ""+ref+ref
        #f.close()
        Genotype_all[vcf_single.name] = Genotype
        linecache.clearcache()


Genotype_all.index = Location

# starting the fisher exact test.
import scipy.stats as stats

controls = Genotype_all.columns[0:Con_count]
patients = Genotype_all.columns[Con_count: (Con_count + Pat_count)]
refine_table = pd.DataFrame(index=["Control_AA","Control_AB","Patient_AA","Patient_AB", "pvalue","oddsratio"])

#f = FloatProgress(min = 0, max = len(Genotype_all)-1)
#display(f)
starttime = datetime.now()
widgets = ['Fisher-exact-test: ', Percentage(), ' ',
           Bar(marker=RotatingMarker()), ' ', ETA(), ' ', FileTransferSpeed()]

pbar = ProgressBar(widgets=widgets)

print("[" + str(starttime) + "]" + "Now calculating the fisher exact test..." )
for i in pbar(range(0,len(Genotype_all))):
    #f.value = i
    ref = Genotype_all.index[i].split(":")[2]
    ref = ref+ref

    Con_AA = 0
    Con_AB = 0

    for con in controls:
        if Genotype_all[con][i] == ref:
            Con_AA += 1
        else:
            Con_AB += 1


    Pat_AA = 0
    Pat_AB = 0
    for pat in patients:
        if Genotype_all[pat][i] == ref:
            Pat_AA += 1
        else:
            Pat_AB += 1

    oddsratio, pvalue = stats.fisher_exact([[Con_AA, Pat_AA],[Con_AB,Pat_AB]])

    refine_table[Genotype_all.index[i]] = [Con_AA, Con_AB, Pat_AA, Pat_AB, pvalue, oddsratio]

pbar.finish()
#f.close()

print("[" + str(datetime.now()) + "]" + "Done, spend " + str(datetime.now() - starttime))

refine_table.T.to_csv(str(p.absolute()) + "/Fisher_test_output.csv")
