# coding: utf-8
# %load ./vcf_genotype_fisher.py
from pathlib import Path
import os
import sys
import argparse
import yaml
import pandas as pd
from datetime import datetime


# argument parse
parser = argparse.ArgumentParser(description="Input the list of control and patient vcf file in yaml type.")
parser.add_argument('-l','--list', type=str, help="the yaml file must put with all the vcf files in the same folder.", required=True, metavar=".yaml")
args = parser.parse_args()
list_path = args.list
p = Path(list_path).parent

print(p.absolute())
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
    print("gathering from ", vcf_single.name)
    with vcf_single.open() as vcf:
        for num, line in enumerate(vcf):
            if line.startswith("#") != True:
                skipline_num = num
                break

    with vcf_single.open() as vcf:
        for num,  line in enumerate(vcf):
            if num >= skipline_num:
                try:
                    block = line.split('\t')
                    chromosome = block[0]
                    position = block[1]
                    ref = block[3]
                    alt = block[4]
                    PK = chromosome + ":" + position + ":" + ref + ":" + alt
                    location.add(PK)
                except:
                    # print("Warning: there are '#' outside off header")
                    pass

Location=[]
                
for i in range(0,len(location)):
    Location.append(location.pop())

Location = sorted(Location)
print("the union snp location table contains ", len(Location), " row.")


chromo = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 
          'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
          'chr20','chr21','chr22','chrX','chrY','chrM']

if not os.path.exists(str(p.absolute()) + '/Snp_index'):
    os.mkdir(str(p.absolute()) + '/Snp_index')
else:
    pass
    print("The SNP_index folder exists. Rewriting.", flush=True)
    
    
    

for c in chromo:
    with open(str(p.absolute()) + '/Snp_index/' + c, 'w') as file:
        for line in Location:
            if line.startswith(c):
                print(line, file=file)

del(Location)
#----------------------------------
## transform the GT in vcf file into Genotype.
#from IPython.html.widgets import FloatProgress
#from IPython.display import display
#from time import sleep
starttime = datetime.now()

for chromo_file in os.listdir(str(p.absolute()) + '/Snp_index'):
    print("[" + str(datetime.now()) + "]" + "Handling", chromo_file)
    chromo_index = []
    with open(str(p.absolute()) + '/Snp_index/' + chromo_file) as f:
        for line in f:
            chromo_index.append(line)
            
    if len(chromo_index) != 0:        
        Genotype_all_chromo = pd.DataFrame()
        for n, vcf_single in enumerate(vcf_list):
            print("\tprocessing " + str(n+1) +"/" + str(len(vcf_list)) + " sample: " + vcf_single.name, flush=True)
            with vcf_single.open() as vcf:

                # define the header length of .vcf, prepare to skip it.
                for num, line in enumerate(vcf):
                    if line.startswith("#") != True:
                        skipline_num = num
                        break

            # build the snp table for each sample first
            with vcf_single.open() as vcf:
                PK_sample = dict()  # dict(vcf_per_line:index)
                Sample = []
                line_count = 0

                for num, line in enumerate(vcf):
                    if num >= skipline_num:
                        try:
                            block = line.split('\t')
                            chromosome = block[0]
                            position = block[1]
                            ref = block[3]
                            alt = block[4]
                            PK = chromosome + ":" + position + ":" +ref + ":" + alt
                            PK_sample[PK] = line_count
                            line_count += 1
                            Sample.append(line)
                        except:
                            # print("Warning: there are '#' outside off header")
                            pass

                #f = FloatProgress(min = 0, max = len(Location)-1)
                #display(f)

                Genotype = ["None"]*len(chromo_index)
                for i in range(0,len(chromo_index)):
                    #f.value = i
                    row = chromo_index[i]
                    if row in PK_sample.keys():   # get the vcf line of sample, and transfer the GT label.
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
                Genotype_all_chromo[vcf_single.name] = Genotype

        Genotype_all_chromo.index = chromo_index

        if not os.path.exists(str(p.absolute()) + '/Genotype_table'):
            os.mkdir(str(p.absolute()) + '/Genotype_table')

        Genotype_all_chromo.to_csv(str(p.absolute()) + '/Genotype_table/' + chromo_file + '.csv')
    else:
        print(chromo_file, " had no snp on it, skip it")
        pass
        
