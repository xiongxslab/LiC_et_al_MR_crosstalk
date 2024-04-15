import sys
import pandas as pd

def Pair_mr():
    d1 = pd.read_csv(qtl1, sep='\t', header=None, usecols=[0,1], names=['Phenotype', 'SNP']).groupby('Phenotype')
    d2 = pd.read_csv(qtl2, sep='\t', header=None, usecols=[0,1], names=['Phenotype', 'SNP']).groupby('Phenotype')
    dict1, dict2 = {}, {}
    for Phe1, s1 in d1:
        dict1[Phe1] = s1['SNP'].to_list()
    for Phe2, s2 in d2:
        dict2[Phe2] = s2['SNP'].to_list()
    with open(fout, 'w+') as f:
        for Phe1, s1 in dict1.items():
            for Phe2, s2 in dict2.items():
                if len(set(s1) & set(s2)) != 0:
                    f.write(Phe1 + '\t' + Phe2 + '\n')
tis = sys.argv[1]
Phe1 = sys.argv[2]
Phe2 = sys.argv[3]
nchr = sys.argv[4]

if Phe1 == "m6A":
    Phe1_s = "m6A-"
else:
    Phe1_s = Phe1

if Phe2 == "m6A":
    Phe2_s = "m6A-"
else:
    Phe2_s = Phe2

qtl1 = '~/Database/' + tis + '_' + Phe1_s + 'QTL_signif/' + tis + '.' + Phe1_s + 'QTL.signif_chr' + nchr + '.txt'
qtl2 = '~/Database/' + tis + '_' + Phe2_s + 'QTL_all/' + tis + '.' + Phe2_s + 'QTL.all_chr' + nchr + '.txt'
fout = '~/2SMR/' + Phe1 + '2' + Phe2 + '/' + tis + '_' + Phe1 + '2' + Phe2 + '_pairs/' + tis + '.' + Phe1 + '2' + Phe2 + '_chr' + nchr + '.txt'

Pair_mr()
