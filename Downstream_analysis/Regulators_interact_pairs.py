import pandas as pd
from itertools import product

def Both_signif_pairs():
    exp_signif_list = pd.read_csv('~/Downstream_analysis/Regulators/Lung.m6A_peak_RBP.m6A2me.signif.list', sep = '\t', header = None)
    out_signif_list = pd.read_csv('~/Downstream_analysis/Regulators/Lung.DNAme_TF.m6A2me.signif.list', sep = '\t', header = None)
    d_pair = pd.read_csv('~/Downstream_analysis/Regulators/Lung.m6A2me_signif.Proten_pairs_sites.txt', sep = '\t', header = 0)
    d_string = pd.read_csv('~/Downstream_analysis/Regulators/Lung.m6A2me.Proten_pairs.signif_String.tsv', sep = '\t', header = 0)

    #Protein order from exp to out
    signif_pair_list = list(product(exp_signif_list[0].to_list(), out_signif_list[0].to_list()))
    d_pair['protein_combine'] = d_pair.apply(lambda x: (x['RBP'], x['TF']), axis = 1)
    #No order
    d_string['string_combine'] = d_string.apply(lambda x: (x['node1'], x['node2']), axis = 1)
    interact_pair_list = d_string['string_combine'].to_list()
    experi_pair_list = d_string.loc[d_string['experimentally_determined_interaction'] > 0, 'string_combine'].to_list()

    res0 = pd.DataFrame(data = None, columns = ['Protein_pair', 'count', 'signif_or_not', 'interact_score', 'experi_score'])
    idx0 = 0
    for idx, sd in d_pair.groupby('protein_combine'):
        Protein_pair = '_'.join([str(i) for i in idx])
        count = sd.shape[0]
        if idx in signif_pair_list:
            signif_or_not = 1
            if idx in interact_pair_list or tuple(reversed(idx)) in interact_pair_list:
                if idx in interact_pair_list:
                    interact_score = d_string[d_string['string_combine'] == idx].reset_index()['combined_score'][0]
                else:
                    interact_score = d_string[d_string['string_combine'] == tuple(reversed(idx))].reset_index()['combined_score'][0]
                if idx in experi_pair_list or tuple(reversed(idx)) in experi_pair_list:
                    if idx in experi_pair_list:
                        experi_score = d_string[d_string['string_combine'] == idx].reset_index()['experimentally_determined_interaction'][0]
                    else:
                        experi_score = d_string[d_string['string_combine'] == tuple(reversed(idx))].reset_index()['experimentally_determined_interaction'][0]
                else:
                    experi_score = 0
            else:
                interact_score = 0
                experi_score = 0
        else:
            signif_or_not = 0
            interact_score = 0
            experi_score = 0

        res0.loc[idx0] = [Protein_pair, count, signif_or_not, interact_score, experi_score]
        idx0 += 1

        res0.to_csv('~/Downstream_analysis/Regulators/Lung.m6A2me.Proten_pairs.signif_String_infos.txt', sep = '\t', index = None)

def Single_side_regulators():
    Protein_signif_list = pd.read_csv('~/Downstream_analysis/Regulators/Lung.m6A_peak_RBP.m6A2me.signif.list', sep = '\t', header = None)
    Regulator_list = pd.read_csv('~/Downstream_analysis/Regulators/DNAme_regulators.list', sep = '\t', header = None)
    Regulator_list = Regulator_list[~Regulator_list[0].str.startswith('#')]
    d_string = pd.read_csv('~/Downstream_analysis/Regulators/Lung.m6A2me.Proten_single_Regulators.signif_String.tsv', sep = '\t', header = 0)

    Protein_Regulator_list = list(product(Protein_signif_list[0].to_list(), Regulator_list[0].to_list()))
    d_string['string_combine'] = d_string.apply(lambda x: (x['node1'], x['node2']), axis = 1)
    interact_list = d_string['string_combine'].to_list()
    experi_list = d_string.loc[d_string['experimentally_determined_interaction'] > 0, 'string_combine'].to_list()

    res0 = pd.DataFrame(data = None, columns = ['Protein_Regulator', 'interact_score', 'experi_score'])
    idx0 = 0
    for i in Protein_Regulator_list:
        Protein_Regulator = '_'.join([str(i0) for i0 in i])
        if i in interact_list or tuple(reversed(i)) in interact_list:
            if i in interact_list:
                interact_score = d_string[d_string['string_combine'] == i].reset_index()['combined_score'][0]
            else:
                interact_score = d_string[d_string['string_combine'] == tuple(reversed(i))].reset_index()['combined_score'][0]
            if i in experi_list or tuple(reversed(i)) in experi_list:
                if i in experi_list:
                    experi_score = d_string[d_string['string_combine'] == i].reset_index()['experimentally_determined_interaction'][0]
                else:
                    experi_score = d_string[d_string['string_combine'] == tuple(reversed(i))].reset_index()['experimentally_determined_interaction'][0]
            else:
                experi_score = 0
        else:
            interact_score = 0
            experi_score = 0
         
        res0.loc[idx0] = [Protein_Regulator, interact_score, experi_score]
        idx0 += 1

    res0 = res0[res0['interact_score'] > 0]
    res0.to_csv('~/Downstream_analysis/Regulators/Lung.m6A2me.Proten_single_Regulators.txt', sep = '\t', index = None)

Both_signif_pairs()
#Single_side_regulators()
