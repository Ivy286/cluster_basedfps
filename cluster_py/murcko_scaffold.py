from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric, GetScaffoldForMol
import pandas as pd
import os


def merge_file(file_input1, file_input2):
    df1 = pd.read_csv(file_input1)
    df2 = pd.read_csv(file_input2)
    df_merge = pd.merge(df1, df2, left_on='drugs', right_on='drug_id').drop(['drugs', 'target_pdb_id', 'targets'], axis=1)
    df_merge = df_merge[['drug_smiles', 'drug_id', 'pdb_id', 'Energy', 'pose_file_name']]
    return df_merge


def cluster_by_Murcko_scaffold(smiles_list, id_list, use_generics=True, remove_linker_side_chain=True):
    scaffold2smiles_dict = {}
    count = 0
    total_count = len(smiles_list)
    for smi_item, id_item in zip(smiles_list, id_list):
        mol = Chem.MolFromSmiles(smi_item)
        count += 1
        try:
            scaffold = GetScaffoldForMol(mol)
            if use_generics:
                scaffold = MakeScaffoldGeneric(scaffold)
                # note that for atom that have more than 5 bonds(like phosphorus)
                # this method would raise an error
            if remove_linker_side_chain:
                scaffold = GetScaffoldForMol(scaffold)
                if use_generics:
                    scaffold = MakeScaffoldGeneric(scaffold)
                    # note that for atom that have more than 5 bonds(like phosphorus)
                    # this method would raise an error

            scaffold_smiles = Chem.MolToSmiles(scaffold)
            if scaffold_smiles not in scaffold2smiles_dict:
                scaffold2smiles_dict[scaffold_smiles] = [smi_item]
            else:
                scaffold2smiles_dict[scaffold_smiles].append(smi_item)
        except Exception:
            print('sth is wrong when extracting scaffold from given smiles:', smi_item)
            continue

        if count % 500 == 0:
            print('finished clustering molecules %d / %d' % (count, total_count))

    smiles2scaffold_dict = {}
    for scaffold_ in scaffold2smiles_dict:
        for smi_item_ in scaffold2smiles_dict[scaffold_]:
            smiles2scaffold_dict[smi_item_] = scaffold_
    return scaffold2smiles_dict, smiles2scaffold_dict


def get_cluster_count(scaffold2smiles_dict, scaffold_list, use_str=True):
    count_list = []
    for smi_each in scaffold_list:
        count = len(scaffold2smiles_dict[smi_each])
        if use_str:
            count = str(count)
        count_list.append(count)
    return count_list


def smis_dataprocess(read_text, smi_col=0, id_col=1):

    # if input_file[-4:] == '.csv':
    #     read_text = pd.read_csv(input_file)
    # if input_file[-4:] == '.smi':
    #     read_text = pd.read_csv(input_file, sep='\t', header=None)
    smis_list = read_text.iloc[:, smi_col].tolist()
    id_list = read_text.iloc[:, id_col].tolist()

    # Get smiles to name dict
    smiles2name = {}
    for name, smi in zip(id_list, smis_list):
        smiles2name[smi] = name

    valid_smiles_list = []
    valid_idx_list = []

    for ii, smiles in zip(id_list, smis_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            valid_smiles_list.append(smiles)
            valid_idx_list.append(ii)
        else:
            print('Invalid smile. Index is:', ii)

    smiles_list = valid_smiles_list
    idx_list = valid_idx_list

    scaffold2smiles_dict, smiles2scaffold_dict = cluster_by_Murcko_scaffold(smis_list, idx_list)

    # Save clustered result to csv
    clustered_smiles_list = []
    clustered_scaffold_list = []
    clustered_names_list = []
    clustered_label_list = []
    for ii, (key, value) in enumerate(smiles2scaffold_dict.items()):
        # print('8'*100)
        # print(ii)
        # print(key)
        # print(value)
        clustered_smiles_list.append(key)
        clustered_scaffold_list.append(value)
        clustered_names_list.append(smiles2name[key])

    df_clustered = pd.DataFrame()
    df_clustered['Smiles'] = clustered_smiles_list
    df_clustered['clustered_scaffold'] = clustered_scaffold_list
    df_clustered['ID'] = clustered_names_list

    # Give each scaffold a number
    scaffold2category = {}
    scaffold_set = set(clustered_scaffold_list)
    for i, scaffold_unique in enumerate(scaffold_set):
        scaffold2category[scaffold_unique] = i
    # print(scaffold2category)

    for ii, (key, value) in enumerate(smiles2scaffold_dict.items()):
        clustered_label_list.append(scaffold2category[value])
    # print(clustered_label_list)

    df_clustered['cluster_label'] = clustered_label_list

    count_list = get_cluster_count(scaffold2smiles_dict, clustered_scaffold_list, use_str=True)

    print('=' * 80)
    print('Length of count_list', len(count_list))
    df_clustered['num_count'] = count_list
    print('Total number of clusters:', len(scaffold_set))
    df_clustered = pd.merge(df_clustered, read_text[['drug_id', 'Energy']], left_on='ID', right_on='drug_id').drop(['drug_id'],axis=1)

    df_clustered.to_csv('./clustered_based_scaffold.csv', index=False)
    return df_clustered


if __name__ == "__main__":

    os.chdir('./')
    df = merge_file('r2_select_top1w.csv', 'r2_select_top1w_docking_result_processed.csv')
    smis_dataprocess(df, smi_col=0, id_col=1)
