from rdkit import Chem
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from voxels import VoxelMol


def read_mol2(mol2file):
    mol2bloks = []
    all_names = []
    with open(mol2file, 'r') as f:
        s = ""
        while True:
            line = f.readline()
            if not line:
                mol2bloks.append(s)
                break
            elif line[0] == "#":
                pass
            elif (line == "@<TRIPOS>MOLECULE\n") and (s != ""):
                mol2bloks.append(s)
                # print("mol2bloks: ",mol2bloks)
                s = ""
                f.seek(f.tell() - len(line))
            else:
                s = s + line
                if s.endswith("@<TRIPOS>MOLECULE\n"):
                    name = f.readline()
                    # print("name: ",name)
                    all_names.append(name.strip())
                    f.seek(f.tell() - len(name))
    return mol2bloks, all_names


def do_voxels(mol2_list, name_list, k=20):
    voxels_list = []
    for mol2, name_ in zip(mol2_list, name_list):
        mol = Chem.MolFromMol2Block(mol2, sanitize=True)
        vmols = VoxelMol(mol)
        voxels_list.append(vmols.voxels.flatten().tolist())
    voxels_array = np.array(voxels_list)

    kmeans_model = KMeans(n_clusters=k, random_state=42).fit(voxels_array)
    return kmeans_model.labels_


def main():
    mol2_path = '151.mol2'
    mol2s, names = read_mol2(mol2_path)
    label_ = do_voxels(mol2s, names)
    smis_list = [Chem.MolToSmiles(Chem.MolFromMol2Block(i, sanitize=True)) for i in mol2s if Chem.MolFromMol2Block(i, sanitize=True)]

    # w = Chem.SDWriter('output_1.sdf')
    # for one, two, three in zip(mol2s, names, label_.tolist()):
    #     mol = Chem.MolFromMol2Block(one, sanitize=True)
    #     print(type(mol))
    #     # mol.SetProp('ID', two)
    #     mol.SetProp('shape_class_', three)
    #     w.write(mol)
    # w.close()
    pd.DataFrame({'Smiles': smis_list, 'ID': names, 'shape_class': label_.tolist()}).to_csv('tmp.csv', index=False)


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    end_time = time.time()
    print('Shape cluster cost time: {:.3f}'.format((end_time-start_time)/60.0))


