# Cluster based on 3D shape

from rdkit import Chem
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn import metrics
from voxels import VoxelMol

from pyGPGO.covfunc import matern32
from pyGPGO.acquisition import Acquisition
from pyGPGO.surrogates.GaussianProcess import GaussianProcess
from pyGPGO.GPGO import GPGO


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


def do_voxels(mol2_list, name_list):
    voxels_list = []
    for mol2, name_ in zip(mol2_list, name_list):
        mol = Chem.MolFromMol2Block(mol2, sanitize=True)
        vmols = VoxelMol(mol)
        rdkitmol = vmols.mol
        vmol = VoxelMol(rdkitmol)
        voxels_list.append(vmol.voxels.flatten().tolist())
    voxels_array = np.array(voxels_list)
    # print(voxels_array.shape)
    # print('============')

    def cluster_shape(n_clusters):

        y_pred = KMeans(n_clusters=int(n_clusters), random_state=0).fit_predict(voxels_array)
        # print(y_pred.shape)
        score = metrics.calinski_harabaz_score(voxels_array, y_pred)
        return score

    cov = matern32()
    gp = GaussianProcess(cov)
    acq = Acquisition(mode='ExpectedImprovement')
    param = {
        'n_clusters': ('int', [3, 5]),
    }

    random_state = 0
    gp_max_iter = 3
    np.random.seed(random_state)
    gpgo = GPGO(gp, acq, cluster_shape, param)
    gpgo.run(max_iter=gp_max_iter, init_evals=2)
    para_result = gpgo.getResult()[0]['n_clusters']
    print('Shape cluster: {}'.format(para_result))
    print('*'*60)

    kmeans_model = KMeans(n_clusters=int(para_result), random_state=0).fit(voxels_array)
    return kmeans_model.labels_


def main():
    mol2_path = 'test.mol2'
    mol2s, names = read_mol2(mol2_path)

    new_mol2s = []
    new_names = []
    for ii, jj in zip(mol2s, names):
        if Chem.MolFromMol2Block(ii):
            new_mol2s.append(ii)
            new_names.append(jj)
    label_ = do_voxels(new_mol2s, new_names)

    # write csv file
    smis_list = [Chem.MolToSmiles(Chem.MolFromMol2Block(i, sanitize=True)) for i in mol2s if Chem.MolFromMol2Block(i, sanitize=True)]
    pd.DataFrame({'Smiles': smis_list, 'ID': new_names, 'shape_class': label_.tolist()}).to_csv('tmp_test.csv', index=False)
    w = Chem.SDWriter('output_1.sdf')

    # write sdf file
    for one, two, three in zip(new_mol2s, new_names, label_.tolist()):
        tmp_mol = Chem.MolFromMol2Block(one)
        tmp_mol.SetProp('shape_class', str(three))
        tmp_mol.SetProp('_Name', str(two))
        w.write(tmp_mol)
    w.close()


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    end_time = time.time()
    print('Shape cluster cost time: {:.3f} min'.format((end_time-start_time)/60.0))


