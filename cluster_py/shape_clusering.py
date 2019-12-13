# coding: utf-8

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors3D
import numpy as np
from sklearn.cluster import KMeans


def shape_clustering(mols, k=20):
    mols_features = []
    for mol in mols:
        if mol is None:
            continue
        descrptors = [
            Descriptors3D.NPR1(mol),
            Descriptors3D.NPR2(mol),
            Descriptors3D.Asphericity(mol),
            Descriptors3D.Eccentricity(mol),
            Descriptors3D.InertialShapeFactor(mol),
            # Descriptors3D.PMI1(mol),
            # Descriptors3D.PMI2(mol),
            # Descriptors3D.PMI3(mol),
            Descriptors3D.RadiusOfGyration(mol),
            Descriptors3D.SpherocityIndex(mol),
            rdMolDescriptors.CalcPBF(mol)]
        mols_features.append(descrptors)
    X = np.array(mols_features)
    
    kmeans_model = KMeans(n_clusters=k, random_state=42).fit(X)
    return kmeans_model.labels_


if __name__ == '__main__':
    mols = Chem.SDMolSupplier('test.sdf', removeHs=False)
    
    k = 2
    labels = shape_clustering(mols, k)

    w = Chem.SDWriter('output.sdf')
    idx = 0
    for mol in mols:
        if mol is None:
            continue
        mol.SetProp('shape_category_of_' + str(k), str(labels[idx]))
        w.write(mol)
        idx += 1
    w.close()

