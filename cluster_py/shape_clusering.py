# coding: utf-8

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors3D
import numpy as np
import pandas as pd
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
    # print(kmeans_model.cluster_centers_)
    # print(kmeans_model.cluster_centers_.shape)
    # print(pd.Series(kmeans_model.labels_).value_counts())
    # print('===============')

    # centers = kmeans_model.cluster_centers_
    # labelss = kmeans_model.labels_
    # print(labelss)
    # y_pred = kmeans_model.fit_predict(X)
    #
    # colors = ['#4EACC5', '#FF9C34', '#4E9A06', 'blue', 'green', 'red', 'black']
    # import matplotlib.pyplot as plt
    # plt.figure()
    # for i in range(len(labelss)):
    #     index_sets = np.where(y_pred == i)
    #     cluster = X[index_sets]
    #     plt.scatter(cluster[:, 0], cluster[:, 1], c=colors[i], marker='.')
    #     plt.plot(centers[i][0], centers[i][1], 'o', markerfacecolor=colors[i], markeredgecolor='k', markersize=6)
    # plt.show()

    return kmeans_model.labels_


if __name__ == '__main__':
    mols = Chem.SDMolSupplier('151.sdf', removeHs=False)
    
    k = 50
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

