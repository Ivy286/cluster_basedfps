# !-*- coding: utf-8 -*-
# 骨架聚类后，从每类骨架中按照PIC50随机分train和test，并画PIC50分布图

import pandas as pd
import numpy as np
import os
import math
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt


def selc_top(file_path, write_path_train, write_path_test):
    df = pd.read_csv(file_path, dtype={'cluster_label': np.int32, 'num_count': np.int32})
    print(df[df.cluster_label == 0].shape)
    label = df.iloc[:, 2]
    # print(df.shape)
    dedup_label = list(set(label.tolist()))
    # print(dedup_label)

    ii = 0
    df_select = pd.DataFrame()
    df_test = pd.DataFrame()
    for i in dedup_label:
        selc_name = df[df['cluster_label'] == i]
        # selc_name = df[label == i]  # 此处label是数据流df.iloc[:, 2]，等于df['cluster_label']

        y_test, y_train = train_test_split(selc_name, train_size=0.18)
        print(y_test.shape, y_train.shape)
        if ii == 0:
            df_select = y_train
            df_test = y_test

        else:
            df_select = pd.concat([df_select, y_train], axis=0)
            df_test = pd.concat([df_test, y_test], axis=0)
        ii += 1

    df_select.to_csv(write_path_train, index=False)
    df_test.to_csv(write_path_test, index=False)
    print('**'*20)
    print(df_select.shape)


if __name__ == '__main__':
    os.chdir('/data/myproject/qsar_validation/fak/qsar_model')
    file_path = '286.csv'

    # write_path_train = 'train.csv'
    # write_path_test = 'test.csv'
    write_path_train = 'defactinib_288_train_mpfps.csv'
    write_path_test = 'defactinib_288_test_del2.csv'
    # selc_top(file_path, write_path_train, write_path_test)

    plt.figure(figsize=(6, 18))
    fontdict = {'family': 'Times New Roman', 'size': 9}
    df1 = pd.read_csv(file_path)
    df2 = pd.read_csv(write_path_train)
    df3 = pd.read_csv(write_path_test)
    y1 = df1['PIC50']

    ax1 = plt.subplot(3, 1, 1)
    n, bins, patches = ax1.hist(y1, bins=10, label='or', normed=True, color='green', alpha=0.7)
    ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax1.set_xticks([5, 6, 7, 8, 9])
    ax1.set_xlim([4, 10])
    plt.xticks(fontproperties='Times New Roman', size=8)
    plt.yticks(fontproperties='Times New Roman', size=8)
    plt.xlabel('PIC50', fontdict)
    plt.ylabel('Frequence', fontdict)
    print(bins)

    y2 = df2['PIC50']
    ax2 = plt.subplot(3, 1, 2)
    ax2.hist(y2, bins=bins, normed=True, alpha=0.7)
    ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax2.set_xticks([5, 6, 7, 8, 9])
    ax2.set_xlim([4, 10])
    plt.xticks(fontproperties='Times New Roman', size=8)
    plt.yticks(fontproperties='Times New Roman', size=8)
    plt.xlabel('PIC50', fontdict)
    plt.ylabel('Frequence', fontdict)

    y3 = df3['PIC50']
    ax3 = plt.subplot(3, 1, 3)
    ax3.hist(y3, bins=bins, normed=True, alpha=0.7)
    ax3.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax3.set_xticks([5, 6, 7, 8, 9])
    ax3.set_xlim([4, 10])
    plt.xticks(fontproperties='Times New Roman', size=8)
    plt.yticks(fontproperties='Times New Roman', size=8)
    plt.xlabel('PIC50', fontdict)
    plt.ylabel('Frequence', fontdict)
    print(len(y1), len(y2), len(y3))

    plt.show()
