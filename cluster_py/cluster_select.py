# select top 6.25 from every cluster

import pandas as pd
import numpy as np
import os
import math


def selc_top(file_path, write_path):
    df = pd.read_csv(file_path, dtype={'cluster label': np.int32, 'num count': np.int32})
    # print(df.head())
    label = df.iloc[:, 3]
    dedup_label = list(set(label.tolist()))

    ii = 0
    df_select = pd.DataFrame()
    for i in dedup_label:
        # selc_name = 'df_' + str(i)
        selc_name = df[label == i].sort_values('energy', ascending=True)
        count = math.ceil(len(selc_name)*0.06)

        # if count < 1:
        #     count = 1
        #     df_i_new = df.iloc[:1, :]
        #
        # else:
        #     df_i_new = df.iloc[: count, :]

        selc_name = selc_name.iloc[: count, :]

        if ii == 0:
            df_select = selc_name

        else:
            df_select = pd.concat([df_select, selc_name], axis=0)
        ii += 1

    print(df_select.shape)
    df_select.to_csv(write_path, index=False)
    print('**'*20)


if __name__ == '__main__':
    os.chdir('/data/myproject/qsar_validation/fak/gen/single_r2_r3dock')
    file_path = 'qsar_r2_top5w_second_clustered_sele4w.csv'
    write_path = 'cluster_sele3k.csv'
    selc_top(file_path, write_path)
