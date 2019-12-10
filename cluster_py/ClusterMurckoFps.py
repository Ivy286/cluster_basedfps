from cluster_m import *
from murcko_scaffold import *
import pandas as pd
import optparse
import sys
import os

##########################################
## Options and defaults
##########################################

# DATAPATH = os.getenv('DATAPATH')
# SAVEDPATH = os.getenv('SAVEDPATH')


def getOptions():
    parser = optparse.OptionParser('python *.py [option]')
    parser.add_option('--input1', dest='input1', help='format: csv', default='')
    parser.add_option('--input2', dest='input2', help='format: csv', default='')
    parser.add_option('--smi_col', dest='smi_col', help='if file is csv, please give smiles columns, default is 0')
    parser.add_option('--id_col', dest='id_col', help='if file is csv, please give ID columns, default is 1')
    parser.add_option('--fp', dest='fp', help='fingerprint type: tp,mc,mo (Topological Fingerprints, MACCS Keys, Morgan Fingerprints), default is mc', default='mc')
    parser.add_option('--radius', dest='radius', help=' the radius of the Morgan fingerprint, default is 2',type='int', default=2)
    parser.add_option('--algorithm', dest='algorithm', help='cluster algorithm :b,m (Butina, Murtagh), default is b', default='b')
    parser.add_option('--cutoff', dest='cutoff', help='distThresh(0-1),elements within this range of each other are considered to be neighbors, needed for Butina cluster algorithm, default is 0.5', type='float', default=0.5)
    parser.add_option('--nclusts', dest='nclusts', help='number of clusters, needed for Murtagh cluster algorithm, default is 1',type='int', default=1)
    parser.add_option('--murtype', dest='Murtype', help='Method for Murtagh:Wards, SLINK, CLINK, UPGMA, needed when Murtagh is set as algorithm, default is Wards', default='Wards')
    parser.add_option('--out', dest='output', help='output sdf file or csv file', default='')
    options, args = parser.parse_args()

    print(options.input1, options.input2, options.output)
    # options.input = DATAPATH
    # options.output = os.path.join(SAVEDPATH, options.output)
    # if options.input1 == '' or options.input2 or options.output == '':
    #     parser.print_help()
    #     print("No input1 or input2 or output is provided")
    #     sys.exit(1)

    return options


def main():
    options = getOptions()
    fpOpdict = {'tp': 'Topological Fingerprints', 'mc': 'MACCS Keys', 'mo': 'Morgan Fingerprints'}
    algOpdict = {'b': 'Butina', 'm': 'Murtagh'}
    options.algorithm = algOpdict[options.algorithm]
    print("fingerprint type: %s" % fpOpdict[options.fp])
    if options.fp == 'mo':
        print("radius: %s" % str(options.radius))
    print("cluster algorithm: %s" % options.algorithm)
    if options.algorithm == "Murtagh":
        print("Murtagh method: %s" % options.Murtype)
        print("Murtagh cluster number set: %s" % options.nclusts)
    elif options.algorithm == "Butina":
        print("cutoff(distThresh) : %s" % options.cutoff)

    print('input file reading...')
    smis_dataprocess(merge_file(options.input1, options.input2), smi_col=0, id_col=1)

    merge_df = pd.read_csv('clustered_based_scaffold.csv')
    merge_df_new = merge_df[merge_df['num_count'] >= 5]

    df_tmp = pd.DataFrame()
    for num, i in enumerate(set(merge_df_new['cluster_label'].tolist())):
        df_item = merge_df[merge_df['cluster_label'] == int(i)]
        input_parse = ChemParse(df_item)
        print('cluster%d fingerprint calculating...' % i)
        input_parse.get_fps(options.fp, options.radius)
        print('cluster%d clustering...' % i)
        fp_cluster = Fingerprint_Cluster(input_parse.fps)
        fp_cluster.distance_matrix()
        fp_cluster.cluster_dict(options.algorithm, options.cutoff, options.Murtype, options.nclusts)
        print('done, output of cluster{} to {}'.format(i, options.output))
        if num == 0:
            out_df, index = input_parse.clusterOutput(options.output, fp_cluster.cdict)
            df_tmp = out_df
        else:
            out_df, index = input_parse.clusterOutput(options.output, fp_cluster.cdict, index_=int(len(fp_cluster.clustdict)))
            df_tmp = pd.concat([df_tmp, out_df], axis=0)

        df_final = pd.concat([df_tmp, merge_df[merge_df['num_count'] < 5]], axis=0)
        df_final = pd.merge(df_final, merge_df, on='ID')  #.drop(['drug_id', ''],axis=1)
        # for c in fp_cluster.clustdict:
        #     print("cluster%s: %s" % (str(c), str(fp_cluster.clustdict[c])))
        df_final.to_csv(options.output, index=False)


if __name__ == "__main__":
    # os.chdir('./')
    # df = merge_file('r2_select_top1w.csv', 'r2_select_top1w_docking_result_processed.csv')
    # smis_dataprocess(df, smi_col=0, id_col=1)
    main()
