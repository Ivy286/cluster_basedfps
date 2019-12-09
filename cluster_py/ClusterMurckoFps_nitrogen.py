from ClusterFps import *
import optparse
import sys
import os

##########################################
## Options and defaults
##########################################

DATAPATH = os.getenv('DATAPATH')
SAVEDPATH = os.getenv('SAVEDPATH')


def getOptions():
    parser = optparse.OptionParser('python *.py [option]')
    parser.add_option('--in', dest='input', help='intput sdf, csv or smi file', default='')
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

    options.input = DATAPATH
    options.output = os.path.join(SAVEDPATH, options.output)
    if options.input == '' or options.output == '':
        parser.print_help()
        print("No input or output is provided")
        sys.exit(1)

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
    input_parse = ChemParse(options.input)
    if vars(options)['smi_col'] is None or vars(options)['id_col'] is None:
        # print('888888888888888888888888888888888')
        input_parse.input_reader()
    else:
        input_parse.input_reader(vars(options)['smi_col'], vars(options)['id_col'])
    print('fingerprint calculating...')
    input_parse.get_fps(options.fp, options.radius)
    print('clustering...')
    fp_cluster = Fingerprint_Cluster(input_parse.fps)
    fp_cluster.distance_matrix()
    fp_cluster.cluster_dict(options.algorithm, options.cutoff, options.Murtype, options.nclusts)
    print('done, output to %s' % options.output)
    input_parse.clusterOutput(options.output, fp_cluster.cdict)
    for c in fp_cluster.clustdict:
        print("cluster%s: %s" % (str(c), str(fp_cluster.clustdict[c])))
    
    
if __name__ == "__main__":
    main()
