# read sdf, smi or csv file and do fingerprint calculation AND clustering
import pandas as pd
import numpy
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina, Murtagh
from rdkit.ML.Cluster import ClusterUtils


class ChemParse(object):
    def __init__(self, input_file):
        self.source = ''
        self.fps = ''
        self.input_file = input_file
        self.ids = ''
        
    def input_reader(self, smi_col='Smiles', id_col='ID'):
        print('input file format', self.input_file[-4:])
        # print('*******'*10)

        if self.input_file[-4:] == '.sdf':
            sdf_mols = []
            for sdf_item in Chem.SDMolSupplier(self.input_file):
                if sdf_item:
                    sdf_mols.append(sdf_item)
            self.source = sdf_mols
        elif self.input_file[-4:] == '.csv':
            csv_mols = []
            ids_new = []
            csv_df = pd.read_csv(self.input_file)
            for csv_item, id_item in zip(csv_df[smi_col], csv_df[id_col]):
                csv_mol = Chem.MolFromSmiles(csv_item)
                if csv_mol:
                    csv_mols.append(csv_mol)
                    ids_new.append(id_item)
            self.source = csv_mols
            self.ids = ids_new
        elif self.input_file[-4:] == '.smi':
            smi_mols = []
            smi_ids = []
            smi_df = pd.read_csv(self.input_file, header=None, sep='\t')
            for smi_item, id_item in zip(smi_df.iloc[:, 0], smi_df.iloc[:, 1]):
                smi_mol = Chem.MolFromSmiles(smi_item)
                if smi_mol:
                    smi_mols.append(smi_mol)
                    smi_ids.append(id_item)
            self.source = smi_mols
            self.ids = smi_ids

        else:
            print('please input correct file: sdf, csv or smi')
        print('valid smiles number: ', len(self.source))
        # print('======'*10)
        
    def get_fps(self, ftype, radius=1):
        fps = []
        if self.source == '' or self.source is None:
            print('self.source in None')
            pass
        else:
            if ftype == 'mo':
                fps = [AllChem.GetMorganFingerprint(m, radius, useFeatures=True) for m in self.source]
            elif ftype == 'tp':
                fps = [FingerprintMols.FingerprintMol(m) for m in self.source]
            elif ftype == 'mc':
                fps = [MACCSkeys.GenMACCSKeys(m) for m in self.source]
        self.fps = fps
        print('calculated fps number: ', len(self.fps))

    def clusterOutput(self, output, cdict):
        if output[-4:] == '.sdf':
            sdfout = Chem.SDWriter(output)

            for index, m in enumerate(self.source):
                classid = cdict[index][0]
                isCentroid = cdict[index][1]
                m.SetProp("fp_cluster", str(classid))
                m.SetProp("isCentroid", isCentroid)
                sdfout.write(m)
            sdfout.close()
        if output[-4:] == '.csv':
            class_id = []
            centroid_id = []
            smiles = []
            for index, m in enumerate(self.source):
                class_id.append(cdict[index][0])
                centroid_id.append(cdict[index][1])
                smiles.append(Chem.MolToSmiles(m))
            pd.DataFrame({'Smiles': smiles, 'fp_cluster': class_id, 'isCentroid': centroid_id}).to_csv(output, index=False)


class Fingerprint_Cluster(object):
    def __init__(self, fps):
        self.fplist = fps
        self.dist = []
        self.cdict = {} 
        self.clustdict = {}

    # generate the distance matrix
    def distance_matrix(self):
        self.dist = []
        nfps = len(self.fplist)
        for i in range(1, nfps):
            sims = DataStructs.BulkTanimotoSimilarity(self.fplist[i], self.fplist[:i])
            self.dist.extend([1-x for x in sims])
        
    # generate cluster dict as {1:[1,2,3],2:[4,5,6]...}
    def cluster_dict(self, algorithm, cutoff=0.5, method='Wards', ncluster=1):
        if algorithm == 'Butina':
            self.ClusterFps_Butina(self.dist, len(self.fplist), cutoff)
        elif algorithm == 'Murtagh':
            self.ClusterFps_Murtagh(self.dist, len(self.fplist), method, ncluster)
        
    def ClusterFps_Butina(self, dists, nfps, cutoff):
        self.cdict = {}
        cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
        for index, eachcs in enumerate(cs): 
            self.clustdict[index+1] = eachcs
            for eachid in eachcs:
                self.cdict[eachid] = [index+1]
                if eachid == eachcs[0]:
                    self.cdict[eachid].append("true")
                else:
                    self.cdict[eachid].append("flase")
            
    def ClusterFps_Murtagh(self, dists, nfps, method, ncluster):
        self.cdict = {}
        cs = None
        if method == 'Wards':
            cs = Murtagh.ClusterData(dists, len(self.fplist), Murtagh.WARDS, isDistData=1)
        elif method == 'SLINK':
            cs = Murtagh.ClusterData(dists, len(self.fplist), Murtagh.SLINK, isDistData=1)
        elif method == 'CLINK':
            cs = Murtagh.ClusterData(dists, len(self.fplist), Murtagh.CLINK, isDistData=1)
        elif method == 'UPGMA':
            cs = Murtagh.ClusterData(dists, len(self.fplist), Murtagh.UPGMA, isDistData=1)

        splitClusts = ClusterUtils.SplitIntoNClusters(cs[0], ncluster)
        # centroids = [ClusterUtils.FindClusterCentroidFromDists(x,dists) for x in splitClusts]
        for index, cluster in enumerate(splitClusts):
            children = cluster.GetPoints() 
            pts = [x.GetData() for x in children]  
            self.clustdict[index+1] = pts 
            for pt in pts:
                self.cdict[pt] = [index + 1] 
                if pt == pts[0]:
                    self.cdict[pt].append("true")
                else:
                    self.cdict[pt].append("flase")
                    




