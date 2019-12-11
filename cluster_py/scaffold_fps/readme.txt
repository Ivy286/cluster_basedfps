############################################
   ##### based on murcko scaffold #####
############################################
python murcko_scaffold
--输入支持类型：csv
--输出支持类型：csv

生成的文件示例：
Smiles	clustered_scaffold	ID	cluster_label	num_count	Energy
[3H]C(C)(C)N(c1cc(CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)c(OC)cn1)S(C)(=O)=O	C1CCC(CCC2CCCC(CC3CCCCC3)C2)CC1	r2_gen_63144	83	5863	-10.37
C=C([11C@@H](C)S(C)(=O)=O)N(C)c1ncc([SH](C)(C)=O)cc1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	C1CCC(CCC2CCCC(CC3CCCCC3)C2)CC1	r2_gen_73290	83	5863	-8.63
[C-]c1cccc(C(C)(C)C)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	C1CCC(CCC2CCCC(CC3CCCCC3)C2)CC1	r2_gen_235206	83	5863	-8.93
[C-]c1cccc(C(=C)S(C)(=O)=O)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	C1CCC(CCC2CCCC(CC3CCCCC3)C2)CC1	r2_gen_46759	83	5863	-9.91
[C]c1cc(C)cc(N(C)S(C)(=O)=O)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	C1CCC(CCC2CCCC(CC3CCCCC3)C2)CC1	r2_gen_198664	83	5863	-10.41


############################################
        ##### based on fps #####
############################################
python ClusterMurckoFps.py --in=cluster_test.csv --out=output.csv --fp=mo --cutoff=0.3
-- 输入支持三种类型：smi，csv，sdf
-- 输出支持两种类型：csv，sdf

生成文件的示例：
Smiles	ID	fp_cluster	isCentroid
[3H]C(C)(C)N(c1cc(CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)c(OC)cn1)S(C)(=O)=O	r2_gen_63144	1	flase
C=C([11C@@H](C)S(C)(=O)=O)N(C)c1ncc([SH](C)(C)=O)cc1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_73290	11	TRUE
[C-]c1cccc(C(C)(C)C)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_235206	1	flase
[C-]c1cccc(C(=C)S(C)(=O)=O)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_46759	1	flase


############################################
##### based on murcko scaffold and fps #####
############################################
Usage: python *.py [option]

Options:
  -h, --help            show this help message and exit
  --input1=INPUT1       format: csv
  --input2=INPUT2       format: csv
  --fp=FP               fingerprint type: tp,mc,mo (Topological Fingerprints,
                        MACCS Keys, Morgan Fingerprints), default is mc
  --radius=RADIUS        the radius of the Morgan fingerprint, default is 2
  --algorithm=ALGORITHM
                        cluster algorithm :b,m (Butina, Murtagh), default is b
  --cutoff=CUTOFF       distThresh(0-1),elements within this range of each
                        other are considered to be neighbors, needed for
                        Butina cluster algorithm, default is 0.5
  --nclusts=NCLUSTS     number of clusters, needed for Murtagh cluster
                        algorithm, default is 1
  --murtype=MURTYPE     Method for Murtagh:Wards, SLINK, CLINK, UPGMA, needed
                        when Murtagh is set as algorithm, default is Wards
  --out=OUTPUT          output sdf file or csv file

命令示例：python ClusterMurckoFps.py --input1=xxx.csv --input2=xxx.csv --out=output.csv --fp=mo --cutoff=0.3  (0.82min/1w 1core)

input1.csv文件示例：
drugs	target_pdb_id	targets	drug_smiles
r2_gen_63144		FAK	[3H]C(C)(C)N(c1cc(CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)c(OC)cn1)S(C)(=O)=O
r2_gen_73290	3BZ3	FAK	C=C([11C@@H](C)S(C)(=O)=O)N(C)c1ncc([SH](C)(C)=O)cc1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F
r2_gen_243538		FAK	C=C1Cc2c(ccc(C)c2CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)N(S(C)(=O)=O)C(C)C1
r2_gen_223315	3BZ3	FAK	C=C1CC(C)(c2ccncc2CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)C1=O

input2.csv文件示例：
drug_id	pdb_id	Energy	pose_file_name
r2_gen_128599	3BZ3	-9.58	r2_gen_128599_3BZ3_best_docking_pose.mol2
r2_gen_35880	3BZ3	-7.83	r2_gen_35880_3BZ3_best_docking_pose.mol2
r2_gen_33877	3BZ3	-6.69	r2_gen_33877_3BZ3_best_docking_pose.mol2
r2_gen_180761	3BZ3	-6.01	r2_gen_180761_3BZ3_best_docking_pose.mol2
r2_gen_154613	3BZ3	-8.44	r2_gen_154613_3BZ3_best_docking_pose.mol2

最后生成的文件有两个：
based_scaffold_cluster.csv	基于骨架的聚类结果
output.csv	最后输出文件，包含基于骨架和分子指纹聚类结果

output生成文件示例：
Smiles	ID	fp_cluster_label	isCentroid	clustered_scaffold	scaffold_cluster_label	num_count	Energy
CNC(=O)c1ccc(Nc2ncc(C(F)(F)F)c(NCC3=C4C(C)=C(F)C(C)=CCN4S(C)=Cc4ccccc43)n2)cc1	r2_gen_287480	1	TRUE	C1CCC(CC2CCCC(CCC3C4CCCCCC4CCC4CCCCC43)C2)CC1	0	1	-6.76
CNC(=O)c1ccc(Nc2ncc(C(F)(F)F)c(NCc3sccc3S(=O)(=O)Cc3cc(C)c(N)nc(=O)c3)n2)cc1	r2_gen_40412	2	TRUE	C1CCCC(CCC2CCCC2CCC2CCCC(CC3CCCCC3)C2)CC1	1	1	-8.69
CNC(=O)c1ccc(Nc2ncc(C(F)(F)F)c(NCc3c(-c4cccnc(=O)c4)ccc(C)c3N(C)S(C)(=O)=O)n2)cc1	r2_gen_290607	4	TRUE	C1CCC(CC2CCCC(CCC3CCCCC3C3CCCCCC3)C2)CC1	2	2	-6.8


