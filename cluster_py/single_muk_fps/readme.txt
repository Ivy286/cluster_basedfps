
##### based on murcko scaffold #####
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



##### based on fps #####
python ClusterMurckoFps.py --in=cluster_test.csv --out=tmptmp.csv --fp=mo --cutoff=0.3
-- 输入支持三种类型：smi，csv，sdf
-- 输出支持两种类型：csv，sdf

生成文件的示例：
Smiles	ID	fp_cluster	isCentroid
[3H]C(C)(C)N(c1cc(CNc2nc(Nc3ccc(C(=O)NC)cc3)ncc2C(F)(F)F)c(OC)cn1)S(C)(=O)=O	r2_gen_63144	1	flase
C=C([11C@@H](C)S(C)(=O)=O)N(C)c1ncc([SH](C)(C)=O)cc1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_73290	11	TRUE
[C-]c1cccc(C(C)(C)C)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_235206	1	flase
[C-]c1cccc(C(=C)S(C)(=O)=O)c1CNc1nc(Nc2ccc(C(=O)NC)cc2)ncc1C(F)(F)F	r2_gen_46759	1	flase


