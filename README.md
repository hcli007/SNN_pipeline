[![Latest PyPI Version](https://img.shields.io/pypi/v/snn-pipeline.svg)](https://pypi.org/project/snn-pipeline/)
[![Python package](https://github.com/hcli007/SNN_pipeline/actions/workflows/python-package.yml/badge.svg)](https://github.com/hcli007/SNN_pipeline/actions/workflows/python-package.yml)
# A pipeline for constructing and analyzing synteny network
A pipeline for constructing and analyzing synteny network.
## Requirements
python>=3.8

MCScanX

diamond>=0.9.14

HMMER>=3.1

python library:

igraph>=0.1.14

networkx>=2.6.3

numpy>=1.21.2

pandas>=1.4.1

- [install](#install)
- [prepare data](#preparation)
- [build species_list](#buildlist)
- [build network_database](#builddatabase)
- [extract target network](#extract)
- [build SNN](#SNN)
- [example](#example)
- [file explanation](#explanation)

## Usage
<a name="install"></a>
### install snn_pipeline
Create a dedicated environment for SNN (Recommended).

`conda create -n snn_pipline`

`conda activate snn_pipline`

Install python (Recommend installing versions after Python 3.8).

`conda install python`

Install snn_pipline.

`pip install snn_pipeline`

You can also manually install this process.

`tar -zxvf snn_pipeline-0.1.0.tar.gz`

`cd snn_pipeline-0.1.0`

`pip install .`

View the installation location of the package.

`pip show snn_pipeline`

Install [MCScanX](https://github.com/wyp1125/MCScanX), [diamond](https://github.com/bbuchfink/diamond), [hmmer](https://github.com/EddyRivasLab/hmmer) and [kofam_scan](https://github.com/takaram/kofam_scan) and **ensure that the above software has been added to the environment variables**.

Try:
`conda install bioconda::diamond`

`conda install bioconda::hmmer`

Test files in 
`/test-data`.


<a name="preparation"></a>
### Step 1: Prepare the bed file and pep file of the species 
Firstly, please prepare bed files for all species and fasta format files for proteins. Determine a **four character abbreviation** for each species, 
for example, _Arabidopsis thaliana_ can be set as `ath_`, and _Amborella trichopoda_ can be set as `atr_`.
The bed file contains the location information of cDNA, and the format requirements are as follows: the first column is the chromosome number, which needs to be preceded by the species abbreviation and the first two characters must be unique to the species (**distinguished by uppercase and lowercase letters**)

bed file of _Arabidopsis thaliana_:

`aThChr1	ath_AT1G01010	3631	5899`

bed file of _Amborella trichopoda_:

`Atrscaffold00001	atr_scaffold00001.1atr	3379	6049`

The second column is the ID number of cDNA, which needs to include the species abbreviation, such as ath_T1G01010. The third column represents the starting position of cDNA, and the fourth column represents the ending position of cDNA. For the fasta file of proteins, please ensure that the protein ID in the file corresponds to the ID in the bed file, and change the file suffix to. pep. Ensure that the bed file and pep file with the same abbreviation, for example, the bed and pep files corresponding to Arabidopsis thaliana are **ath.bed** and **ath.pep**, respectively.

<a name="buildlist"></a>
### Step 2: Create a name list file that includes species abbreviations and taxonomic information
Fill in the names of the pep file and bed file corresponding to the species in the first column of the name list, the abbreviation of the ID in the second column of the name list, and the classification information of the species in the third, fourth, fifth, and sixth columns.
for example:
```
#file_name	ID_abbreviation	Clade	Order	Family	Species_Name
Mdgd	Mdgd	Super-Rosides	Rosales	Rosaceae	Malus x domestica
ppe	ppe_	Super-Rosides	Rosales	Rosaceae	Prunus persica
ath	ath_	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana
Alyr	Alyr	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis lyrata
atr	atr_	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda
```

<a name="builddatabase"></a>
### Step 3: Constructing a species overall collinear network
Run the **synetbuild** module.

`synetbuild -i your_path/species_list -d  your_pep_bed_path -o your_output_path`

You can input some additional options to adjust the parameters for building the network.

`synetbuild -i your_path/species_list -d  your_pep_bed_path -o your_output_path -k 5 -s 5 -m 10 -p 4 -D -T`

If you enter the above parameters, a collinear network will be constructed with the top **5** results **hit** by blastp, Minimum **5** of **Anchors** for a synteny block, and Maximum **10** of Genes allowed as the **GAP** between Anchors as parameters, and **4** threads will be used for blastp. Inputting` -D` will run MCScanX's duplicate_gene_classifier module, and inputting` -T `will run the detectability collineear_tandem_arrays module, which can help you better search for tandem genes. Detailed help information can be viewed through `synetbuild -h`

After the step is executed, you will receive a **SynNet** file, which is the total network file and will be used for subsequent analysis.

<a name="extract"></a>
### Step 4: Extract the networks that interest you
Extracting subnetworks using the **synetfind** module. Take the **SynNet** file obtained in the previous step as input for `-n`, use `-m` input your protein model file , and you will get a network composed of nodes hit by the model.

`synetfind -i your_path/species_list -m your_path/hmm_file -d your_pep_bed_path -n your_path/SynNet_file -o your_output_path`

You can input the threshold of E-value through `-E` to make the search conditions more stringent or relaxed (default is 0.001). 

`synetfind -i your_path/species_list -m your_path/hmm_file -d your_pep_bed_path -n your_path/SynNet_file -o your_output_path -E 0.01`

If you have a specific list of IDs：
`synetfind -l your_path/your_id_list_file -n your_path/SynNet_file -o your_output_path`

Detailed help information can be viewed through `synetfind -h`.

After the step is executed, you will receive a **cleaned-network** file, this is a two column file composed of nodes hit by the model. At the same time, you will receive a **SynNet_f** file containing collinearity information for all species. Subsequent modules require the input of these two files.

<a name="SNN"></a>
### Step 5: Create a synteny neighborhood network
Before creating SNN, you need to format the syn network files.

`synetprefix -n your_path/SynNet_file -o your_output_path`

After running the previous module, you will receive a **prefix file**. 

Create an SNN using the following command:

`synetcontext -i your_path/species_list -e your_path/cleaned-network_file -b your_path/SynNet_f_file -N your_path/prefix_file -d your_pep_bed_path -o your_output_path -S 10 --block_stat --KEGG`

In the example, the `-e` parameter is used, which will extract the node information from the edge file and search for the syntenic genes upstream and downstream of these nodes. Please input the synteny information file for all species with the `-n` option and enter the prefix file with the `-N` option. and enter the maximum range of flanking genes you want to search for in `-S`.

If you have classified genes based on phylogenetic relationships or community structures in a network, or if you have specific nodes of interest, you can use the `-I` instead of `-e` as input:

`synetcontext -i your_path/species_list -I your_path/ID_list_file -n your_path/SynNet_f_file -N your_path/prefix_file -d your_pep_bed_path -o your_output_path -S 10 --block_stat --KEGG`

The `-I` parameter requires you to prepare a customized list of nodes:
```
AlyrAL1G11030
AlyrAL1G11030
AlyrAL1G51660
AlyrAL1G11030
AlyrAL1G35790
```

if you input `--block_stat` , module will statistically analyze the blocks that make up SNN, the following content is from the file **block_len_avg_stat.tsv**. The results represent the block composition of the example SNN, expressed as a matrix of average block lengths:

```
	Malus x domestica	Prunus persica	Arabidopsis thaliana	Arabidopsis lyrata	Amborella trichopoda
Malus x domestica	26.037037037037038	29.246376811594203	15.6	15.773584905660377	10.45
Prunus persica	29.246376811594203	18.166666666666668	15.075	15.72972972972973	12.38888888888889
Arabidopsis thaliana	15.6	15.075	15.666666666666666	32.97959183673469	8.333333333333334
Arabidopsis lyrata	15.773584905660377	15.72972972972973	32.97959183673469	16.1875	8.75
Amborella trichopoda	10.45	12.38888888888889	8.333333333333334	8.75	6.0
```

if you input '--KEGG', the module will use kofamscan to annotate the proteins corresponding to nodes in SNN.
```
#ID	Clade	Order	Family	Species_Name	KEGG_definition	Greedy_cluster	Infomap_cluster
atr_scaffold00040.1atr	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda	None	119	260
ath_AT1G43710	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana	serine decarboxylase [EC:4.1.1.-]	119	260
atr_scaffold00040.6atr	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda	None	120	256
ath_AT1G43650	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana	None	120	256
atr_scaffold00040.18atr	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda	tubby and related proteins	121	255
ath_AT1G43640	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana	tubby and related proteins	121	255
atr_scaffold00040.25atr	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda	None	122	254
ath_AT1G43630	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana	None	122	254
atr_scaffold00040.35atr	Basal-Angiosperm	Amborellales	Amborellaceae	Amborella trichopoda	cinnamyl-alcohol dehydrogenase [EC:1.1.1.195]	123	253
ath_AT1G43620	Super-Rosides	Brassicales	Brassicaceae	Arabidopsis thaliana	sterol 3beta-glucosyltransferase [EC:2.4.1.173]	123	253
```

<a name="example"></a>
Next, we will demonstrate how to construct the SNN for OSC (Oxidosqualene cyclase) in Brassicaceae species. Example datasets and results can be obtained from the site`test_data/Brassicaceae`.
Due to upload limitations, the PEP data in the example has been compressed using zip. Please decompress it before use by running`unzip your_path/Brassicaceae/*.zip`.

First, construct the overall synteny network for the Brassicaceae species.

`synbuild -i your_path/test_data/Brassicaceae/species_list -d your_path/test_data/Brassicaceae/ -o your_path`

The OSC generally contains two domains, SQHop_cyclase_C and SQHop_cyclase_N, so we need to determine which sequences in the dataset contain both of the above domains. 

We can obtain a list of sequence IDs in the dataset that contain the aforementioned domains by running `synetfind` twice.

`synetfind -i your_path/test_data/Brassicaceae/species_list -m your_path/test_data/Brassicaceae/SQHop_cyclase_C.hmm -d your_path/test_data/Brassicaceae/ -n your_path/SynNetBuild/SynNet-k6s5m25 -o your_output_path`
`synetfind -i your_path/test_data/Brassicaceae/species_list -m your_path/test_data/Brassicaceae/SQHop_cyclase_N.hmm -d your_path/test_data/Brassicaceae/ -n your_path/SynNetBuild/SynNet-k6s5m25 -o your_output_path`

The next step is to obtain the intersection of the two lists.

`cat your_path/SynNet_SQHop_cyclase_C.hmm/all.hmm_genelist your_path/SynNet_SQHop_cyclase_N.hmm/all.hmm_genelist > SQHop_cyclase_C_N.merge`

`sort SQHop_cyclase_C_N.merge| uniq -d > SQHop_cyclase_C_N.namelist`

Similarly, obtain the intersection of the hit regions.

`cat your_path/SynNet_SQHop_cyclase_C.hmm/SQHop_cyclase_C.hmm.genelist_SynNet_f your_path/SynNet_SQHop_cyclase_N.hmm/SQHop_cyclase_N.hmm.genelist_SynNet_f > SQHop_cyclase_C_N.merge_SynNet_f`

`SQHop_cyclase_C_N.merge_SynNet_f| uniq -d > SQHop_cyclase_C_N.hmm.genelist_SynNet_f`

Preprocess the SynNet file.

`synetprefix -n your_path/SynNetBuild/SynNet-k6s5m25 -o your_output_path`

Use `SQHop_cyclase_C_N.namelist` as the input list to build OSC-SNN.

`synetcontext -i your_path/test_data/Brassicaceae/species_list -l your_path/SQHop_cyclase_C_N.namelist -b your_path/SynNetBuild/SQHop_cyclase_C_N.hmm.genelist_SynNet_f -N your_path/SynNetTrans/SynNet-k6s5m25.prefix -d your_path/test_data/Brassicaceae/ -o your_output_path -S 10 --block_stat --block_mining --KEGG`

The command will generate multiple folders and files.
<a name="explanation"></a>
**Here is the explanation of the files generated in the process:**

`block_info` file stores the block and collinear gene information of SNN，

`result`file contains information about collinear genes，

`infomap.2col.tsv` file saves the infomap algorithm classification information for each node，

`infomap.size_fq.tsv` file saves the frequency information of each category under the infomap algorithm，

`sp_cluster.tsv` file compiles the species classification information of the nodes, as well as the information from the infomap algorithm and the greedy algorithm，

`sp_KO_cluster.tsv` file adds the results of KEGG annotation based on the classification information，

`SQHop_cyclase_C_N.namelist_k1.nodes.pep` file contains the protein sequences of all nodes in the SNN，

`K_core_result` folder contains points, edges, degree, and greedy classification information generated by different k-values in the K-core algorithm for the network.

The results of the `--block_stat` are saved in the `block_stat` folder,

`block_num_stat.tsv` file statistics the number of collinearity blocks generated between species in the dataset,

`block_len_sum_stat.tsv` file statistics the total length of collinearity blocks generated between species, measured by the number of genes,

`block_len_avg_stat.tsv` file statistics the average length of the blocks generated between species.

The results of `--block_mining` command are saved in the `specific_block` folder,

`block_info.tsv` file contains information about the blocks and the classification information of the genes within synteny blocks,

`block_info_matrix.tsv` file contains the feature matrix derived from the classification information,

`block_info_score.tsv` file contains the block specificity scores calculated from the feature matrix.

The results of `--KEGG` command are saved in `KEGG_annotation folder`, which contains the annotation information for the nodes,

The above result files with `.tsv` suffix are saved in `/test_data/Brassicaceae/result` directory of this site.

In version 0.2.1 of snn_pipeline, a new collinearity mapping feature called synetmapping has been added, which can match synteny information to corresponding chromosomal blocks

For example, we have obtained the specific block score file `block_info_score.tsv` through the `--block_mining` function of `synetcontext`, and after sorting the scores:

```
target_gene	block_id	block_specificity_score
AlyrAL5G19180&ath_AT3G29255	ath_Alyr287	542.2441494
AlyrAL2G25240&ath_AT1G66960	ath_Alyr85	542.2441494
ath_AT3G12070&thh_10021170m	thh_ath244	542.2441494
AlyrAL5G23630&ath_AT3G45130	ath_Alyr286	523.2628036
AlyrAL8G20190&ath_AT5G48010	ath_Alyr492	505.0476368
gra_005G208400&vvi_01000415001	gra_vvi398	475.4308418
AALP_AA3G344300&thh_10022563m	thh_AALP119	430.0449313
AALP_AA3G344300&bol_012476	bol_AALP195	411.0423863
AALP_AA5G078800&AlyrAL5G23630	Alyr_AALP234	404.9909972
AALP_AA5G078800&ath_AT3G45130	ath_AALP207	404.9909972
tca_029480&vvi_01021474001	tca_vvi478	399.884056
Bostr.0556s0559&cru_10016727	cru_Bost513	395.3533108
bnp_BnaA04g13460D&bnp_BnaC04g35590D	bnp1071	392.320109
bol_012476&thh_10022563m	bol_thh187	382.5703676
gra_011G044500&tca_029480	gra_tca938	380.1720883
gra_011G044500&vvi_01021474001	gra_vvi1058	380.1720883
Csat104780146&cru_10016727	cru_Csat228	379.291982
Csat104713532&Csat104751748	Csat883	367.2519706
Bostr.0556s0559&Csat104780146	Csat_Bost1483	359.5800143
bnp_BnaA04g13460D&bra_D01377	bnp_bra673	358.4514262
```

We select the genes ath_AT3G29255, ath_AT1G66960, ath_AT3G45130, and ath_AT5G48010, and examine the synteny hits in the regions where these genes are located:

`synetmapping -i test_data/Brassicaceae/species_list -l test_data/Brassicaceae/high_score_ath.namelist --bed test_data/Brassicaceae/ath.bed -n your_path/SynNetBuild/SynNet-k6s5m25  -o your_output_path -S 50 -p 4`

We will obtain the `stat.matrix.tsv` file for these genes, for example, `ath_AT5G48010_50.stat.matrix.tsv`：

```
ID	Super_Rosids	Super_Asterids	Basal_Eudicots	Monocots	Magnoliids	Basal_Angiosperm
ath_AT5G47530	41	0	0	0	0	0
ath_AT5G47540	54	0	0	0	0	0
ath_AT5G47550	24	0	0	0	0	0
ath_AT5G47560	29	0	0	0	0	0
ath_AT5G47570	30	0	0	0	0	0
ath_AT5G47580	27	0	1	0	0	0
ath_AT5G47590	11	0	0	0	0	0
ath_AT5G47600	15	0	0	0	0	0
ath_AT5G47610	44	0	0	0	0	0
ath_AT5G47620	31	0	0	0	0	0
ath_AT5G47630	19	0	0	0	0	0
ath_AT5G47635	47	0	0	0	0	0
ath_AT5G47640	17	0	0	0	0	0
ath_AT5G47650	33	0	0	0	0	0
ath_AT5G47660	30	0	0	0	0	0
ath_AT5G47670	4	0	0	0	0	0
ath_AT5G47680	28	0	0	0	0	0
ath_AT5G47690	31	0	0	0	0	0
ath_AT5G47700	12	0	0	0	0	0
ath_AT5G47710	25	0	0	0	0	0
ath_AT5G47720	23	0	0	0	0	0
ath_AT5G47730	32	0	0	0	0	0
ath_AT5G47740	30	0	0	0	0	0
ath_AT5G47750	26	0	0	0	0	0
ath_AT5G47760	25	0	0	0	0	0
ath_AT5G47770	32	0	0	0	0	0
ath_AT5G47780	15	0	0	0	0	0
ath_AT5G47790	12	0	0	0	0	0
ath_AT5G47800	10	0	0	0	0	0
ath_AT5G47810	11	0	0	0	0	0
ath_AT5G47820	23	0	0	0	0	0
ath_AT5G47830	19	0	0	0	0	0
ath_AT5G47840	18	0	0	0	0	0
ath_AT5G47850	9	0	0	0	0	0
ath_AT5G47860	16	0	1	0	0	0
ath_AT5G47870	19	0	1	0	0	0
ath_AT5G47880	21	0	0	0	0	0
ath_AT5G47890	21	0	1	0	0	0
ath_AT5G47900	16	0	1	0	0	0
ath_AT5G47910	20	0	1	0	0	0
ath_AT5G47920	23	0	1	0	0	0
ath_AT5G47928	0	0	0	0	0	0
ath_AT5G47930	23	0	1	0	0	0
ath_AT5G47940	15	0	0	0	0	0
ath_AT5G47950	3	0	0	0	0	0
ath_AT5G47960	7	0	0	0	0	0
ath_AT5G47970	7	0	0	0	0	0
ath_AT5G47980	3	0	0	0	0	0
ath_AT5G47990	2	0	0	0	0	0
ath_AT5G48000	3	0	0	0	0	0
ath_AT5G48010	3	0	0	0	0	0
ath_AT5G48020	17	0	0	0	0	0
ath_AT5G48030	25	0	0	0	0	0
ath_AT5G48040	18	0	0	0	0	0
ath_AT5G48050	3	0	0	0	0	0
ath_AT5G48060	16	0	0	0	0	0
ath_AT5G48070	11	0	0	0	0	0
ath_AT5G48090	15	0	0	0	0	0
ath_AT5G48100	18	0	0	0	0	0
ath_AT5G48110	0	0	0	0	0	0
ath_AT5G48120	14	0	0	0	0	0
ath_AT5G48130	13	0	1	0	0	0
ath_AT5G48140	45	0	0	0	0	0
ath_AT5G48150	28	0	1	0	0	0
ath_AT5G48160	38	0	1	0	0	0
ath_AT5G48170	30	0	1	0	0	0
ath_AT5G48175	5	0	0	0	0	0
ath_AT5G48180	42	0	1	0	0	0
ath_AT5G48190	0	0	0	0	0	0
ath_AT5G48200	0	0	0	0	0	0
ath_AT5G48205	0	0	0	0	0	0
ath_AT5G48210	12	0	0	0	0	0
ath_AT5G48220	23	0	0	0	0	0
ath_AT5G48230	30	0	1	0	0	0
ath_AT5G48240	30	0	1	0	0	0
ath_AT5G48250	52	0	1	0	0	0
ath_AT5G48270	35	0	1	0	0	0
ath_AT5G48280	0	0	0	0	0	0
ath_AT5G48290	36	0	1	0	0	0
ath_AT5G48300	19	0	1	0	0	0
ath_AT5G48310	23	0	1	0	0	0
ath_AT5G48320	1	0	0	0	0	0
ath_AT5G48330	14	0	1	0	0	0
ath_AT5G48335	40	0	0	0	0	0
ath_AT5G48340	21	0	0	0	0	0
ath_AT5G48350	1	0	0	0	0	0
ath_AT5G48360	50	0	0	0	0	0
ath_AT5G48370	22	0	0	0	0	0
ath_AT5G48375	2	0	0	0	0	0
ath_AT5G48380	23	0	0	0	0	0
ath_AT5G48385	26	0	0	0	0	0
ath_AT5G48390	16	0	0	0	0	0
ath_AT5G48400	38	0	0	0	0	0
ath_AT5G48410	0	0	0	0	0	0
ath_AT5G48420	2	0	0	0	0	0
ath_AT5G48430	10	0	0	0	0	0
ath_AT5G48440	22	0	0	0	0	0
ath_AT5G48450	24	0	0	0	0	0
ath_AT5G48460	24	0	0	0	0	0
ath_AT5G48470	25	0	0	0	0	0
ath_AT5G48480	32	0	0	0	0	0
```
