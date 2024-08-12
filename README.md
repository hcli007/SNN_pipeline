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

Detailed help information can be viewed through `synetfind -h`.

After the step is executed, you will receive a **cleaned-network** file, this is a two column file composed of nodes hit by the model. At the same time, you will receive a **SynNet_f** file containing collinearity information for all species. Subsequent modules require the input of these two files.

<a name="SNN"></a>
### Step 5: Create a synteny neighborhood network
Before creating SNN, you need to format the syn network files.

`synetprefix -n your_path/SynNet_file -o your_output_path`

After running the previous module, you will receive a **prefix file**. 

Create an SNN using the following command:

`synetcontext -i your_path/species_list -e your_path/cleaned-network_file -n your_path/SynNet_f_file -N your_path/prefix_file -d your_pep_bed_path -o your_output_path -S 10 --block_stat --KEGG`

Or provide a custom node ID list:

`synetcontext -i your_path/species_list -I your_path/ID_list_file -n your_path/SynNet_f_file -N your_path/prefix_file -d your_pep_bed_path -o your_output_path -S 10 --block_stat --KEGG`
Enter the maximum range of flanking genes you want to search for in `- S`.

if you input `--block_stat` , module will statistically analyze the blocks that make up SNN:

```
	Malus x domestica	Prunus persica	Arabidopsis thaliana	Arabidopsis lyrata	Amborella trichopoda
Malus x domestica	26.037037037037038	29.246376811594203	15.6	15.773584905660377	10.45
Prunus persica	29.246376811594203	18.166666666666668	15.075	15.72972972972973	12.38888888888889
Arabidopsis thaliana	15.6	15.075	15.666666666666666	32.97959183673469	8.333333333333334
Arabidopsis lyrata	15.773584905660377	15.72972972972973	32.97959183673469	16.1875	8.75
Amborella trichopoda	10.45	12.38888888888889	8.333333333333334	8.75	6.0
```

if you input '--KEGG', the module will use kofamscan to annotate the proteins corresponding to nodes in SNN.
