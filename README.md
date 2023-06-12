# PGCM - a comparative study of Pan-Genome Construction Methods
This repository contains code and notebooks related to the manuscript titled "The effect of methodological considerations on the construction of gene-based plant pan-genomes" (currently under revision for GBE). Please note that the code is intended for research reproducibility, not as standalone scripts or pipelines.  

## Instructions for reproducing the analyses described in the manuscript

### Step 0: Clone the repository
```
git clone https://github.com/MayroseLab/PGCM.git
```

### Step 1: Download and prepare data
#### 1.1 - reference, assembly and annotation data
Start by downloading the input data required for pan-genome construction.  
For _A. thaliana_:
```
wget https://datadryad.org/stash/downloads/file_stream/2337370 -O A_thaliana_PG_construct_data.tar.gz
tar -zxvf A_thaliana_PG_construct_data.tar.gz
```
For soybean:
```
wget https://datadryad.org/stash/downloads/file_stream/2337371 -O soybean_PG_construct_data.tar.gz
tar -zxvf soybean_PG_construct_data.tar.gz
```
For data related to the meta-analysis of pan-genomes:
```
wget https://datadryad.org/stash/downloads/file_stream/2347427 -O PG_meta.tar.gz
tar -zxvf PG_meta.tar.gz
```
#### 1.2 - sequence data
Download the FastQ files according to the identifiers in Supplementary Table S2 of the manuscript. For _A. thaliana_ you can use FTP or [Kingfisher](https://github.com/wwood/kingfisher-download) to download from ENA. For soybean, you need to download from GSA.  
Once you obtain the FastQ files, subsample them to the required sequencing depth. E.g., to produce 50x data for the _A. thaliana_ accession An-1:
```
zcat ERR3624579_1.fastq.gz | head -118811880 | gzip > RESULT/per_sample/An-1/data/ERR3624579_1.fastq.gz
zcat ERR3624579_2.fastq.gz | head -118811880 | gzip > RESULT/per_sample/An-1/data/ERR3624579_2.fastq.gz
```

### Step 2: Construct pan-genomes
Pan-genomes need to be constructed using the software [Panoramic](https://github.com/MayroseLab/Panoramic).
#### 2.1 - install and configure Panoramic
Follow the instructions [here](https://github.com/MayroseLab/Panoramic/wiki/Panoramic-setup) to setup Panoramic. Make sure it runs (use the test runs as explained) before proceeding.
#### 2.2 - configure the Panoramic runs
To construct each of the pan-genomes detailed in Supplementary Table S1, you will need to run Panoramic with the right configuration. To do that, copy the relevant yml and TSV file from the `panoramic_conf` directory in this repository, and modify it to match the paths on your filesystem. In addition, you will need the files `LQ_samples_info.tsv` and `HQ_samples_info.tsv`, which you do not need to modify. Further information regarding Panoramic configurations can be found [here](https://github.com/MayroseLab/Panoramic/wiki/Running-Panoramic).
#### 2.3 - prepare sequencing data
You will need to copy the subsampled FastQ files to specific locations to have Panoramic use them. Start by preparing the directory structure (e.g. for _A. thaliana_):
```
mkdir -p RESULT/per_sample/{An-1,C24,Cvi-0,Eri,Kyo,Ler,Sha}/data
```
Next, copy/link the FastQ files for each sample into the corresponding directories. Make sure you keep the naming convention accordint to the identifiers in `LQ_samples_info.tsv`. For example, for An-1:
```
ln ERR3624579_1.fastq.gz RESULT/per_sample/An-1/data
ln ERR3624579_2.fastq.gz RESULT/per_sample/An-1/data
```
  
If you are constructing one of the "HQ\_asm" (high-quality assemblies) or "hybrid" pan-genomes, you will need to provide Panoramic with assemblies rather than reas. To do that:
```
mkdir -p RESULT/per_sample/
```
Then, for each sample create an assembly directory and place the assembly fasta there. For example, for An-1 in the HQ_asm pan-genome:
```
mkdir -p RESULT/per_sample/An-1/RG_assembly_ERR3624579/ragtag_output/
ln An-1.chr.all.v2.0.chr0.fasta RESULT/per_sample/An-1/RG_assembly_ERR3624579/ragtag_output/ragtag.scaffolds.fasta
ln An-1.chr.all.v2.0.chr0.fasta RESULT/per_sample/An-1/RG_assembly_ERR3624579/ragtag_output/ragtag.correct.fasta
touch RESULT/per_sample/An-1/RG_assembly_ERR3624579/ragtag_output/ragtag.scaffolds.agp
```
#### 2.4 - run Panoramic
Once your configurations and data are ready, you can run Panoramic to construct the pan-genome. [This page](https://github.com/MayroseLab/Panoramic/wiki/Running-Panoramic) explains how to do that.

### Step 3: Compare pan-genomes
Pairwise comparisons between pan-genomes are performed using a dedicated Snakemake pipeline.  
All required code and configurations can be found under the `compare_pan_genomes` directory in this repository.  
To run a specific comparison:
Activate the conda environment you used for running panoramic:
```
conda activate panoramic
```
Copy the relevant configurations from the repository, e.g.:
```
cp PGCM/compare_pan_genomes/conf/At_DN_x50_vs_MTP_x50/{conf.yml,pan_genomes.tsv} ./
```
Modify the paths to match your filesystem, then use snakemake to run the comparison:
```
snakefile="PGCM/compare_pan_genomes/compare_pan_genomes.snakefile"
snakemake -s $snakefile --configfile conf.yml --latency-wait 60 --use-conda -p -j 10
```
Repeat this for all pairwise comparisons.

### Step 4: Analyze nonreference pan-genes
Nonreference genes detected in some of the pan-genomes were analyzed to determine the level of reliability.  
To reproduce this analysis, use the code and configurations under `analyze_nonref`. Usage instructions are the same aas in Step 3, except the configurations and snakefile are different.

### Step 5: Summarize the analysis
Once you have completed all the analyses described above, you can summarize the results and create the figures by running the Jupyter notebooks found under `notebooks`.  
Start by creating the required conda environment using conda or mamba:
```
mamba env create -f PGCM/notebooks/jupyterlab.yml
```
Activate the environment:
```
conda activate PGCM_notebooks
```
Start the Jupyter server:
```
jupyter lab
```
Open the notebooks through the Jupyter GUI, modify the paths to match your filesystem and run.
