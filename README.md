# ALFATClust - <span style="color:yellow">AL</span>ignment-<span style="color:yellow">F</span>ree <span style="color:yellow">A</span>daptive <span style="color:yellow">T</span>hreshold <span style="color:yellow">Clust</span>ering
## Overview
Biological sequence clustering tool with dynamic threshold for individual clusters. Suitable for clustering multiple groups of homologous sequences.

## Citation
Chiu, J.K.H., Ong, R.TH. Clustering biological sequences with dynamic sequence similarity threshold. *BMC Bioinformatics* 23, 108 (2022). https://doi.org/10.1186/s12859-022-04643-9

## Release update
- 2022-10-14:<br/>
    - Fixed the false reporting of errors raised by certain E notation values in the form of xEy where x is integer.
- 2022-09-26:<br/>
    - Added error handling for sequence distance estimation process.
- 2022-06-28:<br/>
    - Replaced a Python library function that is not supported by Mac OS.
    - Updated the sequence cluster evaluation module to fix the runtime error that is raised during direct execution in host machine.
    - Added the conda environment file.

## Sequence file requirements
The input sequence file must be:
1. Consisting of either DNA or protein sequences
2. In FASTA format

## Installation
- **Method 1: Docker**<br/>
    ALFATClust is available as a Docker package, from which a Docker image can be built and then use it to create a Docker container as a virtual environment. Details of Docker and its installation can be found [here](https://www.docker.com/).<br/><br/>
    **Step 1:** The Docker image can be built via either the repository URL or the local directory:<br/>
    - **Option 1: Build with repository URL**<br/>
        The following command builds a Docker image managed by the host Docker engine:
        > docker build -t *\<image name>* github<span>.com/phglab/ALFATClust</span><br/>

        The Docker image built will be named as *\<image name>*. Users may name their own images such as "phglab/alfatclust".<br/>
    - **Option 2: Build locally**<br/>
        Image can also be built after cloning or downloading the ALFATClust repository to local directory:
        > docker build -t *\<image name>* *\<path of dockerfile>*<br/>

        *\<path of dockerfile>* locates the file path of "Dockerfile". If the current path is the root directory of ALFATClust (i.e. the original parent directory of Dockerfile), just specify "." for it.<br/>

    **Step 2:** Once the image is built, a Docker container can be created from it:<br/>
    > docker run -it --mount type=bind,src=*\<path of host data directory>*,dst=*\<path of container data directory>* --name *\<container name>* *\<image name>*

    A directory (specified in *\<path of host data directory>*) on the host machine can be mounted to the Docker container as *\<path of container data directory>*. An absolute (full) path is recommended. Users may also name their own containers such as "alfatclust". *\<image name>* is the Docker image named in step 1.

    Refer to [here](https://docs.docker.com/engine/reference/commandline/service_create/#add-bind-mounts-volumes-or-memory-filesystems) for more mounting options.

    **Step 3:** To start an existing container:
    > docker start -ai *\<container name>*

    *\<container name>* is the Docker container named in step 2.

    **Step 4:** To exit the container, run the following command in the terminal running it:
    > exit

- **Method 2: conda**<br/>
    Machines having [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed can set up a conda environment to run ALFATClust.<br/>

    **Step 1:** Update conda to the latest version, e.g. for Miniconda it can be updated via terminal:
    > conda update conda

    **Step 2:** Create a new conda environment for ALFATClust execution using the environment file "alfatclust-conda.yml":
    > conda env create -n *\<environment name>* --file *\<path of alfatclust-conda.yml>*

    The conda environment created will be named as *\<environment name>*. Users may name their own environments such as "alfatclust_env". *\<alfatclust-conda.yml>* is located under the root directory of ALFATClust.

    **Step 3:** Activate the created conda environment by:
    > conda activate *\<environment name>*

    *\<environment name>* is the environment named in step 2.

    **Step 4:** To deactivate (exit) the conda environment, run the following command:
    > conda deactivate

- **Method 3: Direct execution in host**<br/>
    The source codes of ALFATClust are under the directory "main". Simply copy the contents in the "main" folder to a local folder. Users may consider adding the path of this local folder to PATH variable. Also, make sure the following tools and libraries are properly installed and can be invoked by ALFATClust. The version tested is indicated in parentheses.

    - Python runtime:
        - Python 3 (3.7.1)

    - Python packages:
        - numpy (1.20.1)
        - scipy (1.6.0)
        - pandas (1.2.2)
        - biopython (1.78)
        - python-igraph (0.8.3)
        - leidenalg (0.8.3)

    - Third-party tool:
        - Mash (2.2.2)
        - MMseqs2 (12.113e3)

Mash [1] can be installed using apt in Ubuntu; an alternative is to download its source codes (requires compilation) or binaries from [here](https://github.com/marbl/Mash/releases). MMseqs2 [2] is used for pre-clustering only. Make sure they are included in the system path.

## Usage
***Command***
- Docker:
    > alfatclust [optional arguments] -i *\<sequence file path>* -o *\<output cluster file path>*

    Note: Both full and relative file paths are accepted.

- conda/Direct execution (assuming the current working directory is the root directory of ALFATClust):
    > ./alfatclust.py [optional arguments] -i *\<sequence file path>* -o *\<output cluster file path>*

    Note: When the current working directory is somewhere else, locate "alfatclust.py" using a (full/relative) path instead.

***Mandatory arguments***

| Argument name                             | Description                                                |
| ----------------------------------------- | ---------------------------------------------------------- |
| -i/--input *\<sequence file path>*        | (full/relative) input DNA/protein sequence FASTA file path |
| -o/--output *\<output cluster file path>* | (full/relative) output sequence cluster file path          |

***Optional arguments***

| Argument name                             | Description [default value]                                                                               |
| ----------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| -e/--evaluate *\<cluster eval file path>* | evaluate the clusters and export the evaluation results to (full/relative) *\<cluster eval file path>*    |
| -l/--lower *\<lower>*                     | set the lower bound of the sequence distance estimate (resolution parameter) to *\<lower>* [0.75]         |
| -d/--step *\<step>*                       | set the step size of the sequence distance estimate range to *\<step>* [0.025]                            |
| -p/--precluster                           | always run pre-clustering                                                                                 |
| -k/--kmer *\<kmer>*                       | set the Mash kmer size parameter to *\<kmer>* [DNA: 17; protein: 9]                                       |
| -s/--sketch *\<sketch>*                   | set the Mash sketch size parameter to *\<sketch>* [2000]                                                  |
| -m/--margin *\<margin>*                   | ignore any Mash distance above 1 - max(*\<lower>* - *\<margin>*, 0) [0.2]                                 |
| -f/--filter *\<filter>*                   | discard a Mash distance when its shared hash ratio is below *\<filter>*, NOT recommended                  |
| -t/--thread *\<thread>*                   | set the number of threads to *\<thread>* (for Mash and cluster evaluation only) [all available CPU cores] |
| -S/--seed *\<seed>*                       | set the seed value to *\<seed>*                                                                           |
| -h/--help                                 | show help message and exit                                                                                |

***Evaluation report***

The evaluation report consists of the following columns:

| Column name                         | Description                                                                     |
| ----------------------------------- | ------------------------------------------------------------------------------- |
| Cluster Id                          | Cluster Id for the non-singleton cluster                                        |
| No. of sequences                    | Number of sequences in the cluster                                              |
| Average sequence identity           | Cluster average pairwise sequence identity with respect to the center sequence* |
| Min. sequence identity              | Cluster minimum pairwise sequence identity with respect to the center sequence* |
| Center sequence                     | Representative center sequence selected for cluster                             |
| Sequence for min. sequence identity | Sequence showing the lowest pairwise sequence identity with the center sequence |

*Sequence identity = number of matched nucleotides or amino acids / (alignment length - terminal gaps)

***Configuration file***

The configuration file "settings.cfg" is located under directory "/usr/local/bin/phglab/alfatclust" inside the Docker container, or under the same host directory as the main Python script "<span>alfatclust.py</span>". It consists of the default values for the following parameters organized into various categories:

- *EstimatedSimilarity*
    - *High*: Upper bound of the sequence distance estimate (resolution parameter) range
    - *Low*: Lower bound of the sequence distance estimate range
    - *StepSize*: Step size of the sequence distance estimate range during clustering

- *Threshold*
    - *Precluster*: No. of sequences above which pre-clustering is performed to partition sequences

- *DNAMash*
    - *Kmer*: DNA k-mer size for Mash
    - *Sketch*: DNA sketch size for Mash

- *ProteinMash*
    - *Kmer*: Protein k-mer size for Mash
    - *Sketch*: Protein k-mer size for Mash

- *NoiseFilter*
    - *Margin*: Any pairwise Mash distance *d* > 1 - max(*\<lower bound>* - *\<margin>*, 0) is regarded as noise and discarded

- *DNAEvaluation* (for sequence alignment during DNA cluster evaluation)
    - *MatchScore*: Score for a nucleotide match
    - *MismatchPenalty*: Penalty for a nucleotide mismatch
    - *GapOpeningPenalty*: Penalty for opening a gap
    - *GapExtensionPenalty*: Penalty for extending a gap

- *ProteinEvaluation* (for sequence alignment during protein cluster evaluation)
    - *ScoreMatrix*: Score matrix used for amino acid matching (refer to Biopython documentation for the built-in score matrices available)
    - *GapOpeningPenalty*: Penalty for opening a gap
    - *GapExtensionPenalty*: Penalty for extending a gap

## Sample datasets
The sample datasets are available in folder *sample_datasets*, which includes:

1. Antimicrobial resistance (AMR) gene datasets (data sources: ARG-ANNOT [3], CARD [4-6], and ResFinder [7])
    **argdit_nt_06feb2020_full.fa** (DNA)
    **argdit_aa_06feb2020_full.fa** (protein)

2. Non-AMR plasmid gene dataset (data source: PLSDB [8])
    **plasmid_genes_20191017.fa** (DNA)

## References
[1] Ondov, B. D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132.<br/>
[2] Steinegger, M. and J. Söding. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026.<br/>
[3] Gupta, S. K., et al. (2014). ARG-ANNOT, a New Bioinformatic Tool To Discover Antibiotic Resistance Genes in Bacterial Genomes. Antimicrobial Agents and Chemotherapy, 58(1), 212-220.<br/>
[4] Alcock, B. P., et al. (2019). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic Acids Research, 48(D1), D517-D525.<br/>
[5] Jia, B., et al. (2017). CARD 2017: expansion and model-centric curation of the comprehensive antibiotic resistance database. Nucleic Acids Research, 45(D1), D566-D573.<br/>
[6] McArthur, A. G., et al. (2013). The Comprehensive Antibiotic Resistance Database. Antimicrobial Agents and Chemotherapy, 57(7), 3348-3357.<br/>
[7] Zankari, E., et al. (2012). Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy, 67(11), 2640-2644.<br/>
[8] Galata, V., et al. (2018). PLSDB: a resource of complete bacterial plasmids. Nucleic Acids Research, 47(D1), D195-D202.
