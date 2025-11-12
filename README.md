# CYP-Annotator

## Overview

CYP-Annotator is a Python pipeline designed for the automatic, functional annotation of Cytochrome P450 (CYP) sequences. It processes "Subject" datasets (genomes or transcriptomes) to identify CYP candidates by comparing them against a curated set of "Bait" (reference) sequences (https://github.com/k-georgi/CYP_Annotator/tree/main/Data).

The workflow performs the following steps:

1.  **Candidate Search**: Identifies initial candidates using BLASTp (default) or HMMER (`--use_hmmer y`).
2.  **Phylogenetic Classification**: Classifies candidates into "Ingroup" (putative CYPs) and "Outgroup" based on their phylogenetic proximity to the bait and outgroup references.
3.  **Ortholog Assignment**: Assigns ingroup candidates to known CYP families and subfamilies by calculating patristic distances (phylogenetic tree distance) to the reference baits.
4.  **Domain Analysis**: Screens candidates for the presence of conserved CYP domains using an HMM profile (`--hmm_domains`).
5.  **Motif Analysis**: Checks for conserved protein motifs via static definitions (`--protein_motifs`) and HMM profiles (`--hmm_motifs`).
6.  **Expression Analysis**: If expression data is provided (`--expression`), it is used to filter out low-expressed genes or potential pseudogenes (`--min_avg_tpm`).
7.  **Paralog Collapsing**: Optionally identifies and collapses paralog groups (`--collapse y`), retaining a single representative sequence.
8.  **Summary Generation**: Compiles all findings into a comprehensive summary table (`08_summary.txt`).
9.  **Export**: Generates tree visualizations (if `ete3` is installed) and an interactive HTML overview of all results.

This tool adapts functions and logic from the **MYB\_annotator** ([doi: 10.1186/s12864-022-08452-5](http://dx.doi.org/10.1186/s12864-022-08452-5)) and **bHLH\_annotator** ([doi: 10.1186/s12864-023-09877-2](https://doi.org/10.1186/s12864-023-09877-2)).

## Setup

### Installation of the dependencies

The following dependencies are necessary for the execution of the pipeline:

  * [Python 3](https://www.python.org/):
      * [dendropy](https://dendropy.readthedocs.io/en/main/): `pip install dendropy`
      * [pandas](https://pandas.pydata.org/docs/index.html): `pip install pandas`
      * [matplotlib](https://matplotlib.org/stable/index.html): `pip install matplotlib`
      * [seaborn](https://seaborn.pydata.org/): `pip install seaborn`
      * [ete3](http://etetoolkit.org/): `pip install ete3` (Required for tree visualization)
  * [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/): `sudo apt install ncbi-blast+`
  * [HMMER](http://hmmer.org/documentation.html): `conda install -c bioconda hmmer`
  * [MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html): `sudo apt install mafft`
  * [MUSCLE](https://www.drive5.com/muscle/): (Installation recommended)
  * [FastTree](http://www.microbesonline.org/fasttree/#Install): `sudo apt-get install -y fasttree`
  * [RAxML-NG](https://github.com/amkozlov/raxml-ng): (Precompiled binaries recommended)
  * [IQ-TREE](http://www.iqtree.org/): (Precompiled binaries recommended)

The CYP\_Annotator can be cloned from GitHub:

```batch
git clone https://github.com/k-georgi/CYP_Annotator
cd CYP_Annotator
```

## Usage

The pipeline is executed via the command line, using **exactly one** of the following input methods:

```bash
# Method 1: Using a data folder
python3 CYP_Annotator.py --data <PATH_TO_DATA_FOLDER>

# Method 2: Using explicit paths
python3 CYP_Annotator.py --baits <PATH> --subject <PATH> --baits_info <PATH>

# Method 3: Using a file collection CSV
python3 CYP_Annotator.py --file_collection <PATH_TO_FILE_COLLECTION_CSV>
```

-----

### Optional arguments

#### Optional arguments regarding input/output

|Command|Description|Default
|--|--|--
|`--outgroup <PATH>` |Path to outgroup FASTA file (optional)|`None`
|`--hmm <PATH>` |Path to bait HMM file|`None`
|`--hmm_domains <PATH>` |Path to hmm domains file (optional)|`None`
|`--hmm_motifs <PATH>` |Path to hmm motifs file (optional)|`None`
|`--protein_motifs <PATH>` |Path to protein motifs file (optional)|`None`
|`--expression <PATH>` |Path to expression matrix file (optional)|`None`
|`--metadata <PATH>` |Path to expression metadata file (optional)|`None`
|`--output_folder <STR>` |Name of the output folder (optional)|`CYP_Annotator_Output`
|`--processed_input_folder <STR>` |Name of the processed input folder (optional)|`None`
|`--name <STR>` |STRING\_USED\_AS\_PREFIX\_IN\_FILENAMES|`""`
|`--trim_names <y/n>` |Trim sequence names at first space or tab (y/n)|`y`

#### Optional arguments for tool adjustments

|Command|Description|Default
|--|--|--
|`--mode_aln <STR>` |Tool used for multiple alignments (`mafft/muscle`)|`mafft`
|`--mode_tree <STR>` |Tool used for tree construction (fasttree/raxml/iqtree)|`fasttree`
|`--mafft <STR>` |MAFFT command|`mafft`
|`--muscle <STR>` |MUSCLE command|`muscle`
|`--fasttree <STR>` |Fasttree command|`fasttree`
|`--raxml <STR>` |RAXML command|`raxml-ng`
|`--iqtree <STR>` |IQ-TREE command|`iqtree`
|`--blastp <STR>` |PATH\_TO\_AND\_INCLUDING\_BINARY\_BLASTp|`blastp`
|`--makeblastdb <STR>` |PATH\_TO\_AND\_INCLUDING\_BINARY\_MAKEBLASTDB|`makeblastdb`
|`--hmmsearch <STR>` |PATH\_TO\_HMMSEARCH|`hmmsearch`

#### Optional arguments for search and classification

|Command|Description|Default
|--|--|--
|`--use_hmmer <y/n>` |Use HMMER for initial candidate search (y/n)|`n`
|`--simcutp <FLOAT>` |BLASTP\_SIMILARITY\_CUTOFF|`40.00`
|`--poscutp <INT>` |BLASTP\_POSSIBLE\_HIT\_NUMBER\_PER\_BAIT\_CUTOFF|`100`
|`--lencutp <INT>` |BLASTP\_MIN\_LENGTH\_CUTOFF|`200`
|`--bitcutp <INT>` |BLASTP\_BITSCORE\_CUTOFF|`80`
|`--filterdomain <y/n>` |DOMAIN\_FILTER\_FOR\_CLASSIFICATION (y/n)|`n`
|`--minscore <FLOAT>` |MINIMAL\_SCORE to be considered ingroup|`0.5`
|`--numneighbours <INT>` |NUMBER\_OF\_NEIGHBOURS\_FOR\_CLASSIFICATION|`24`
|`--neighbourdist <FLOAT>` |NEIGHBOUR\_DISTANCE|`5.0`
|`--minneighbours <INT>` |MINIMAL\_NUMBER\_OF\_NEIGHBOURS|`0`

#### Optional arguments for orthologs and paralogs

|Command|Description|Default
|--|--|--
|`--static_pd <y/n>` |Ortholog assignment with static thresholds (y/n)|`n`
|`--threshold_factor <FLOAT>` |Factor for adding deviation/IQR to mean/median in dynamic threshold calculation|`0.5`
|`--subfamily_threshold <FLOAT>` |Theshold for patristic distance considering orthologs|`1.1`
|`--family_threshold <FLOAT>` |Theshold for patristic considering further neighbours|`2.7`
|`--individual_tree <y/n>` |Create individual tree with specific bait sequences (y/n)|`n`
|`--bait_column <STR>` |Optional: column name for bait filtering|`Evidence`
|`--bait_keyword <STR>` |Keyword for bait filtering|`Literature`
|`--ortholog_prefix <STR>` |Prefix for ortholog filtering|`All`
|`--individual_ortholog_prefix <STR>`|Prefix for individual ortholog filtering|`None`
|`--collapse <y/n>` |Reduce in-paralogs to one representative|`y`
|`--paralogdist <FLOAT>` |Distance of paralogs in masking step|`10.0`

#### Optional arguments for functional analysis (Domains, Motifs, Expression)

|Command|Description|Default
|--|--|--
|`--domain_Score <FLOAT>` |c-Evalue for hmm domain integration|`100`
|`--motif_cEvalue <FLOAT>` |c-Evalue for hmm motif integration|`0.01`
|`--min_avg_tpm <FLOAT>` |Average tpm for genes to be considered expressed|`1.0`
|`--min_single_tpm <FLOAT>` |Single tpm for genes to be considered expressed|`5.0`
|`--min_paralog_tpm <FLOAT>` |Min tpm for paralog conservation|`1.0`

#### Optional arguments regarding performance

|Command|Description|Default
|--|--|--
|`--cpu_max <INT>` |Max CPUs|`4`
|`--parallel <y/n>` |Run classification in parallel mode (y/n)|`y`
|`--num_process_candidates <INT>`|Max number of candidates per ingroup/outgroup classification|`200`

-----

### Adjustment of input data files

The pipeline relies on several key input files. You can provide them using one of the three methods described in the **Usage** section.

|File|Argument|Description|Format
|--|--|--|--
|Baits|`--baits`|FASTA file containing the reference (bait) sequences, including ingroup (CYPs) and outgroup sequences.|`.fasta` / `.fa`
|Baits Info|`--baits_info`|**Crucial CSV file.** Must contain a header. The first column `ID` must match the FASTA headers in the baits file. Other columns, like `Family`, `Subfamily`, and `Evidence`, are required for annotation and filtering.|`.csv`
|Subject|`--subject`|FASTA file(s) or folder(s) containing the sequences (genome, transcriptome) to be annotated. Can be CDS or PEP.|`.fasta` / `.fa`
|Data Folder|`--data`|A single folder containing `baits.fasta`, `baits_info.csv`, and `subject.fasta` (or a `subjects/` subdirectory). Optional files like `hmm_domains.hmm` can also be placed here.|Folder
|File Collection|`--file_collection`|A CSV file specifying paths to all other inputs. Overrides all other path arguments if used.|`.csv`

## Requirements

Python3, dendropy, pandas, matplotlib, seaborn, ete3, BLAST+, HMMER, MAFFT or MUSCLE, FastTree or RAxML-NG or IQ-TREE

## Reference

> Georgi K. 'Development of a bioinformatics tool for automatic, functional annotation of plant cytochromes P450.'

Functions used in this script are in part taken from:

  * **MYB\_annotator** ([doi: 10.1186/s12864-022-08452-5](http://dx.doi.org/10.1186/s12864-022-08452-5))
  * **bHLH\_annotator** ([doi: 10.1186/s12864-023-09877-2](https://doi.org/10.1186/s12864-023-09877-2))
