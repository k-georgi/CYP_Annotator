# region Information
# Skripte aus bHLH-Annotator, MYB-Annotator und Tree-Annotator werden genutzt (Referenzen ergänzen)
"""
Region containing general information about the script.
"""

""""
Functions used in this script are in part taken from bHLH annotator and MYB annotator
"""

__version__ = "v0.1"
__usage__ = """
                    python3 CYP_Annotator.py
                    --data <PATH_TO_DATA_FOLDER> 
                    | --baits <PATH_TO_BAITS> --subject <PATH_TO_SUBJECT> --baits_info <PATH_TO_BAITS_INFO>
                    | --file_collection <PATH_TO_FILE_COLLECTION_CSV>

                    optional:
                    --outgroup <PATH_TO_OUTGROUP_FASTA>
                    --hmm <PATH_TO_BAIT_HMM_FILE>
                    --hmm_domains <PATH_TO_HMM_DOMAINS_FILE>
                    --hmm_motifs <PATH_TO_HMM_MOTIFS_FILE>
                    --protein_motifs <PATH_TO_PROTEIN_MOTIFS_FILE>
                    --expression <EXPRESSION_MATRIX_CSV>
                    --metadata <METADATA_CSV_FOR_EXPRESSION>

                    --output_folder <OUTPUT_FOLDER>
                    --processed_input_folder <PROCESSED_INPUT_FOLDER>
                    --name <STRING_USED_AS_PREFIX_IN_FILENAMES>
                    --trim_names <TRIM_SEQUENCE_NAMES>(y/n)[n]

                    --cpu_max <MAX_NUMBER_OF_CPUS>[4]
                    --parallel <PARALLEL_CLASSIFICATION_MODE>(y/n)[y]
                    --num_process_candidates <CANDIDATES_PER_PARALLEL_RUN>[200]

                    --use_hmmer <USE_HMMER_FOR_INITIAL_SEARCH>(y/n)[n]

                    --mode_aln <ALIGNMENT_TOOL>(mafft|muscle)[mafft]
                    --mode_tree <TREE_TOOL>(fasttree|raxml)[fasttree]
                    --mafft <PATH_TO_MAFFT>[mafft]
                    --muscle <PATH_TO_MUSCLE>[muscle]
                    --fasttree <PATH_TO_FASTTREE>[fasttree]
                    --raxml <PATH_TO_RAXML>[raxml-ng]

                    --blastp <PATH_TO_BLASTP>[blastp]
                    --makeblastdb <PATH_TO_MAKEBLASTDB>[makeblastdb]
                    --hmmsearch <PATH_TO_HMMSEARCH>[hmmsearch]

                    --simcutp <BLASTP_SIMILARITY_CUTOFF>[40.0]
                    --poscutp <MAX_HITS_PER_BAIT>[100]
                    --lencutp <MIN_HIT_LENGTH>[80]
                    --bitcutp <MIN_BLASTP_BITSCORE>[60]

                    --filterdomain <DOMAIN_FILTER_FOR_CLASSIFICATION>(y/n)[y]
                    --minscore <MIN_SCORE_FOR_INGROUP>[0.5]
                    --numneighbours <NUMBER_OF_NEIGHBOURS>[20]
                    --neighbourdist <NEIGHBOUR_DISTANCE>[5.0]
                    --minneighbours <MIN_NEIGHBOURS>[0]

                    --static_pd <USE STATIC PATRISTIC DISTANCE THRESHOLDS>(y/n)[n]
                    --threshold_factor <FACTOR FOR DYNAMIC THRESHOLD CALCULATION>[0.2]
                    --subfamily_threshold <ORTHOLOG_DISTANCE_THRESHOLD>[0.5]
                    --family_threshold <NEIGHBOUR_DISTANCE_THRESHOLD>[2.0]
                    --individual_tree <GENERATE_TREE_FOR_EACH_BAIT>(y/n)[y]

                    --bait_column <COLUMN_FOR_BAIT_FILTERING>[Evidence]
                    --bait_keyword <KEYWORD_FOR_BAIT_FILTERING>[Literature]
                    --ortholog_prefix <PREFIX_FOR_ORTHOLOG_OUTPUT>[selected]
                    --individual_ortholog_prefix <PREFIX_FOR_INDIVIDUAL_ORTHOLOGS>[None]

                    --collapse <REDUCE_PARALOGS_TO_REPRESENTATIVES>(y/n)[y]
                    --paralogdist <DISTANCE_THRESHOLD_FOR_PARALOGS>[10.0]
                    --min_paralog_tpm <MIN_TPM_FOR_PARALOGS>[1.0]

                    --domain_Score <HMM_DOMAIN_C-EVALUE>[100]
                    --motif_cEvalue <MOTIF_C-EVALUE>[0.01]

                    --min_avg_tpm <MAX_AVG_TPM_FOR_PSEUDOGENES>[1.0]
                    --min_single_tpm <MAX_SINGLE_TPM_FOR_PSEUDOGENES>[5.0]

                    bug reports and feature requests: kgeo@uni-bonn.de
                    """

#endregion

# region Import
"""
Region containing all module imports and dependency checks.
"""

import argparse
import csv
import dendropy
import datetime
import importlib.util
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pandas as pd
import random
import re
import seaborn as sns
import shutil
import subprocess
import sys
import time

from collections import defaultdict
from operator import itemgetter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

try:
	import hashlib
except ImportError:
	pass

Tree, TreeStyle, NodeStyle, TextFace = None, None, None, None
ete3_available = False

if sys.version_info < (3, 13): # cgi module removed in Python 3.13
    try:
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace
        ete3_available = True # Set flag to True only if import succeeds
    except ImportError:
        print("Warning: ete3 library found but could not be imported. Tree plotting with ete3 disabled.")

else:
    print(f"Warning: Python version {sys.version_info.major}.{sys.version_info.minor} is >= 3.13. "
          "The 'cgi' module required by ete3's web components is removed. "
          "Skipping ete3 import to avoid errors. Tree plotting with ete3 disabled.")
# endregion

# region Constants
"""
Region containing relevant constants.
"""

FORBIDDEN_CHARS = [";", ":", "(", ")", "="]
FASTA_EXTENSIONS = {'.fasta', '.fa'}

GENETIC_CODE = {
    'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
    'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
    'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
    'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
    'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
    'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
    'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
    'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
    'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
    'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
    'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
    'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
    'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
}

BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4
}
#endregion

# region Input
"""
Region handling command-line argument parsing and validation for sequence analysis.
Includes input file paths, tool configurations, and runtime parameters.
"""

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Comparison of sequence groups in FASTA files. "
    "Either provide: (1) a Data folder OR (2) both --baits and --subject."
)

# Data folder (optional, but required if baits/subject not provided)
parser.add_argument(
    "--data",
    type=str,
    default=None,
    help="Path to the Data folder (required if --baits/--subject not provided)"
)

# File Collection (optional, but required if neither baits/subject nor data provided)
parser.add_argument(
    "--file_collection",
    type=str,
    default=None,
    help="Path to the file collection csv file(required if --baits/--subject or data not provided)"
)

# Baits and Subject (optional, but both needed if data is missing)
parser.add_argument(
    "--baits",
    type=str,
    nargs="+",  # Accepts one or more paths
    default=None,
    help="Path to baits (FASTA file or folder). Required if --data is not provided."
)
parser.add_argument(
    "--baits_info",
    type=str,
    default=None,
    help="Path to baits info file (optional)"
)
parser.add_argument(
    "--subject",
    type=str,
    nargs="+",  # Accepts one or more paths
    default=None,
    help="Path to subject (FASTA file(s) or folder(s)). Required if --data is not provided."
)

# Optional files
parser.add_argument(
    "--outgroup",
    type=str,
    default=None,
    help="Path to outgroup FASTA file (optional)"
)
parser.add_argument(
    "--hmm",
    type=str,
    default=None,
    help="Path to bait HMM file"
)
parser.add_argument(
    "--hmm_domains",
    type=str,
    default=None,
    help="Path to hmm domains file (optional)"
)
parser.add_argument(
    "--hmm_motifs",
    type=str,
    default=None,
    help="Path to hmm motifs file (optional)"
)
parser.add_argument(
    "--protein_motifs",
    type=str,
    default=None,
    help="Path to protein motifs file (optional)"
)

# Optional Arguments
parser.add_argument(
    "--trim_names",
    type=str,
    default="y",
    choices=['y', 'n', 'yes', 'no'],
    help="Trim sequence names at first space or tab (y/n, default: y)"
)
parser.add_argument(
    "--name",
    type=str,
    default="",
    help="STRING_USED_AS_PREFIX_IN_FILENAMES"
)
parser.add_argument(
    "--cpu_max",
    type=int,
    default=4,
    help="Max CPUs (default: 4)"
)

# Folders
parser.add_argument(
    "--output_folder",
    type=str,
    default="CYP_Annotator_Output",
    help="Name of the output folder (optional)"
)
parser.add_argument(
    "--processed_input_folder",
    type=str,
    default=None,
    help="Name of the processed input folder (optional)"
)

# Search
parser.add_argument(
    "--use_hmmer",
    type=str,
    default="n",
    choices=['y', 'n', 'yes', 'no'],
    help="Use HMMER for initial candidate search (y/n, default: n)"
)
parser.add_argument(
    "--blastp",
    type=str,
    default="blastp",
    help="PATH_TO_AND_INCLUDING_BINARY_BLASTp"
)
parser.add_argument(
    "--makeblastdb",
    type=str,
    default="makeblastdb",
    help="PATH_TO_AND_INCLUDING_BINARY_MAKEBLASTDB"
)
parser.add_argument(
    "--hmmsearch",
    type=str,
    default="hmmsearch",
    help="PATH_TO_HMMSEARCH"
)
parser.add_argument(
    "--simcutp",
    type=float,
    default=40.00,
    help="BLASTP_SIMILARITY_CUTOFF (default = 40.00)"
)
parser.add_argument(
    "--poscutp",
    type=int,
    default=100,
    help="BLASTP_POSSIBLE_HIT_NUMBER_PER_BAIT_CUTOFF (default = 100)"
)
parser.add_argument(
    "--lencutp",
    type=int,
    default=200,
    help="BLASTP_MIN_LENGTH_CUTOFF (default = 80)"
)
parser.add_argument(
    "--bitcutp",
    type=int,
    default=80,
    help="BLASTP_BITSCORE_CUTOFF (default = 60)"
)

# Classification
parser.add_argument(
    "--parallel",
    type=str,
    default="y",
    choices=['y', 'n', 'yes', 'no'],
    help="Run classification in parallel mode (y/n, default: y)"
)
parser.add_argument(
    "--num_process_candidates",
    type=int,
    default=200,
    help="Max number of candidates per ingroup/outgroup classification (default: 200)"
)
parser.add_argument(
    "--mode_aln",
    type=str,
    default="mafft",
    choices=['mafft', 'muscle'],
    help="Tool used for multiple alignments (mafft/muscle, default: mafft)"
)
parser.add_argument(
    "--mode_tree",
    type=str,
    default="fasttree",
    choices=['fasttree', 'raxml', 'iqtree'],
    help="Tool used for tree construction (fasttree/raxml/iqtree, default: fasttree)"
)
parser.add_argument(
    "--mafft",
    type=str,
    default="mafft",
    help="MAFFT command"
)
parser.add_argument(
    "--muscle",
    type=str,
    default="muscle",
    help="MUSCLE command"
)
parser.add_argument(
    "--fasttree",
    type=str,
    default="fasttree",
    help="Fasttree command"
)
parser.add_argument(
    "--raxml",
    type=str,
    default="raxml-ng",
    help="RAXML command"
)
parser.add_argument(
    "--iqtree",
    type=str,
    default="iqtree",
    help="IQ-TREE command"
)
parser.add_argument(
    "--numneighbours",
    type=int,
    default=24,
    help="NUMBER_OF_NEIGHBOURS_FOR_CLASSIFICATION (default: 10)"
)
parser.add_argument(
    "--neighbourdist",
    type=float,
    default=5.0,
    help="NEIGHBOUR_DISTANCE (default: 5.0)"
)
parser.add_argument(
    "--minneighbours",
    type=int,
    default=0,
    help="MINIMAL_NUMBER_OF_NEIGHBOURS (default: 0)"
)
parser.add_argument(
    "--minscore",
    type=float,
    default=0.5,
    help="MINIMAL_SCORE to be considered ingroup (default: 0.5)"
)
parser.add_argument(
    "--filterdomain",
    type=str,
    default="n",
    choices=['y', 'n', 'yes', 'no'],
    help="DOMAIN_FILTER_FOR_CLASSIFICATION (y/n, default: y)"
)

# Orthologs
parser.add_argument(
    "--static_pd",
    type=str,
    default="n",
    choices=['y', 'n', 'yes', 'no'],
    help="Ortholog assignment with static thresholds (y/n, default: n)"
)
parser.add_argument(
    "--threshold_factor",
    type=float,
    default=0.5,
    help="Factor for adding deviation/IQR to mean/median in dynamic threshold calculation (default: 0.2)"
)
parser.add_argument(
    "--subfamily_threshold",
    type=float,
    default=1.1,
    help="Theshold for patristic distance considering orthologs (default: 1.1)"
)
parser.add_argument(
    "--family_threshold",
    type=float,
    default=2.7,
    help="Theshold for patristic considering further neighbours (default: 2.2)"
)
parser.add_argument(
    "--individual_tree",
    type=str,
    default="n",
    choices=['y', 'n', 'yes', 'no'],
    help="Create individual tree with specific bait sequences (y/n, default: y)"
)
parser.add_argument(
    '--bait_column',
    type=str,
    default='Evidence',
    help='Optional: column name for bait filtering (default: Species)')
parser.add_argument(
    '--bait_keyword',
    type=str,
    default='Literature',
    help='Keyword for bait filtering (default: Arabidopsis thaliana)')
parser.add_argument(
    '--ortholog_prefix',
    type=str,
    default='All',
    help='Prefix for ortholog filtering (default: selected)')
parser.add_argument(
    '--individual_ortholog_prefix',
    type=str,
    default=None,
    help='Prefix for individual ortholog filtering (default: Arabidopsis)')

# Domain Check
parser.add_argument(
    "--domain_Score",
    type=float,
    default=100,
    help="c-Evalue for hmm motif integration (default: 100)"
)

# Motif Check
parser.add_argument(
    "--motif_cEvalue",
    type=float,
    default=0.01,
    help="c-Evalue for hmm motif integration (default: 0.01)"
)

# Expression
parser.add_argument(
    "--expression",
    type=str,
    default=None,
    help="Path to expression matrix file (optional)"
)
parser.add_argument(
    "--metadata",
    type=str,
    default=None,
    help="Path to expression metadata file (optional)"
)
parser.add_argument(
    "--min_avg_tpm",
    type=float,
    default=1.0,
    help="Average tpm that pseudogenes may not exceed (default: 1.0)"
)
parser.add_argument(
    "--min_single_tpm",
    type=float,
    default=5.0,
    help="Single tpm that pseudogenes may not exceed (default: 5.0)"
)

# Paralogs
parser.add_argument(
    "--collapse",
    type=str,
    default="y",
    choices=['y', 'n', 'yes', 'no'],
    help="Reduce in-paralogs to one representative"
)
parser.add_argument(
    "--paralogdist",
    type=float,
    default=10.0,
    help="Distance of paralogs in masking step (default: 10.0)"
)
parser.add_argument(
    "--min_paralog_tpm",
    type=float,
    default=1.0,
    help="Min tpm for paralog conservation (default: 1.0)"
)

def validate_args(args):
    """Validates inputs with XOR logic: Either data OR file_collection OR (baits + baits_info)"""
    has_data = args.data is not None
    has_file_collection = args.file_collection is not None
    has_baits_info = args.baits is not None and args.baits_info is not None
    
    # Count how many input methods are provided
    input_methods = sum([has_data, has_file_collection, has_baits_info])
    
    # Check that exactly one input method is selected
    if input_methods != 1:
        raise ValueError(
            "Exactly one of the following must be provided:\n"
            "1. --data folder\n"
            "2. --file_collection CSV file\n"
            "3. Both --baits and --baits_info\n"
            f"Current inputs: data={args.data}, file_collection={args.file_collection}, "
            f"baits={args.baits}, baits_info={args.baits_info}"
        )

    # Validate based on selected input method
    if has_data:
        if not os.path.isdir(args.data):
            raise ValueError(f"Data folder does not exist: {args.data}")
            
    elif has_file_collection:
        if not os.path.isfile(args.file_collection):
            raise ValueError(f"File collection CSV does not exist: {args.file_collection}")
        if not args.file_collection.lower().endswith('.csv'):
            raise ValueError("File collection must be a CSV file")
            
    else:  # has_baits_info
        # Validate baits
        if not os.path.isfile(args.baits):
            raise ValueError(f"Baits path must be a file: {args.baits}")
        if not args.baits.lower().endswith(('.fasta', '.fa')):
            raise ValueError("Baits file must be a FASTA file (.fasta/.fa)")
            
        # Validate baits_info
        if not os.path.isfile(args.baits_info):
            raise ValueError(f"Baits info file does not exist: {args.baits_info}")
        if not args.baits_info.lower().endswith('.csv'):
            raise ValueError("Baits info file must be a CSV file")

    return args
# endregion

# region Preparation
"""
Region including utility functions for checking tool availability, collecting input files,
preparing and processing FASTA/motif files, managing file paths, and generating metadata documentation.
"""

def prepare_input_files(args: 'argparse.Namespace') -> Dict[str, any]:
    """
    Finds all raw input files, processes them (clean, translate, copy),
    and returns a dictionary with paths to the processed files.
    This function combines the logic of collect_files and process_all_files.
    """
    # 1. Setup and Path Initialization
    output_dir = args.processed_input_folder or os.path.join(args.output_folder, "Processed_Input/")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    raw_paths = {}
    processed_paths = {}

    # 2. Find Raw File Paths (from data folder and explicit args)
    if args.data:
        data_path = Path(args.data)
        if not data_path.is_dir():
            raise ValueError(f"Data folder does not exist: {args.data}")
        
        subjects_dir = data_path / "subjects"
        raw_paths['subject_files'] = collect_fasta_files([str(p) for p in [data_path / "subject.fasta", subjects_dir] if p.exists()])
        raw_paths['bait_file'] = str(data_path / "baits.fasta") if (data_path / "baits.fasta").exists() else None
        raw_paths['baits_info_file'] = str(data_path / "baits_info.csv") if (data_path / "baits_info.csv").exists() else None
        raw_paths['protein_motifs_file'] = find_case_insensitive_file(data_path, "protein_motifs.txt")
        raw_paths['hmm_motifs_file'] = find_case_insensitive_file(data_path, "hmm_motifs.hmm")
        raw_paths['hmm_domains_file'] = find_case_insensitive_file(data_path, "hmm_domains.hmm")
        raw_paths['expression_file'] = find_case_insensitive_file(data_path, "expression.txt")
        raw_paths['metadata_file'] = find_case_insensitive_file(data_path, "metadata.csv")

    # Explicit paths override data folder paths
    if args.subject:
        raw_paths['subject_files'] = collect_fasta_files([args.subject])
    if args.baits:
        raw_paths['bait_file'] = args.baits
    if args.baits_info:
        raw_paths['baits_info_file'] = args.baits_info
    if args.protein_motifs:
        raw_paths['protein_motifs_file'] = args.protein_motifs
    if args.hmm_motifs:
        raw_paths['hmm_motifs_file'] = args.hmm_motifs
    if args.hmm_domains:
        raw_paths['hmm_domains_file'] = args.hmm_domains
    if args.expression:
        raw_paths['expression_file'] = args.expression
    if args.metadata:
        raw_paths['metadata_file'] = args.metadata

    # 3. Validate Required Files
    if not raw_paths.get('subject_files'):
        raise ValueError("No valid subject files found")
    if not raw_paths.get('bait_file') or not Path(raw_paths['bait_file']).is_file():
        raise ValueError("No valid bait file found")
    if not raw_paths.get('baits_info_file') or not Path(raw_paths['baits_info_file']).is_file():
        raise ValueError("No valid baits info file found")

    # 4. Process Files
    # Process Baits and create subgroups
    baits_info, info_headers = read_baits_info(raw_paths['baits_info_file'])
    processed_paths['baits_info'] = baits_info
    processed_paths['baits_info_headers'] = info_headers
    
    bait_results = process_baits_file(
        raw_paths['bait_file'],
        baits_info,
        output_dir,
        args.trim_names
    )
    processed_paths.update({
        'baits_with_outgroup_path': bait_results[0],
        'baits_no_outgroup_path': bait_results[1],
        'transcript_baits_path': bait_results[2],
        'protein_baits_path': bait_results[3],
        'literature_baits_path': bait_results[4],
        'outgroup_path': bait_results[5]
    })

    # Process Subjects
    processed_paths['subject_files_paths'] = [
        process_single_file(sf, Path(sf).stem, output_dir, args.trim_names) 
        for sf in raw_paths['subject_files']
    ]

    # Process Optional Files
    if raw_paths.get('protein_motifs_file'):
        processed_paths['protein_motifs_path'] = copy_file(raw_paths['protein_motifs_file'], "protein_motifs", output_dir)
    if raw_paths.get('hmm_motifs_file'):
        processed_paths['hmm_motifs_path'] = copy_file(raw_paths['hmm_motifs_file'], "hmm_motifs", output_dir)
    if raw_paths.get('hmm_domains_file'):
        processed_paths['hmm_domains_path'] = copy_file(raw_paths['hmm_domains_file'], "hmm_domains", output_dir)
    if raw_paths.get('baits_info_file'):
        processed_paths['baits_info_path'] = copy_file(raw_paths['baits_info_file'], "baits_info", output_dir)
    if raw_paths.get('expression_file'):
        processed_paths['expression_path'] = process_expression_file(raw_paths['expression_file'], "expression", output_dir, args.trim_names)
    if raw_paths.get('metadata_file'):
        processed_paths['metadata_path'] = copy_file(raw_paths['metadata_file'], "metadata", output_dir)

    return processed_paths

def check_required_tools() -> Dict[str, bool]:
    """Check installation of required bioinformatics tools."""
    tools = {
        'MAFFT': bool(shutil.which('mafft')),
        'MUSCLE': bool(shutil.which('muscle')),
        'FASTTREE': bool(shutil.which('fasttree')),
        'BLAST': bool(shutil.which('blastp')),
        'HMMER': bool(shutil.which('hmmsearch')),
        'RAXML': bool(shutil.which('raxmlHPC')) or bool(shutil.which('raxml-ng'))
    }

    errors = []
    if not (tools['MAFFT'] or tools['MUSCLE']):
        errors.append("ERROR: Neither MAFFT nor MUSCLE is installed. At least one alignment tool is required.")
    if not (tools['FASTTREE'] or tools['RAXML']):
        errors.append("ERROR: Neither FASTTREE nor RAxML is installed. At least one phylogeny tool is required.")
    if not (tools['BLAST'] or tools['HMMER']):
        errors.append("ERROR: Neither BLAST nor HMMER is installed. At least one search tool is required.")

    if errors:
        print("\n".join(errors))
        print("\nInstallation guide:")
        print("- MAFFT: http://mafft.cbrc.jp/alignment/software/")
        print("- MUSCLE: https://www.drive5.com/muscle/")
        print("- FASTTREE: https://anaconda.org/bioconda/fasttree")
        print("- RAxML: https://cme.h-its.org/exelixis/web/software/raxml/")
        print("- BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download")
        print("- HMMER: http://hmmer.org/")
        sys.exit(1)

    return tools

def process_baits_file(
    bait_file: str,
    baits_info: Dict[str, Dict[str, str]],
    output_dir: str = "Processed_Input",
    trim_names: str = 'n'
) -> Tuple[str, str, str, str, str, str]:
    """
    Process bait file and create group-specific files based on pre-read baits_info.
    Returns tuple of paths: (all_baits, no_outgroup, transcript, protein, literature, outgroup)
    """
    output_dir_path = Path(output_dir)
    all_baits_path = process_single_file(bait_file, "baits_with_outgroup", output_dir, trim_names)
    sequences = read_fasta(all_baits_path)

    bait_groups = {seq_id: subdict.get("Evidence", "-") for seq_id, subdict in baits_info.items()}

    output_files = {
        'no_outgroup': output_dir_path / "clean_baits.fasta", 'transcript': output_dir_path / "transcript_baits.fasta",
        'protein': output_dir_path / "protein_baits.fasta", 'literature': output_dir_path / "literature_baits.fasta",
        'outgroup': output_dir_path / "outgroup.fasta"
    }

    if all(f.exists() for f in output_files.values()):
        print("Group-specific bait files already exist, skipping creation.")
        return tuple([str(p) for p in [all_baits_path] + list(output_files.values())])

    grouped_sequences = {key: {} for key in output_files}
    grouped_sequences['no_outgroup'] = {}

    for seq_id, seq in sequences.items():
        group = bait_groups.get(seq_id, "").lower()
        
        if "outgroup" in group:
            grouped_sequences['outgroup'][seq_id] = seq
        else:
            grouped_sequences['no_outgroup'][seq_id] = seq
            if "transcript" in group: grouped_sequences['transcript'][seq_id] = seq
            if "protein" in group: grouped_sequences['protein'][seq_id] = seq
            if "literature" in group: grouped_sequences['literature'][seq_id] = seq
            
    for group_name, seq_dict in grouped_sequences.items():
        if seq_dict and group_name in output_files:
            write_fasta(seq_dict, output_files[group_name])

    return tuple([str(p) for p in [all_baits_path] + list(output_files.values())])

def read_baits_info(baits_info_file: str) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """Read baits_info CSV and return ID-to-attribute dictionary and header list."""
    baits_info, headers = {}, []
    try:
        with open(baits_info_file, 'r', newline='') as f:
            reader = csv.reader(f)
            headers = next(reader)
            for row in reader:
                if row:
                    seq_id = row[0].strip()
                    baits_info[seq_id] = {headers[i].strip(): row[i].strip() for i in range(1, len(row))}
    except Exception as e:
        raise ValueError(f"Invalid baits info file {baits_info_file}: {e}")
    return baits_info, headers

def process_single_file(
    input_file: str, output_name: str, output_dir: str = "Processed_Input", trim_names: str = 'n'
) -> str:
    """Process a single FASTA file: clean, translate, write output and mapping."""
    output_dir_path = Path(output_dir)
    output_fasta = output_dir_path / f"{output_name}.fasta"
    mapping_file = output_dir_path / f"{output_name}_mapping.txt"

    if output_fasta.exists() and mapping_file.exists():
        print(f"{input_file} already processed, skipping.")
        return str(output_fasta)

    sequences = load_and_process_sequences(input_file, trim_names)
    cleaned_sequences, sequence_mapping = {}, {}
    
    for original_id, sequence in sequences.items():
        cleaned_id = clean_sequence_id(original_id, trim_names)
        cleaned_sequences[cleaned_id] = sequence
        sequence_mapping[original_id] = cleaned_id

    write_fasta(cleaned_sequences, output_fasta)
    write_mapping_table(sequence_mapping, mapping_file)
    return str(output_fasta)

def load_and_process_sequences(file_path: str, trim_names: str = 'n') -> Dict[str, str]:
    """Load FASTA sequences and translate if nucleotide."""
    sequences = read_fasta(file_path)
    if sequences and is_nucleotide_sequence(next(iter(sequences.values()))):
        sequences = translate_sequences(sequences)
    if trim_names.upper() in ['Y', 'YES']:
        sequences = {clean_sequence_id(seq_id, trim_names): seq for seq_id, seq in sequences.items()}
    return sequences

def read_fasta(file_path: str) -> Dict[str, str]:
    """Read FASTA file and return sequences dictionary."""
    sequences = {}
    current_id = None
    current_seq = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line.upper())

            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)

        if not sequences:
            raise ValueError(f"FASTA file is empty: {file_path}")
        
        return sequences

    except Exception as e:
        raise ValueError(f"Error reading FASTA file {file_path}: {str(e)}")

def is_nucleotide_sequence(sequence: str) -> bool:
    """Determine if sequence is nucleotide or peptide."""
    clean_seq = re.sub(r'[\- ]', '', sequence[:100])
    return all(c in "ACGTN" for c in clean_seq)

def translate_sequences(sequences: Dict[str, str]) -> Dict[str, str]:
#     """Translate nucleotide sequences to peptide sequences."""
#     return {seq_id: str(seq(re.sub(r'[\- ]', '', seq)).translate()) for seq_id, seq in sequences.items()}

    translated_sequences = {}
        
    for seq_id, nuc_sequence in sequences.items():
        
        # 1. Sequenz säubern (wie im Originalcode)
        #    .upper() ist wichtig, da die codon_table Großbuchstaben erwartet.
        cleaned_seq = re.sub(r'[\- ]', '', nuc_sequence).upper()
        
        peptide_seq = []
        
        # 2. Sequenz in 3er-Schritten (Codons) durchlaufen
        #    Der range-Stop stellt sicher, dass wir nur vollständige Codons nehmen
        #    (z.B. bei einer Länge von 7 wird nur bis Index 6 (Pos 0, 1, 2, 3, 4, 5) gegangen)
        for i in range(0, len(cleaned_seq) - (len(cleaned_seq) % 3), 3):
            codon = cleaned_seq[i:i+3]
            
            # 3. Codon übersetzen
            #    .get(codon, 'X') gibt 'X' zurück, falls das Codon nicht
            #    gefunden wird (z.B. 'NNN' oder andere ambigue Basen).
            amino_acid = GENETIC_CODE.get(codon, 'X')
            
            # 4. Bei Stop-Codon die Translation beenden
            if amino_acid == '*':
                break
            
            peptide_seq.append(amino_acid)
        
        # Das Dictionary mit der übersetzten Sequenz füllen
        translated_sequences[seq_id] = "".join(peptide_seq)
        
    return translated_sequences

def clean_sequence_id(seq_id: str, trim_names: str = 'n') -> str:
    """Clean sequence ID by removing forbidden characters and optionally trimming."""
    if trim_names.upper() in ['Y', 'YES']:
        seq_id = seq_id.split()[0]
    seq_id = re.sub(r'[():,;]', '-', seq_id)
    return seq_id.encode("ascii", "ignore").decode()

def write_fasta(sequences: Dict[str, str], output_path: Path) -> None:
    """Write sequences to FASTA file."""
    if output_path.exists() and output_path.stat().st_size > 0: return
    with open(output_path, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")

def write_mapping_table(mapping: Dict[str, str], output_path: Path) -> None:
    """Write sequence ID mapping table."""
    if output_path.exists() and output_path.stat().st_size > 0: return
    with open(output_path, 'w') as f:
        f.write("InitialID\tCleanID\n")
        for original_id, cleaned_id in mapping.items():
            f.write(f"{original_id}\t{cleaned_id}\n")

def process_expression_file(input_file: str, output_name: str, output_dir: str, trim_names: str) -> str:
    """Process expression file, cleaning IDs in the first column."""
    output_path = Path(output_dir) / f"{output_name}.txt"
    if output_path.exists() and output_path.stat().st_size > 0: return str(output_path)
    try:
        with open(input_file, 'r') as infile, open(output_path, 'w') as outfile:
            outfile.write(infile.readline()) # header
            for line in infile:
                parts = line.strip().split('\t')
                if parts:
                    parts[0] = clean_sequence_id(parts[0], trim_names)
                    outfile.write('\t'.join(parts) + '\n')
        return str(output_path)
    except Exception as e:
        raise ValueError(f"Failed to process expression file {input_file}: {e}")

def copy_file(input_file: str, output_name: str, output_dir: str) -> str:
    """Copy file to output directory with new name."""
    output_dir_path = Path(output_dir)
    suffix = ".hmm" if "hmm" in output_name else ".txt"
    output_path = output_dir_path / f"{output_name}{suffix}"
    if output_path.exists() and output_path.stat().st_size > 0: return str(output_path)
    try:
        shutil.copy(input_file, output_path)
        return str(output_path)
    except Exception as e:
        raise ValueError(f"Failed to copy {input_file} to {output_path}: {e}")

def collect_fasta_files(paths: List[str]) -> List[str]:
    """Collect FASTA files from given paths (files or directories)."""
    files = []
    for path_str in paths:
        path = Path(path_str)
        if path.is_file():
            files.append(str(path))
        elif path.is_dir():
            files.extend([str(f) for f in path.glob("*.fa*")])
    return [f for f in files if Path(f).suffix.lower() in {'.fa', '.fasta', '.fna', '.faa'}]

def find_case_insensitive_file(directory: Path, filename: str) -> Optional[str]:
    """Find file with case-insensitive matching."""
    for f in directory.iterdir():
        if f.is_file() and f.name.lower() == filename.lower():
            return str(f)
    return None

def read_file_collection_csv(csv_path: str) -> Dict[str, str]:
    """
    Read file collection CSV and return paths as dictionary.
    CSV format: First row headers, second row paths.
    Returns dictionary with keys: Data, Baits, Baits_Info, etc.
    """
    paths = {}
    try:
        with open(csv_path, 'r', newline='') as f:
            reader = csv.reader(f)
            headers = next(reader)
            values = next(reader)
            
            expected_columns = [
                'Data', 'Baits', 'Baits_Info', 'HMM_Domains', 
                'HMM_Motifs', 'Protein_Motifs', 'Expression', 'Metadata'
            ]
            
            if not all(col in headers for col in expected_columns):
                raise ValueError("CSV file is missing required columns")
            
            for header, value in zip(headers, values):
                if header in expected_columns and value.strip():
                    path = Path(value.strip())
                    if not path.is_absolute():
                        path = Path.cwd() / path
                    paths[header] = str(path.resolve())
    
    except Exception as e:
        raise ValueError(f"Error reading file collection CSV {csv_path}: {str(e)}")
    
    return paths

def load_name_mapping_table(mapping_table_file: str) -> Dict[str, str]:
    """
    Loads subject name mapping table from a tab-delimited file.
    Returns dictionary mapping clean IDs to original IDs.
    """
    mapping_table = {}
    try:
        with open(mapping_table_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    # Format: original_id clean_id
                    mapping_table[parts[1]] = parts[0]
    except Exception as e:
        raise ValueError(f"Failed to load mapping table {mapping_table_file}: {str(e)}")
    
    return mapping_table
    
def md5_calculator( input_file ):
    """Calculates MD5 checksum of a file."""
    
    with open( input_file, "rb" ) as f:
        content = f.read()
    try:
        return hashlib.md5( content ).hexdigest()
    except NameError:
        return "n/a"

def generate_documentation_file(
    doc_file, fam_bait_seq_file_all, fam_bait_seq_file, fam_info_file,
    hmm_domains_path, hmm_motifs_path, protein_motifs_path, output_folder, subject_files,
    search, mode_aln, mode_tree, blastp, makeblastdb, hmmsearch, cpu_max,
    mafft, muscle, raxml, fasttree,
    bitscore_cutoff_p, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p,
    min_score_cutoff, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff,
    dist_cutoff_factorB, parallel_mode, num_process_candidates,
    name, keepnames, collapse):
    """Generate TSV file documenting all used files, parameters, and tool versions."""

    with open(doc_file, "w") as out:
        # Column header
        out.write("Description\tInfo\tMD5\n")

        # Citation
        out.write("Citation\tPlease cite 'Georgi K. and Pucker B. (2025). 'Development of a bioinformatics tool for automatic, functional annotation of plant cytochromes P450.' when using CYP_Annotator.py.\tn/a\n")

        # Version
        out.write(f"CYP_Annotator.py version\t{__version__}\tn/a\n")

        # Files with MD5 checksums
        out.write(f"{name} bait file\t{fam_bait_seq_file_all}\t{md5_calculator(fam_bait_seq_file_all)}\n")
        out.write(f"{name} bait info file\t{fam_info_file}\t{md5_calculator(fam_info_file)}\n")
        out.write(f"Thinned bait file\t{fam_bait_seq_file}\t{md5_calculator(fam_bait_seq_file)}\n")

        for path in subject_files:
            out.write(f"Subject FASTA file\t{path}\t{md5_calculator(path)}\n")

        # Optional: HMM and Protein Motifs
        if hmm_domains_path:
            out.write(f"HMM domains file\t{hmm_domains_path}\t{md5_calculator(hmm_domains_path)}\n")
        else:
            out.write("HMM domains file\tn/a\tn/a\n")

        if hmm_motifs_path:
            out.write(f"HMM motifs file\t{hmm_motifs_path}\t{md5_calculator(hmm_motifs_path)}\n")
        else:
            out.write("HMM motifs file\tn/a\tn/a\n")

        if protein_motifs_path:
            out.write(f"Protein motifs file\t{protein_motifs_path}\t{md5_calculator(protein_motifs_path)}\n")
        else:
            out.write("Protein motifs file\tn/a\tn/a\n")

        # Output folder
        out.write(f"Output folder\t{output_folder}\tn/a\n")

        # Parameters
        out.write(f"Tool for initial candidate selection\t{search}\tn/a\n")
        out.write(f"Tool for alignment\t{mode_aln}\tn/a\n")
        out.write(f"Tool for tree construction\t{mode_tree}\tn/a\n")
        out.write(f"Maximal CPUs\t{cpu_max}\tn/a\n")
        out.write(f"Type of input\t{'CDS' if collapse else 'PEP'}\tn/a\n")
        out.write(f"Prefix of output file names\t{name}\tn/a\n")
        out.write(f"Keep sequence names (--keepnames)\t{keepnames}\tn/a\n")
        out.write(f"Collapse paralogs (--collapse)\t{collapse}\tn/a\n")

        # Tool paths
        out.write(f"blastp path\t{blastp}\tn/a\n")
        out.write(f"makeblastdb path\t{makeblastdb}\tn/a\n")
        out.write(f"hmmsearch path\t{hmmsearch}\tn/a\n")
        out.write(f"mafft path\t{mafft}\tn/a\n")
        out.write(f"muscle path\t{muscle}\tn/a\n")
        out.write(f"raxml path\t{raxml}\tn/a\n")
        out.write(f"fasttree path\t{fasttree}\tn/a\n")

        # BLAST Filters
        out.write(f"Minimal bitscore cutoff\t{bitscore_cutoff_p}\tn/a\n")
        out.write(f"Minimal BLASTp similarity cutoff\t{similarity_cutoff_p}\tn/a\n")
        out.write(f"Max number of BLASTp hits\t{possibility_cutoff_p}\tn/a\n")
        out.write(f"Minimal BLASTp hit length\t{length_cutoff_p}\tn/a\n")

        # Classification
        out.write(f"Minimal ingroup score\t{min_score_cutoff}\tn/a\n")
        out.write(f"Neighbourhood size\t{neighbour_cutoff}\tn/a\n")
        out.write(f"Branch length cutoff factor\t{mean_factor_cutoff}\tn/a\n")
        out.write(f"Minimal number of neighbours\t{min_neighbour_cutoff}\tn/a\n")
        out.write(f"Distance cutoff paralog masking\t{dist_cutoff_factorB}\tn/a\n")

        # Parallelization
        out.write(f"Parallelization enabled\t{parallel_mode}\tn/a\n")
        out.write(f"Max candidates per process\t{num_process_candidates}\tn/a\n")

        # Tool versions (only two columns!)
        try:
            muscle_version = subprocess.check_output(f"{muscle} -help", shell=True, stderr=subprocess.PIPE).decode("utf-8")[4:20]
            out.write(f"Muscle version\t{muscle_version}\n")
        except Exception:
            out.write("Muscle version detection failed\tn/a\n")

        try:
            mafft_version = subprocess.check_output(f"{mafft} --version", stderr=subprocess.STDOUT, shell=True).decode("utf-8").strip()
            out.write(f"MAFFT version\t{mafft_version}\n")
        except Exception:
            out.write("MAFFT version detection failed\tn/a\n")

        out.write("FastTree version\tPLEASE_ADD_MANUALLY\n")

        try:
            raxml_version = subprocess.check_output(f"{raxml} --version", shell=True, stderr=subprocess.PIPE).decode("utf-8").strip()
            out.write(f"RAxML version\t{raxml_version}\n")
        except Exception:
            out.write("RAxML version detection failed\tn/a\n")

        try:
            hmmsearch_version_raw = subprocess.check_output(f"{hmmsearch} -h", shell=True, stderr=subprocess.PIPE).decode("utf-8")
            hmmsearch_version = hmmsearch_version_raw.split("#")[2].strip()
            out.write(f"hmmsearch version\t{hmmsearch_version}\n")
        except Exception:
            out.write("hmmsearch version detection failed\tn/a\n")

#endregion

# region Search
"""
Region providing functions to run search für CYP candidates based on either BLAST or HMM.
"""

def run_blast_search(query_file, subject_file, output_file, blast_db_folder, cpu_max=1, makeblastdb="makeblastdb", blastp="blastp"):
    """Run BLAST search and return results file path"""
    blast_db = os.path.join(blast_db_folder, "blastdb")

    # Create BLAST database (silently)
    p = subprocess.Popen(
        args=[
            makeblastdb,
            "-in", subject_file,
            "-out", blast_db,
            "-dbtype", "prot"
        ],
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL
    )
    p.communicate()

    # Run BLAST search (silently)
    p = subprocess.Popen(
        args=[
            blastp,
            "-query", query_file,
            "-db", blast_db,
            "-out", output_file,
            "-outfmt", "6",
            "-evalue", "0.001",
            "-num_threads", str(cpu_max)
        ],
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL
    )
    p.communicate()

    return output_file
    
def run_hmmsearch(hmm_file: str, subject_file: str, output_file: str, hmmsearch: str = "hmmsearch", domtblout: bool = False) -> str:
    """Run HMMER search and return results file path"""
    
    if domtblout:
        args = [hmmsearch, "--domtblout", output_file, hmm_file, subject_file]
    else:
        args = [hmmsearch, "--tblout", output_file, hmm_file, subject_file]

    # Run silently by redirecting stdout and stderr
    p = subprocess.Popen(
        args=args,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    p.communicate()

    return output_file
    
def load_blast_results(blast_result_file, similarity_cutoff, possibility_cutoff, 
                        length_cutoff, bitscore):
    """Load and filter BLAST results"""
    valid_blast_hits = {}
    
    with open(blast_result_file, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            if (float(parts[2]) > similarity_cutoff and 
                float(parts[3]) > length_cutoff and 
                float(parts[-1]) >= bitscore):
                
                hit = {'gene': parts[0], 'score': float(parts[-1])}
                if parts[1] in valid_blast_hits:
                    valid_blast_hits[parts[1]].append(hit)
                else:
                    valid_blast_hits[parts[1]] = [hit]
    
    # Reduce hits to given number of possibilities
    final_hits = {}
    for key, hits in valid_blast_hits.items():
        sorted_hits = sorted(hits, key=lambda x: x['score'], reverse=True)
        genes = []
        for hit in sorted_hits:
            if hit['gene'] not in genes and len(genes) < possibility_cutoff:
                genes.append(hit['gene'])
        final_hits[key] = genes
    
    return final_hits

def load_hmmsearch_results(result_file):
    """Load HMMER search results"""
    results = {}
    with open(result_file, "r") as f:
        for line in f:
            if line[0] != "#":
                parts = line.strip().split()
                if parts:
                    results[parts[0]] = None
    return results

def hmm_build(fasta_file: str, hmm_output: str, alignment_tool: str = "mafft") -> str:
    """Build HMM from FASTA file with alignment preprocessing"""
    # Temporary files
    temp_dir = os.path.dirname(hmm_output)
    aln_file = os.path.join(temp_dir, "temp_alignment.fasta")
    trimmed_aln = os.path.join(temp_dir, "trimmed_alignment.fasta")
    
    # Create alignment
    try:
        # Simple alignment (no parallel processing needed)
        if alignment_tool == "mafft":
            with open(aln_file, 'w') as out_file:
                subprocess.run(
                    ["mafft", "--quiet", "--thread", "1", fasta_file],
                    check=True,
                    stdout=out_file,
                    stderr=subprocess.PIPE
                )
        else:  # MUSCLE
            subprocess.run(
                ["muscle", "-align", fasta_file, "-output", aln_file, "-threads", "1"],
                check=True,
                stderr=subprocess.PIPE
            )
        
        # Trim alignment
        alignment_trimming(aln_file, trimmed_aln)
        
        # Build HMM
        subprocess.run(
            ["hmmbuild", hmm_output, trimmed_aln],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # Clean up temporary files
        os.remove(aln_file)
        os.remove(trimmed_aln)
        
        return hmm_output
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"HMM build failed: {e.stderr.decode()}")
# endregion

# region Tree
"""
Region that handles alignment creation and tree construction.
Functions to process alignments and treefiles are provided as well.
"""

def create_alignments(num_par: int, tree_output_folder: str, name: str, 
                        number: int, aln_input_file: str, aln_file: str, 
                        mode_aln: str, mafft: str, muscle: str, cpu_par: int, 
                        cpu_max: int) -> None:
    """
    Creates alignments for all candidate chunks in parallel if chosen.
    """
    alignment_processes = []
    cpu_use = 0
    
    for i in range(num_par):
        aln_input = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_input_file}")
        aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_file}")
        
        while cpu_use + cpu_par > cpu_max:
            time.sleep(2)
            # Clean up finished processes and recalculate CPU usage
            active_processes = []
            for p in alignment_processes:
                if p.poll() is None: 
                    active_processes.append(p)
                else:
                    # Check for errors in completed processes
                    if p.returncode != 0:
                        print(f"Warning: Alignment process failed with return code {p.returncode}")
            alignment_processes = active_processes
            cpu_use = len(alignment_processes) * cpu_par
        
        if not os.path.isfile(aln):
            try:
                if mode_aln == "mafft":
                    cmd = [mafft, "--quiet", "--thread", str(cpu_par), aln_input]
                    with open(aln, 'w') as out_file:
                        p = subprocess.Popen(cmd, stdout=out_file, stderr=subprocess.PIPE)
                else:
                    cmd = [muscle, "-align", aln_input, "-output", aln, "-threads", str(cpu_par)]
                    p = subprocess.Popen(cmd, stderr=subprocess.PIPE)
                
                alignment_processes.append(p)
                cpu_use += cpu_par
                #print(f"Started alignment process {i+1}")
                
            except Exception as e:
                print(f"Error starting alignment {i+1}: {str(e)}")
    
    # Wait for all alignment processes to complete
    for p in alignment_processes:
        returncode = p.wait()
        if returncode != 0:
            print(f"Warning: Alignment process failed with return code {returncode}")
    
    #print("All alignment processes completed")

def alignment_trimming(aln_file: str, cln_aln_file: str, occupancy: float = 0.1) -> None:
    """
    Trims alignment by removing columns with low occupancy.
    """
    alignment = load_alignment(aln_file, {})
    
    if len(list(alignment.keys())) > 0:
        # identify valid residues in aligned sequences (columns with sufficient occupancy)
        valid_index = []
        for idx, aa in enumerate(list(alignment.values())[0]):
            counter = 0
            for key in list(alignment.keys()):
                if alignment[key][idx] != "-":
                    counter += 1
            if counter / float(len(list(alignment.keys()))) > occupancy:
                valid_index.append(idx)
            
        # generate new sequences
        with open(cln_aln_file, "w") as out:
            for key in list(alignment.keys()):
                seq = alignment[key]
                new_seq = []
                for idx in valid_index:
                    new_seq.append(seq[idx])
                out.write(">" + key + '\n' + "".join(new_seq) + '\n')
    else:
        with open(cln_aln_file, "w") as out:
            out.write("")

def load_alignment(aln_file: str, tmp_mapping: Dict[str, str]) -> Dict[str, str]:
    """
    Loads alignment from FASTA file and dict for header mapping.   
    Returns dictionary with sequence headers as keys and sequences as values.
    """
    sequences = {}
    
    with open(aln_file) as f:
        header = f.readline()[1:].strip()
        try:
            header = tmp_mapping[header]
        except KeyError:
            pass
        seq = []
        line = f.readline()
        while line:
            if line[0] == '>':
                sequences.update({header: "".join(seq)})
                header = line.strip()[1:]
                try:
                    header = tmp_mapping[header]
                except KeyError:
                    pass
                seq = []
            else:
                seq.append(line.strip())
            line = f.readline()
        sequences.update({header: "".join(seq)})
    return sequences

def parallel_tree_constructor(num_process_candidates: int, tree_output_folder: str, 
                            aln_candidate_file: str, aln_input_file: str, 
                            aln_file: str, cln_aln_file: str, bait_file: str, 
                            candidate_file: str, name: str, number: int, 
                            mode_aln: str, mode_tree: str, mafft: str, muscle: str, 
                            raxml: str, fasttree: str, iqtree: str, cpu_max: int, parallel_mode: str,
                            cand_color: str = "green", bait_color: str = "red") -> List[str]:

    """
    Handles parallel construction of alignments and phylogenetic trees.
    Returns a list of paths to the generated tree files.
    """
    # Input validation
    if not os.path.isfile(bait_file):
        raise FileNotFoundError(f"Bait sequence file not found: {bait_file}")
    if not os.path.isfile(candidate_file):
        raise FileNotFoundError(f"Candidate file not found: {candidate_file}")
    
    try:
        candidates = read_fasta(candidate_file)
        baits = read_fasta(bait_file)
        cand_ids = set(candidates.keys())
        bait_ids = set(baits.keys()) 
        if not candidates:
            raise ValueError("No candidates found in input file")
        
        # Calculate parallelization parameters
        candidates_number = min(num_process_candidates, len(candidates))
        num_par = math.ceil(len(candidates) / candidates_number)
        cpu_par = math.floor(cpu_max / num_par) if parallel_mode.upper() in ["YES", "Y"] else cpu_max
        
        # Adjust CPU allocation
        if parallel_mode.upper() in ["YES", "Y"] and cpu_par < 8:
            cpu_par = max(cpu_max, 1)
        
        if parallel_mode.upper() in ["NO", "N"]:
            num_par = 1
            candidates_number = len(candidates)
        
        # Create output directory
        Path(tree_output_folder).mkdir(parents=True, exist_ok=True)
        
        #print(f"Processing {len(candidates)} candidates in {num_par} parallel chunks")
        
        # Step 1: Create candidate files and combine with bait sequences
        candidate_keys = list(candidates.keys())
        # Shuffle candidate files
        random.shuffle(candidate_keys)

        for i in range(num_par):
            start = i * candidates_number
            end = min((i + 1) * candidates_number, len(candidates))
            
            cand_file = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_candidate_file}")
            aln_input = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_input_file}")

            if not os.path.isfile(aln_input) or os.path.getsize(aln_input) == 0:
                # Create candidate file
                with open(cand_file, "w") as out:
                    for c in candidate_keys[start:end]:
                        out.write(f'>{c}\n{candidates[c]}\n')
                
                # Combine bait and candidate files
                with open(aln_input, 'w') as outfile:
                    with open(bait_file, 'r') as infile1:
                        outfile.write(infile1.read())
                    with open(cand_file, 'r') as infile2:
                        outfile.write(infile2.read())
                #print(f"Combined bait and candidate file created: {aln_input}")
            else:
                print(f"Skipping creation of {aln_input}, file already exists and is not empty.")
        
        # Step 2: Create alignments using the Alignment class
        create_alignments(
            num_par=num_par,
            tree_output_folder=tree_output_folder,
            name=name,
            number=number,
            aln_input_file=aln_input_file,
            aln_file=aln_file,
            mode_aln=mode_aln,
            mafft=mafft,
            muscle=muscle,
            cpu_par=cpu_par,
            cpu_max=cpu_max
        )
        
        # Step 3: Create phylogenetic trees with improved process management
        tree_processes = []
        trees = []
        cpu_use = 0
        
        for i in range(num_par):
            aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_file}")
            cln_aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{cln_aln_file}")
            
            # Wait for CPU resources if needed
            while cpu_use + cpu_par > cpu_max:
                time.sleep(2)
                # Clean up finished processes and recalculate CPU usage
                active_processes = []
                for p in tree_processes:
                    if p.poll() is None:  # Process still running
                        active_processes.append(p)
                    else:
                        # Check for errors in completed processes
                        if p.returncode != 0:
                            print(f"Warning: Tree construction process failed with return code {p.returncode} for process {p.args}")
                            stderr_output = p.stderr.read().decode() if p.stderr else "No stderr output"
                            print(f"Stderr: {stderr_output}")
                tree_processes = active_processes
                cpu_use = len(tree_processes) * cpu_par
            
            # Perform alignment trimming using the Alignment class
            if not os.path.isfile(cln_aln):
                try:
                    alignment_trimming(aln, cln_aln, occupancy=0.01)
                except Exception as e:
                    print(f"Alignment trimming failed for chunk {i+1}: {str(e)}")
                    continue
            else:
                print(f"Skipping alignment trimming for {cln_aln}, file already exists.")
            
            # Tree construction
            tree_file = "" # Initialisiere tree_file hier
            try:
                if mode_tree == "raxml":
                    prefix = os.path.join(tree_output_folder, f"{name}{number}{i}_RAxML_tree")
                    tree_file = f"{prefix}.raxml.bestTree"
                    
                    if not os.path.isfile(tree_file):
                        cmd = [raxml, "--all", "--threads", str(cpu_par), "--model", "LG+G8+F", 
                                "--msa", cln_aln, "--prefix", prefix]
                        print(f"Executing RAxML command: {' '.join(cmd)}")
                        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                        tree_processes.append(p)
                        cpu_use += cpu_par
                        print(f"Started RAxML tree construction {i+1} with PID {p.pid}")
                    else:
                        print(f"Skipping RAxML for {tree_file}, file already exists.")

                elif mode_tree == "fasttree":
                    tree_file = os.path.join(tree_output_folder, f"{name}{number}{i}_FastTree_tree.tre")
                    
                    if not os.path.isfile(tree_file):
                        cmd = [fasttree, "-wag", "-nopr", cln_aln]
                        #print(f"Executing FastTree command: {' '.join(cmd)}")
                        with open(tree_file, 'w') as outfile:
                            p = subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE)
                        tree_processes.append(p)
                        cpu_use += cpu_par
                        #print(f"Started FastTree construction {i+1} with PID {p.pid}")
                    #else:
                        #print(f"Skipping FastTree for {tree_file}, file already exists.")
                
                elif mode_tree == "iqtree":
                    prefix = os.path.join(tree_output_folder, f"{name}{number}{i}_IQTREE_tree")
                    tree_file = f"{prefix}.treefile"  # IQ-TREE writes the final ML tree here
                    if not os.path.isfile(tree_file):
                        # IQ-TREE typical options:
                        # -s: alignment, -m: model, -T: threads, -pre: output prefix
                        #cmd = [iqtree, "-s", cln_aln, "-m", "LG+F+G", "-T", str(cpu_par), "-pre", prefix]
                        cmd = [iqtree, "-s", cln_aln, "-m", "LG+F+G", "-T", "AUTO", "-pre", prefix]
                        print(f"Executing IQ-TREE command: {' '.join(cmd)}")  # all comments in English
                        # Run asynchronously; IQ-TREE writes to files (.log, .iqtree, .treefile)
                        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                        tree_processes.append(p)
                        cpu_use += cpu_par
                        #print(f"Started IQ-TREE construction {i+1} with PID {p.pid}")

                else:
                    print(f"Unsupported tree mode: {mode_tree}. Skipping tree construction for chunk {i+1}.")
                    continue # Springe zur nächsten Iteration

                if tree_file: # Füge nur hinzu, wenn ein tree_file Name generiert wurde
                    trees.append(tree_file)
                        
            except Exception as e:
                print(f"Error starting tree construction {i+1}: {str(e)}")
                import traceback
                traceback.print_exc()

        # Wait for all tree construction processes to complete
        for p in tree_processes:
            try:
                stdout, stderr = p.communicate(timeout=3600)
                if p.returncode != 0:
                    print(f"Warning: Tree construction process failed with return code {p.returncode}")
                    print(f"Stdout: {stdout.decode()}")
                    print(f"Stderr: {stderr.decode()}")
            except subprocess.TimeoutExpired:
                print(f"Process {p.pid} timed out. Terminating...")
                p.kill()
                stdout, stderr = p.communicate()
                print(f"Stdout after kill: {stdout.decode()}")
                print(f"Stderr after kill: {stderr.decode()}")
            except Exception as e:
                print(f"Error while waiting for process {p.pid}: {str(e)}")
        
        #print(f"Successfully started {len(trees)} phylogenetic tree processes (some may have failed).")
        
        # Filter out trees that weren't actually created due to errors
        existing_trees = [tree for tree in trees if os.path.isfile(tree) and os.path.getsize(tree) > 0]
        #print(f"Found {len(existing_trees)} actually existing tree files.")

        return existing_trees
        
    except Exception as e:
        print(f"Error in parallel_tree_constructor: {str(e)}")
        import traceback
        traceback.print_exc()
        raise

def replace_spaces_in_tree(tree_file):
    """
    Replace spaces with underscores in leaf names of a Newick tree file.
    Returns path to the modified tree file.
    """
    with open(tree_file, 'r') as f:
        tree_str = f.read()
    
    # Pattern explanation:
    # ([^()\[\]:,]+) - captures leaf names (anything that's not tree structure characters)
    # (:\d*\.?\d*)? - optionally captures branch lengths that might follow
    modified_tree = re.sub(r'([^()\[\]:,]+)(:\d*\.?\d*)?', 
                          lambda m: m.group(1).replace(' ', '_') + (m.group(2) if m.group(2) else ''), 
                          tree_str)
    
    with open(tree_file, 'w') as f:
        f.write(modified_tree)
    
    return tree_file

# endregion

# region CandidateAssignment
"""
Region that provides functions for candidate assignment using phylogenetic trees.
Collection of functions to divide CYP candidates into ingroup and outgroup groups
as well as assigning ingroup CYPs to orthologs.
"""
def create_in_out_anno(baits_path: str, outgroup_path: str) -> Tuple[List[str], List[str]]:
    """Create in_list and out_list from bait and outgroup FASTA files.
    Returns tuple containing in_list and out_list.
    """
    # Read bait sequences
    try:
        baits = read_fasta(baits_path)
        in_list = list(baits.keys())
    except Exception as e:
        raise ValueError(f"Failed to process bait sequences: {str(e)}")
    
    # Initialize empty outgroup list
    out_list = []
    
    # Process outgroup if provided
    if outgroup_path:
        try:
            outgroup = read_fasta(outgroup_path)
            out_list = list(outgroup.keys())
        except Exception as e:
            print(f"Warning: Could not process outgroup sequences: {str(e)}")
    
    return in_list, out_list

def load_in_out_classification_file(tmp_result_table):
    """Load family member classification from file"""
    cyp_classification = {}
    with open(tmp_result_table, "r") as f:
        f.readline()  # remove header
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            cyp_classification.update({parts[1]: float(parts[4])})
            line = f.readline()
    return cyp_classification

def split_into_ingroup_and_outgroup(tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff, hmmsearch_results):
    """Split subject sequences into ingroup and outgroup based on reference baits"""
    
    # preparation of data structure
    groups_around_ref_gene = {}
    for gene in (in_list + out_list):
        groups_around_ref_gene.update({gene: []})
    
    # find node objects of reference genes
    tree = dendropy.Tree.get_from_path(tree_file, "newick")
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
    my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
    
    # Find available reference sequences in the tree
    ref_node_objects = {}
    available_nodes = set(node.label for node in tree.taxon_namespace)
    
    for node in tree.taxon_namespace:
        if node.label:
            node.label = node.label.replace(" ", "_")
            if node.label in groups_around_ref_gene:
                ref_node_objects.update({node.label: node})
    
    # Create lists of actually present reference genes
    present_in_list = [gene for gene in in_list if gene in ref_node_objects]
    present_out_list = [gene for gene in out_list if gene in ref_node_objects]
    
    # Log missing reference sequences
    missing_in_genes = set(in_list) - set(present_in_list)
    missing_out_genes = set(out_list) - set(present_out_list)
    
    if missing_in_genes:
        print(f"Warning: The following ingroup reference sequences are missing from tree {os.path.basename(tree_file)}: {', '.join(missing_in_genes)}")
    
    if missing_out_genes:
        print(f"Warning: The following outgroup reference sequences are missing from tree {os.path.basename(tree_file)}: {', '.join(missing_out_genes)}")
    
    # Check if we have any reference sequences at all
    if not present_in_list and not present_out_list:
        print(f"Error: No reference sequences found in tree {os.path.basename(tree_file)}. Cannot perform classification.")
        return {}
    
    if not present_in_list:
        print(f"Warning: No ingroup reference sequences found in tree {os.path.basename(tree_file)}.")
    
    if not present_out_list:
        print(f"Warning: No outgroup reference sequences found in tree {os.path.basename(tree_file)}.")
    
    # Build reference gene nodes list and dictionary
    ref_gene_nodes = []
    ref_gene_nodes_dict_to_check = {}
    
    for gene in (present_in_list + present_out_list):
        ref_gene_nodes.append(ref_node_objects[gene])
        ref_gene_nodes_dict_to_check.update({ref_node_objects[gene]: None})
    
    results = {}
    for i, t1 in enumerate(tree.taxon_namespace):
        try:
            ref_gene_nodes_dict_to_check[t1]
        except KeyError:  # only run analysis for non-reference sequences
            path_distances = []
            patristic_distances = {}
            
            # Calculate distance to all other sequences in tree
            for t2 in tree.taxon_namespace:
                try:
                    path_distance = pdm.path_edge_count(t1, t2)
                    patr_distance = pdm.patristic_distance(t1, t2)
                    path_distances.append({"key": t2.label, "val": path_distance})
                    patristic_distances.update({t2.label: patr_distance})
                except Exception as e:
                    # Handle potential distance calculation errors
                    print(f"Warning: Could not calculate distance between {t1.label} and {t2.label}: {str(e)}")
                    continue
            
            if not path_distances:
                # If no distances could be calculated, skip this sequence
                print(f"Warning: No distances calculated for sequence {t1.label}, skipping.")
                continue
            
            in_counter = 0
            out_counter = 0
            cand_counter = 0
            
            sorted_distances = sorted(path_distances, key=itemgetter("val"))
            neighbor_limit = min([len(path_distances), neighbour_cutoff])
            
            for each in sorted_distances[:neighbor_limit]:
                seq_label = each["key"]
                
                # Check if we have patristic distance for this sequence
                if seq_label not in patristic_distances:
                    continue
                    
                patr = patristic_distances[seq_label]
                
                # Exclude outliers based on mean distance
                if patr < mean_factor_cutoff * my_mean_nearest_taxon_distance:
                    if seq_label in present_in_list:
                        in_counter += 1
                    elif seq_label in present_out_list:
                        out_counter += 1
                    else:
                        cand_counter += 1
            
            # Handle HMM search results
            hmm = '-' if len(hmmsearch_results) == 0 else "yes" if t1.label in hmmsearch_results else "no"
            
            # Calculate score based on available reference sequences
            total_ref_neighbors = in_counter + out_counter
            
            if total_ref_neighbors > min_neighbour_cutoff:
                score = float(in_counter) / total_ref_neighbors if total_ref_neighbors > 0 else 0.0
                results.update({
                    t1.label: {
                        "score": score,
                        "in": in_counter,
                        "out": out_counter,
                        'hmm': hmm,
                        "tree": os.path.basename(tree_file)
                    }
                })
            else:
                results.update({
                    t1.label: {
                        "score": 0.0,
                        "in": in_counter,
                        "out": out_counter,
                        "hmm": hmm,
                        "tree": os.path.basename(tree_file)
                    }
                })
    
    #print(f"Successfully classified {len(results)} sequences from tree {os.path.basename(tree_file)}")
    return results

def candidates_baits_assignment(
    baits_path,
    tree_file,
    member_candidates,
    bait_groups,
    baits_info,
    static_pd="Y",
    factor=1.0,
    subfamily_threshold=1.1,
    family_threshold=2.7
):
    """
    Assigns candidates to baits based on patristic distances.
    Returns:
        subfamily_per_ref, family_per_ref
        (Thresholds are stored in each assignment entry)
    """
    ref_members = read_fasta(baits_path)
    subfamily_per_ref = {gene: [] for gene in ref_members.keys()}
    family_per_ref = {gene: [] for gene in ref_members.keys()}

    tree = dendropy.Tree.get_from_path(tree_file, "newick")
    tree.reroot_at_midpoint()
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)

    ref_node_objects = {}
    new_node_objects = {}
    for node in tree.taxon_namespace:
        if node.label is not None:
            node.label = node.label.replace(" ", "_")
        if node.label in ref_members:
            ref_node_objects[node.label] = node
        elif node.label in member_candidates:
            new_node_objects[node.label] = node

    ref_nodes = [ref_node_objects.get(g) for g in ref_members.keys() if g in ref_node_objects]
    candidate_nodes = [new_node_objects.get(g) for g in member_candidates]

    static_mode = str(static_pd).strip().upper() in ("Y", "YES")

    for candidate_node in candidate_nodes:
        distances = []
        for ref_node in ref_nodes:
            try:
                edge_dist = pdm.path_edge_count(candidate_node, ref_node)
                patr_dist = pdm.patristic_distance(candidate_node, ref_node)
                distances.append({
                    'ref_id': ref_node.label,
                    'group': bait_groups.get(ref_node.label, 'unknown'),
                    'edges': edge_dist,
                    'patr': patr_dist
                })
            except KeyError:
                continue

        # primary: patristic distance; tie-breaker: fewer edges preferred
        distances.sort(key=lambda x: (x['patr'], x['edges']))
        # distances.sort(key=lambda x: x['patr'])

        for dist in distances:
            patr = dist.get('patr')
            if patr is None:
                continue

            ref_id = dist['ref_id']
            # --- Threshold determination ---
            if static_mode:
                sub_thr = subfamily_threshold
                fam_thr = family_threshold
            else:
                sub_thr = subfamily_threshold
                fam_thr = family_threshold

                if ref_id in baits_info:
                    ref_info = baits_info[ref_id]

                    # --- Subfamily threshold ---
                    norm_flag = str(ref_info.get("Subfamily_Normal", "")).strip().lower()
                    if norm_flag == "y":
                        mean_val = ref_info.get("Subfamily_MeanDist")
                        std_val = ref_info.get("Subfamily_StdDist")
                        if mean_val is not None and mean_val != '' and std_val is not None and std_val != '':
                            sub_thr = (float(mean_val) + float(factor) * float(std_val))
                            if sub_thr > subfamily_threshold * 1.5:
                                sub_thr = subfamily_threshold * 1.5
                    elif norm_flag == "n":
                        median_val = ref_info.get("Subfamily_MedianDist")
                        iqr_val = ref_info.get("Subfamily_IQRDist")
                        if median_val is not None and median_val != '' and iqr_val is not None and iqr_val != '':
                            sub_thr = float(median_val) + float(factor) * float(iqr_val)
                            if sub_thr > subfamily_threshold * 1.5:
                                sub_thr = subfamily_threshold * 1.5

                    # --- Family threshold ---
                    fam_flag = str(ref_info.get("Family_Normal", "")).strip().lower()
                    if fam_flag == "y":
                        mean_val = ref_info.get("Family_MeanDist")
                        std_val = ref_info.get("Family_StdDist")
                        if mean_val is not None and mean_val != '' and std_val is not None and std_val != '':
                            fam_thr = float(mean_val) + float(factor) * float(std_val)
                            if fam_thr > family_threshold * 1.5:
                                fam_thr = family_threshold * 1.5
                            # elif fam_thr < sub_thr:
                            #     fam_thr = sub_thr
                            elif fam_thr < 2.0:
                                fam_thr = 2.0
                    elif fam_flag == "n":
                        median_val = ref_info.get("Family_MedianDist")
                        iqr_val = ref_info.get("Family_IQRDist")
                        if median_val is not None and median_val != '' and iqr_val is not None and iqr_val != '':
                            fam_thr = float(median_val) + float(factor) * float(iqr_val)
                            if fam_thr > family_threshold * 1.5:
                                fam_thr = family_threshold * 1.5
                            # elif fam_thr < sub_thr:
                            #     fam_thr = sub_thr
                            elif fam_thr < 2.0:
                                fam_thr = 2.0

            # --- Assignment ---
            assignment_data = {
                'candidate_id': candidate_node.label,
                'edges': dist['edges'],
                'patr': patr,
                'sub_thr': sub_thr,
                'fam_thr': fam_thr
            }

            if patr <= sub_thr:
                subfamily_per_ref.setdefault(ref_id, []).append(assignment_data)
            elif sub_thr < patr <= fam_thr:
                family_per_ref.setdefault(ref_id, []).append(assignment_data)

    return subfamily_per_ref, family_per_ref

def filter_baits_by_info(baits_info: Dict[str, Dict[str, str]], column: str, keyword: str) -> List[str]:
    """
    Filter bait IDs based on a given column and keyword in baits_info.
    Returns a list of bait sequence IDs.
    """
    return [
        seq_id for seq_id, info in baits_info.items()
        if column in info and keyword.lower() in info[column].lower()
    ]

def perform_candidate_tree_analysis(
    prefix: str,
    first_prefix: str,
    second_prefix: str,
    bait_fasta_path: str,
    clean_members: Dict[str, str],
    clean_members_file: str,
    baits_info: Dict[str, Dict[str, str]],
    subject_name_mapping_table: Dict[str, str],
    bait_groups: Dict[str, str],
    args,
    supplement_folder: str,
    tree_folder: str,
    group_around_ref_file: str,
    ref_mapping_file: str,
    filtered_fasta_file: str,
    linneage_specific_fasta: str
) -> None:
    """
    Performs phylogenetic analysis and candidate assignment based on 
    patristic distances to bait/reference sequences.

    This function handles tree construction (if not already done), assigns candidates
    to baits based on distance thresholds, writes groupings and assignments to 
    supplementary files, and extracts relevant sequences for final tree building.

    Returns:
    - filtered_tree_path: path to the first tree used for candidate-bait assignment
    - final_tree_path: path to the final tree with selected candidates and references
    """

    # Prepare file paths for intermediate alignment and tree steps
    aln_prefix = f"{args.name}{first_prefix}_"
    output_tree_folder = os.path.join(supplement_folder, f"{aln_prefix}first_tree/")
    aln_candidate_file = f"{aln_prefix}alignment_candidates.fasta"
    aln_input_file = f"{aln_prefix}alignment_input.fasta"
    aln_file = f"{aln_input_file}.aln"
    cln_aln_file = f"{aln_input_file}.aln.cln"

    # If the first tree does not exist, construct it
    if not (os.path.isdir(tree_folder) and any(f.startswith(first_prefix + "_first") for f in os.listdir(tree_folder))):
        tree_file = parallel_tree_constructor(
            len(clean_members), output_tree_folder, aln_candidate_file, aln_input_file,
            aln_file, cln_aln_file, bait_fasta_path, clean_members_file, f"{first_prefix}_first_", "",
            args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree, args.iqtree,
            args.cpu_max, "N"
        )
        filtered_tree_path = replace_spaces_in_tree(tree_file[0])
        shutil.copy(filtered_tree_path, tree_folder)
    else:
        # Use previously generated tree if available
        filtered_tree_path = os.path.join(tree_folder, f"{args.name}{first_prefix}_first_0_FastTree_tree.tre")

    # Assign candidates to references using patristic distance thresholds
    subfamily_mapping, family_mapping = candidates_baits_assignment(
        bait_fasta_path,
        filtered_tree_path,
        clean_members.keys(),
        bait_groups,
        baits_info,
        static_pd=args.static_pd,         
        factor=args.threshold_factor,                
        subfamily_threshold=args.subfamily_threshold,
        family_threshold=args.family_threshold
    )

    # Obatin group information
    groups = build_sequence_groups(
        baits_fasta_path=bait_fasta_path,
        clean_members=clean_members,
        subfamily_mapping=subfamily_mapping,
        family_mapping=family_mapping,
    ) 

    # Four group labels
    label_map_assignments = {
        "group1_baits": "Baits",
        "group2_subfamily": "Subfamily",
        "group3_family_only": "FamilyOnly",
        "group4_unassigned": "Unassigned",
    }
    palette_assignments = {
        "group1_baits": "#1f77b4",
        "group2_subfamily": "#2ca02c",
        "group3_family_only": "#ff7f0e",
        "group4_unassigned": "#7f7f7f",
    }

    write_itol_annotation(
        groups=groups,
        output_folder=tree_folder,
        dataset_label=f"{args.name}_{prefix}",
        filename_suffix=prefix,
        label_map=label_map_assignments,
        palette=palette_assignments
    )

    # Write groupings around each bait to file (only if not already created)
    if not os.path.isfile(group_around_ref_file) or os.path.getsize(group_around_ref_file) == 0:
        with open(group_around_ref_file, "w") as out:
            out.write(
                "Bait_Member\tBait_Evidence_Level\t"
                "Subfamily_Candidate_Members\tSubfamily_Thresholds\t"
                "Family_Candidate_Members\tFamily_Thresholds\n"
            )
            all_ref_genes = set(subfamily_mapping.keys()) | set(family_mapping.keys())
            for ref_gene in sorted(all_ref_genes):
                if ref_gene not in bait_groups:
                    continue

                # Subfamily members & thresholds
                subfamily_entries = subfamily_mapping.get(ref_gene, [])
                subfamily_names = [
                    subject_name_mapping_table[c['candidate_id']] for c in subfamily_entries
                ]
                subfamily_thresholds = [
                    str(c.get('sub_thr', '')) for c in subfamily_entries
                ]

                # Family members & thresholds
                family_entries = family_mapping.get(ref_gene, [])
                family_names = [
                    subject_name_mapping_table[c['candidate_id']] for c in family_entries
                ]
                family_thresholds = [
                    str(c.get('fam_thr', '')) for c in family_entries
                ]

                group = bait_groups.get(ref_gene, 'unknown')
                out.write(
                    f"{ref_gene}\t{group}\t"
                    f"{'; '.join(subfamily_names)}\t{'; '.join(subfamily_thresholds)}\t"
                    f"{'; '.join(family_names)}\t{'; '.join(family_thresholds)}\n"
                )

    # Write detailed mapping of candidates to subfamily/family assignments
    if not os.path.isfile(ref_mapping_file) or os.path.getsize(ref_mapping_file) == 0:
        with open(ref_mapping_file, "w") as out:
            header = (
                "Candidate_Member\tOriginal_ID\t"
                "Subfamily_Level_Ortholog\tSubfamily_Evidence_Level\tSubfamily_Family\tSubfamily_Subfamily\t"
                "Subfamily_Edge_Distance\tSubfamily_Patristic_Distance\tSubfamily_Threshold\t"
                "Family_Level_Ortholog\tFamily_Evidence_Level\tFamily_Family\tFamily_Subfamily\t"
                "Family_Edge_Distance\tFamily_Patristic_Distance\tFamily_Threshold\n"
            )
            out.write(header)

            all_candidates = set()
            candidate_to_subfamily = {}
            candidate_to_family = {}

            for ref, candidates in subfamily_mapping.items():
                for cand in candidates:
                    cid = cand['candidate_id']
                    all_candidates.add(cid)
                    candidate_to_subfamily.setdefault(cid, []).append({**cand, 'ref_id': ref})

            for ref, candidates in family_mapping.items():
                for cand in candidates:
                    cid = cand['candidate_id']
                    all_candidates.add(cid)
                    candidate_to_family.setdefault(cid, []).append({**cand, 'ref_id': ref})

            ref_name = ""
            for new_gene in sorted(all_candidates):
                subfamily_assign = sorted(candidate_to_subfamily.get(new_gene, []), key=lambda x: x['patr'])
                family_assign = sorted(candidate_to_family.get(new_gene, []), key=lambda x: x['patr'])

                max_len = max(len(subfamily_assign), len(family_assign))
                for i in range(max_len):
                    sub_ortho = subfamily_assign[i] if i < len(subfamily_assign) else {}
                    fam_ortho = family_assign[i] if i < len(family_assign) else {}

                    sub_id = sub_ortho.get('ref_id', '-')
                    sub_grp = bait_groups.get(sub_id, '-')
                    sub_fam = baits_info.get(sub_id, {}).get('Family', '-')
                    sub_subfam = baits_info.get(sub_id, {}).get('Subfamily', '-')
                    sub_edges = str(sub_ortho.get('edges', ''))
                    sub_patr = str(sub_ortho.get('patr', ''))
                    sub_thr = str(sub_ortho.get('sub_thr', ''))

                    fam_id = fam_ortho.get('ref_id', '-')
                    fam_grp = bait_groups.get(fam_id, '-')
                    fam_fam = baits_info.get(fam_id, {}).get('Family', '-')
                    fam_subfam = baits_info.get(fam_id, {}).get('Subfamily', '-')
                    fam_edges = str(fam_ortho.get('edges', ''))
                    fam_patr = str(fam_ortho.get('patr', ''))
                    fam_thr = str(fam_ortho.get('fam_thr', ''))

                    out.write("\t".join([
                        new_gene if (new_gene != ref_name and i == 0) else " " * len(new_gene),
                        subject_name_mapping_table.get(new_gene, '') if i == 0 else "",
                        sub_id, sub_grp, sub_fam, sub_subfam, sub_edges, sub_patr, sub_thr,
                        fam_id, fam_grp, fam_fam, fam_subfam, fam_edges, fam_patr, fam_thr
                    ]) + "\n")

                ref_name = new_gene

    # Extract candidate sequences used for final tree construction
    final_candidates = set()
    for mapping in [subfamily_mapping, family_mapping]:
        for ref_members in mapping.values():
            for entry in ref_members:
                final_candidates.add(entry['candidate_id'])

    
    #filtered_seqs = {seq_id: clean_members[seq_id] for seq_id in final_candidates if seq_id in clean_members}
    filtered_seqs = {seq_id: sequence for seq_id, sequence in clean_members.items() 
                 if seq_id in final_candidates}
    
    write_fasta(filtered_seqs, filtered_fasta_file)

    # Reversed filter (keeping seq_ids NOT IN final_candidates)
    linneage_specific_seqs = {seq_id: sequence for seq_id, sequence in clean_members.items() 
            if seq_id not in final_candidates}
    
    if linneage_specific_fasta != "":
        write_fasta(linneage_specific_seqs, linneage_specific_fasta)

    # Prepare file paths for final tree
    final_tree_output_folder = os.path.join(supplement_folder, f"{aln_prefix}final_tree/")
    final_aln_candidate_file = f"{aln_prefix}final_alignment_candidates.fasta"
    final_aln_input_file = f"{aln_prefix}final_alignment_input.fasta"
    final_aln_file = f"{final_aln_input_file}.aln"
    final_cln_aln_file = f"{final_aln_input_file}.aln.cln"

    # Construct final tree if not already done
    if not (os.path.isdir(tree_folder) and any(f.startswith(second_prefix + "_final") for f in os.listdir(tree_folder))):
        tree_file = parallel_tree_constructor(
            len(clean_members), final_tree_output_folder, final_aln_candidate_file, final_aln_input_file,
            final_aln_file, final_cln_aln_file, bait_fasta_path, filtered_fasta_file, f"{second_prefix}_final_", "",
            args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree, args.iqtree,
            args.cpu_max, "N"
        )
        final_tree_path = replace_spaces_in_tree(tree_file[0])
        shutil.copy(final_tree_path, tree_folder)
    else:
        final_tree_path = os.path.join(tree_folder, f"{args.name}{second_prefix}_final_0_FastTree_tree.tre")

    return filtered_tree_path, final_tree_path

def build_sequence_groups(
    baits_fasta_path: str,
    clean_members: Dict[str, str],
    subfamily_mapping: Dict[str, List[Dict]],
    family_mapping: Dict[str, List[Dict]],
) -> Dict[str, set[str]]:
    """
    Creates four groups:
      - group1_baits: all bait sequence IDs (FASTA headers)
      - group2_subfamily: all clean_members assigned to at least one subfamily
      - group3_family_only: all clean_members assigned to a family but not to a subfamily
      - group4_unassigned: all clean_members without any assignment
    Important: Groups are based on clean_members keys (candidates) and baits; each group is always created.
    """
    # All bait IDs
    bait_records = read_fasta(baits_fasta_path)
    group1_baits = set(bait_records.keys())

    # Collect candidate assignments
    sub_assigned = set()
    fam_assigned = set()

    for entries in subfamily_mapping.values():
        for entry in entries:
            cid = entry.get('candidate_id')
            if cid is not None:
                sub_assigned.add(cid)

    for entries in family_mapping.values():
        for entry in entries:
            cid = entry.get('candidate_id')
            if cid is not None:
                fam_assigned.add(cid)

    # Filter to existing clean_members
    clean_ids = set(clean_members.keys())

    # Group 2: Subfamily
    group2_subfamily = sub_assigned & clean_ids

    # Group 3: Family-only (in family, but not in subfamily)
    group3_family_only = (fam_assigned - sub_assigned) & clean_ids

    # Group 4: Unassigned (clean_members minus all assigned)
    assigned_any = (group2_subfamily | group3_family_only)
    group4_unassigned = clean_ids - assigned_any

    return {
        "group1_baits": group1_baits,
        "group2_subfamily": group2_subfamily,
        "group3_family_only": group3_family_only,
        "group4_unassigned": group4_unassigned,
    } 

def write_itol_annotation(
    groups: Dict[str, set[str]],
    output_folder: str,
    dataset_label: str,
    filename_suffix: str,
    label_map: Optional[Dict[str, str]] = None,
    palette: Optional[Dict[str, str]] = None,
    strip_width: int = 50,
    margin: int = 10,
) -> str:
    """
    Create an iTOL annotation.
    - groups: Dict[group_key -> set[seqIDs]]
    - label_map: optional mapping from group_key -> iTOL label (e.g., {"group1_baits": "Baits"})
    - palette: optional mapping from group_key -> hex color
    - strip_width, margin: display parameters
    Returns: path to the generated file.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Fallback labels from group keys if none are provided
    if label_map is None:
        label_map = {k: k for k in groups.keys()}

    # Fallback colors: deterministically generated if none are provided
    default_colors = [
        "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]
    if palette is None:
        palette = {}
    color_cycle = iter(default_colors)
    colors = {}
    for k in groups.keys():
        colors[k] = palette.get(k, next(color_cycle, "#000000"))

    out_path = os.path.join(output_folder, f"ITOL_Annotation_{filename_suffix}.txt")
    with open(out_path, "w") as fh:
        fh.write("DATASET_COLORSTRIP\n")
        fh.write("SEPARATOR\tTAB\n")
        fh.write(f"DATASET_LABEL\t{dataset_label}\n")
        fh.write("COLOR\t#000000\n")
        fh.write(f"STRIP_WIDTH\t{strip_width}\n")
        fh.write(f"MARGIN\t{margin}\n")
        fh.write("SHOW_INTERNAL\t0\n")
        fh.write("DATA\n")

        for group_key, ids in groups.items():
            color = colors[group_key]
            label = label_map.get(group_key, group_key)
            for seqid in sorted(ids):
                node_id = seqid
                fh.write(f"{node_id}\t{color}\t{label}\n")

    return out_path

#endregion

# region Domain Check
"""
Region to check sequences for domain hits using HMMER.
"""
def domain_check(
    sequences_file: str, 
    domains_file: str, 
    hmmresult: str, 
    hmmoutput: str, 
    hmmsearch: str = "hmmsearch", 
    hmm_score: float = 100.0
) -> Tuple[Dict, Dict, List]:
    """Screen sequences for domain hits using HMMER
    Returns three tuples (results_dict, best_domain_matches, domain_names).
    """
    # Run hmmsearch with domain table output
    run_hmmsearch(
        hmm_file=domains_file,
        subject_file=sequences_file,
        output_file=hmmresult,
        hmmsearch=hmmsearch,
        domtblout=True
    )

    # Read and process results
    hmmcheck = pd.read_csv(hmmresult, delim_whitespace=True, comment="#", header=None)
    hmmcheck = hmmcheck[hmmcheck[13] > hmm_score].sort_values(11).drop_duplicates([0, 3])

    domain_names = sorted(set(hmmcheck[3]))
    seqs = read_fasta(sequences_file)

    results = {}
    best_domain_match = {}

    for key in sorted(seqs.keys()):
        results[key] = {}
        hits = hmmcheck[hmmcheck[0] == key]
        
        best_evalue = float('inf')
        best_domain = "-"
        
        for dom in domain_names:
            match_rows = hits[hits[3] == dom]
            if not match_rows.empty:
                coords = match_rows.iloc[0]
                match_seq = seqs[key][coords[17]-1:coords[18]]
                results[key][str(dom)] = match_seq

                # Check for best domain match
                current_evalue = coords[11]
                if current_evalue < best_evalue:
                    best_evalue = current_evalue
                    best_domain = str(dom)
            else:
                results[key][str(dom)] = ""

        best_domain_match[key] = best_domain

    domain_names = [str(d) for d in domain_names]
    return results, best_domain_match, domain_names

def load_ortholog_family(ref_mapping_file):
    """
    Loads a mapping of candidate sequences to their assigned CYP family
    based on a reference mapping file.
    Returns a dictionary mapping candidate sequence IDs to family names.
    """
    candidate_to_family = {}
    
    # Return empty dictionary if file does not exist
    if not os.path.isfile(ref_mapping_file):
        return candidate_to_family

    with open(ref_mapping_file) as f:
        header = f.readline().strip().split("\t")

        # Determine relevant column indices
        idx_new_member = header.index("Candidate_Member")
        idx_ortho_family = header.index("Subfamily_Family") if "Subfamily_Family" in header else None
        idx_neigh_family = header.index("Family_Family") if "Family_Family" in header else None

        current_new_member = None

        for line in f:
            cols = line.strip().split("\t")
            if not any(cols):
                continue

            # Extract candidate ID from the first column (may span multiple rows)
            new_member = cols[idx_new_member].strip() if idx_new_member < len(cols) else ""
            if new_member:
                current_new_member = new_member
            if not current_new_member:
                continue

            # Extract family assignments from the respective columns
            ortho_family = cols[idx_ortho_family].strip() if idx_ortho_family is not None and idx_ortho_family < len(cols) else ""
            neigh_family = cols[idx_neigh_family].strip() if idx_neigh_family is not None and idx_neigh_family < len(cols) else ""

            # Prefer subfamily assignment unless it is missing or marked as invalid
            family = ortho_family if (ortho_family and ortho_family != "None" and ortho_family != "-") else neigh_family

            # Assign family if a valid value was found
            if family:
                candidate_to_family[current_new_member] = family

    return candidate_to_family

#endregion

# region Motif Check
"""
Region providing functions for static as well as HMM-based motif checks.
"""
def static_motif_check(fasta_file: str, motifs_path: str) -> Dict[str, Dict[str, str]]:
    """
    Checks protein sequences for the presence of specified motifs,
    considering motif order and positional boundaries.
    Returns a tuple containing a dictionary mapping sequence IDs to dictionaries of motif hits
    and a list of motif names in the order they appear in the input file.
    """

    # Load sequences and motif definitions
    sequences = read_fasta(fasta_file)

    motifs = []
    motif_names_list = []
    with open(motifs_path, 'r') as mfile:
        for line in mfile.readlines()[1:]:  # Skip header
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    name, pattern, _, min_p, max_p = parts[:5]
                    motifs.append({
                        'name': name,
                        'pattern': pattern,
                        'min_pos': int(min_p),
                        'max_pos': int(max_p)
                    })
                    motif_names_list.append(name)

    # Evaluate each sequence
    results = {}
    for seq_id, sequence in sequences.items():
        # Initialize result structure for this sequence
        seq_result = {motif['name']: "" for motif in motifs}

        # Precompute all allowed matches per motif
        all_matches = {}
        for motif in motifs:
            name = motif['name']
            pattern = motif['pattern']
            min_pos = motif['min_pos']
            max_pos = motif['max_pos']

            compiled_pattern = compile_motif_pattern(pattern)
            all_matches[name] = [
                (m.start(), m.end(), m.group())
                for m in compiled_pattern.finditer(sequence)
                if min_pos <= m.start() + 1 <= max_pos
            ]

        # Backtracking to find valid motif combinations in correct order
        valid_paths = []

        def backtrack(current_path: List[Tuple[str]], last_end_pos: int, remaining_motifs: List[Dict]):
            if not remaining_motifs:
                valid_paths.append(current_path)
                return

            current_motif = remaining_motifs[0]
            motif_name = current_motif['name']

            # Option 1: skip this motif
            backtrack(current_path + [(motif_name, "")], last_end_pos, remaining_motifs[1:])

            # Option 2: try using each compatible match
            for start, end, match in all_matches[motif_name]:
                if start >= last_end_pos:
                    backtrack(current_path + [(motif_name, match)], end, remaining_motifs[1:])

        backtrack([], 0, motifs)

        # Choose the best path (maximum number of matched motifs)
        if valid_paths:
            best_path = max(valid_paths, key=lambda p: sum(1 for _, m in p if m))
            for motif_name, match in best_path:
                if match:
                    seq_result[motif_name] = match

        results[seq_id] = seq_result

    return (results, motif_names_list)

def compile_motif_pattern(pattern: str) -> re.Pattern:
    """Compile a motif pattern into a regular expression."""
    # Standardize pattern format
    pattern = pattern.replace('(', '[').replace(')', ']').replace('/', '')
    
    # Convert to regex
    regex_pattern = []
    for char in pattern:
        if char == '[':
            regex_pattern.append('[')
        elif char == ']':
            regex_pattern.append(']')
        elif char == 'X':
            regex_pattern.append('.')
        else:
            regex_pattern.append(re.escape(char))
    
    return re.compile(''.join(regex_pattern))

def motif_check(
    sequences_file: str, 
    motifs_file: str, 
    hmmresult: str, 
    hmmoutput: str, 
    hmmsearch: str = "hmmsearch", 
    cEvalue: float = 0.001
) -> Tuple[Dict, List]:
    """Screen sequences for motifs using HMMER
        
    Returns a tuple of a dictionary mapping sequence IDs to dictionaries of motif hits
    and a list of motif names in the order they appear in the input file.
    """
    # Run hmmsearch with domain table output
    run_hmmsearch(
        hmm_file=motifs_file,
        subject_file=sequences_file,
        output_file=hmmresult,
        hmmsearch=hmmsearch,
        domtblout=True
    )

    # Process results
    hmmcheck = pd.read_csv(hmmresult, delim_whitespace=True, comment="#", header=None)
    hmmcheck = hmmcheck[hmmcheck[11] < cEvalue].sort_values(11).drop_duplicates([0, 3])
    
    motifs = list(set(hmmcheck[3]))
    seqs = read_fasta(sequences_file)    
    
    results = {}
    for key in sorted(seqs.keys()):       
        results[key] = {}
        hits = hmmcheck[hmmcheck[0] == key]
        for motif_id in motifs:
            try:
                coords = hits[hits[3] == motif_id].values[0]
                match = seqs[key][coords[17]-1:coords[18]]
                results[key][str(motif_id)] = match          
            except IndexError:
                results[key][str(motif_id)] = ""
    
    motifs = [str(m) for m in motifs]    
    return results, motifs

def create_motif_from_baits(
    protein_motifs_path: str,
    CYP_source: str,
    output_dir: str,
    mafft_available: bool,
    muscle_available: bool,
    flanklength: int = 5,
    mafft_path: str = "mafft",
    muscle_path: str = "muscle",
    hmmbuild_path: str = "hmmbuild"
) -> Tuple[str, Dict[str, int]]:
    """
    Creates a combined HMM profile from motif definitions and bait sequences.

    Steps:
    1. Read motif definitions (name, pattern, score threshold, positional constraints).
    2. Read bait sequences from FASTA file.
    3. Search for motifs using BLOSUM62-based scoring.
    4. Extract matched motif regions with flank extensions.
    5. Perform multiple sequence alignment (MAFFT or MUSCLE).
    6. Build HMM profile using `hmmbuild`.
    7. Merge all HMM profiles into one file.
    8. Clean up temporary files.

    Returns a tuple conataining a path to the final combined HMM file
    and a dictionary mapping motif names to their core lengths.
    """

    def _generate_combinations(motif: str) -> List[str]:
        parts, i = [], 0
        while i < len(motif):
            if motif[i] == '(':
                j = motif.find(')', i)
                parts.append(motif[i+1:j].split('/'))
                i = j + 1
            else:
                parts.append([motif[i]])
                i += 1
        combinations = ['']
        for part in parts:
            combinations = [c + option for c in combinations for option in part]
        return combinations

    def _score_window(window: str, motif: str) -> float:
        score = 0
        for a, b in zip(window, motif):
            if b == 'X':
                continue
            score += BLOSUM62.get((a, b), BLOSUM62.get((b, a), -4))
        return score

    temp_dir = os.path.join(output_dir, "temp_motif_build")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    try:
        # Step 1: Read motif definitions
        motifs = []
        with open(protein_motifs_path, 'r') as mfile:
            for line in mfile.readlines()[1:]:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        name, pattern, thr, min_p, max_p = parts[:5]
                        mlen = len(pattern.replace('(', '').replace(')', '').replace('/', ''))
                        motifs.append((name, pattern, float(thr), int(min_p), int(max_p), mlen))

        bait_sequences = read_fasta(CYP_source)
        motif_records = {name: [] for name, *_ in motifs}

        # Step 2: Search motifs in sequences
        for seq_id, seq in bait_sequences.items():
            seq = seq.upper()
            current_pos = 0
            for name, pattern, threshold, min_pos, max_pos, mlen in motifs:
                best_score, best_seq, best_start = -float('inf'), '', -1
                motif_variants = _generate_combinations(pattern)

                search_start = max(current_pos, min_pos - 1)
                search_end = min(len(seq) - mlen, max_pos - mlen + 1)

                flank_left = 2
                flank_right = 2 if mlen % 2 == 0 else 3
                target_len = mlen + flank_left + flank_right

                for i in range(search_start, search_end):
                    window = seq[i:i+mlen]
                    for variant in motif_variants:
                        sc = _score_window(window, variant)
                        if sc >= threshold and sc > best_score:
                            best_score, best_seq, best_start = sc, window, i

                if best_start != -1:
                    start_flank = max(0, best_start - flank_left)
                    end_flank = min(len(seq), best_start + mlen + flank_right)
                    frag = seq[start_flank:end_flank]
                    header = f"{seq_id}_{name}_{best_start+1}-{best_start+mlen}"
                    motif_records[name].append((header, frag))
                    current_pos = best_start + mlen

        # Step 3–6: Align and build HMMs
        individual_hmm_files = []
        for name, records in motif_records.items():
            if len(records) < 2:
                continue

            fasta_path = os.path.join(temp_dir, f"{name}.fasta")
            aln_path = os.path.join(temp_dir, f"{name}.aln")
            hmm_path = os.path.join(temp_dir, f"{name}.hmm")

            with open(fasta_path, 'w') as f:
                for header, seq in records:
                    f.write(f">{header}\n{seq}\n")

            if mafft_available:
                cmd = [mafft_path, "--auto", "--quiet", fasta_path]
                with open(aln_path, 'w') as out_file:
                    subprocess.run(cmd, stdout=out_file, check=True)
            elif muscle_available:
                cmd = [muscle_path, "-super5", fasta_path, "-output", aln_path]
                subprocess.run(cmd, check=True, capture_output=True)
            else:
                raise RuntimeError("No alignment tool (MAFFT or MUSCLE) is available.")

            cmd = [hmmbuild_path, hmm_path, aln_path]
            subprocess.run(cmd, check=True, capture_output=True)

            individual_hmm_files.append(hmm_path)

        # Step 7: Combine all HMMs into one file
        if not individual_hmm_files:
            raise RuntimeError("No HMMs could be built. Check motif definitions and bait sequences.")

        combined_hmm_path = os.path.join(output_dir, "motifs_from_baits.hmm")
        with open(combined_hmm_path, 'wb') as outfile:
            for hmm_file in individual_hmm_files:
                with open(hmm_file, 'rb') as infile:
                    outfile.write(infile.read())

        return combined_hmm_path, {name: mlen for name, *_ , mlen in motifs}

    except subprocess.CalledProcessError as e:
        print(f"Error during subprocess execution: {e.cmd}")
        print(f"Return Code: {e.returncode}")
        print(f"Stderr: {e.stderr.decode() if e.stderr else 'N/A'}")
        raise RuntimeError(f"A required tool failed to execute: {e.cmd}")
    
    finally:
        # Step 8: Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

def parse_motif_names_lengths(hmm_file: str) -> tuple[list[str], dict[str, int]]:
    """
    Reads HMM file and extracts names and lengths of motifs.
    Returns tuple of names and lengths.
    """
    names = []
    lengths = {}
    current_motif = None
    with open(hmm_file) as f:
        for line in f:
            if line.startswith("NAME"):
                current_motif = line.strip().split()[1]
                names.append(current_motif)
            elif line.startswith("LENG") and current_motif:
                lengths[current_motif] = int(line.strip().split()[1])
    return names, lengths

def write_motif_output(
    summary_file,
    results_file,
    motif_data,
    motif_names,
    raw_lengths,
    flanks,
    final_members,
    subject_name_mapping_table
):
    """
    Writes output files for motifs.
    """
    with open(summary_file, "w") as out1, open(results_file, "w") as out2:
        # Write headers
        out1.write("RefMember\t" + "\t".join(motif_names) + "\n")
        out2.write("RefMember\t" + "\t".join(motif_names) + "\n")

        for candidate in sorted(final_members.keys()):
            # First column: mapped name from subject_name_mapping_table
            new_line_details = [subject_name_mapping_table[candidate]]
            new_line_summary = [subject_name_mapping_table[candidate]]

            for mot in motif_names:
                seq = motif_data[candidate].get(mot, "")
                motif_len = raw_lengths.get(mot, 8)  # Default motif length is 8
                flank_left, flank_right = flanks.get(mot, (0, 0))

                if seq:
                    # If a match is found, remove the flanking residues
                    if len(seq) >= flank_left + flank_right:
                        trimmed_seq = seq[flank_left:len(seq) - flank_right]
                    else:
                        trimmed_seq = seq  # Fallback: use sequence as-is
                else:
                    trimmed_seq = "#" * motif_len  # No hit: use placeholder

                new_line_details.append(trimmed_seq)
                new_line_summary.append("YES" if seq else "NO")

            # Write to files
            out1.write("\t".join(new_line_summary) + "\n")
            out2.write("\t".join(new_line_details) + "\n")

# endregion

# region ExpressionMatrix
"""
Region to handle expression data, filtering candidates für expression values and metadata if provided.
"""
def Candidate_Expression(expression_matrix, candidates, min_avg_tpm, min_single_tpm):
    """
    Reads the expression matrix and returns a list of expressed CYP candidates.
    """
    expressed = {}
    potential_pseudogenes = {}

    with open(expression_matrix, 'r') as f:
        lines = f.readlines()

    if not lines:
        print("Warning: Expression matrix is empty")
        return [], expressed, potential_pseudogenes

    # Get header
    header = lines[0].strip()
    filtered_lines = [header]
    
    candidate_ids = set(candidates.keys())
    found_candidates = set()

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        gene_id = parts[0]

        if gene_id in candidate_ids:
            try:
                # Convert expression values to floats (excluding gene ID)
                tpm_values = list(map(float, parts[1:]))
                avg_tpm = sum(tpm_values) / len(tpm_values)
                max_tpm = min(tpm_values)

                if avg_tpm >= min_avg_tpm or max_tpm >= min_single_tpm:
                    expressed[gene_id] = candidates[gene_id]
                    filtered_lines.append(line)
                else:
                    potential_pseudogenes[gene_id] = candidates[gene_id]

                found_candidates.add(gene_id)
            except ValueError:
                print(f"Warning: Non-numeric TPM values for {gene_id}, skipping this entry.")

    missing_candidates = candidate_ids - found_candidates
    if missing_candidates:
        print(f"Warning: {len(missing_candidates)} candidates not found in expression matrix")

    return filtered_lines, expressed, potential_pseudogenes

def count_commas(line):
    """
    Count the number of commas, ignoring commas within double quotes.
    """
    in_quotes = False
    comma_count = 0
    for char in line:
        if char == '"':
            in_quotes = not in_quotes
        elif char == ',' and not in_quotes:
            comma_count += 1
    return comma_count

def use_metadata(filtered_expression_file, metadata_path, metadata_expression_file, metadata_only_expression_file):
    """"
    Use metadata to filter and annotate an expression matrix.
    """
    metadata_dict = {}
    metadata_ids = set()
    header_positions = {}

    unused_metadata_path = os.path.join(os.path.dirname(filtered_expression_file), "unused_metadata.txt")
    unused_lines = []

    # Step 1: Read metadata file
    with open(metadata_path, 'r', newline='') as meta_file:
        all_lines = meta_file.readlines()

    if len(all_lines) < 2:
        print("Metadata file is empty or has no data rows.")
        return

    meta_header = next(csv.reader([all_lines[0]]))
    
    # Look for 'tissue' and all columns starting with 'condition'
    for idx, col_name in enumerate(meta_header):
        if col_name == "tissue" or col_name.startswith("condition"):
            header_positions[col_name] = idx

    # Determine the expected number of commas (based on first data row)
    expected_commas = count_commas(all_lines[1])

    # Process metadata rows
    for line in all_lines[1:]:
        if count_commas(line) != expected_commas:
            unused_lines.append(line.strip())
            continue
        row = next(csv.reader([line]))
        sample_id = row[0]
        metadata_ids.add(sample_id)
        metadata_dict[sample_id] = [row[idx] for col, idx in header_positions.items()]

    # Write inconsistent metadata lines (wrong comma count) to file
    if unused_lines:
        with open(unused_metadata_path, 'w') as out_unused:
            out_unused.write("These metadata lines had inconsistent comma counts (outside of quotes) and were ignored:\n")
            for l in unused_lines:
                out_unused.write(l + '\n')
        print(f"Skipped {len(unused_lines)} inconsistent metadata lines. See: {unused_metadata_path}")

    # Step 2: Read the filtered expression matrix
    with open(filtered_expression_file, 'r') as f:
        lines = f.readlines()

    if not lines:
        print("Warning: Filtered expression matrix is empty")
        return

    header_parts = lines[0].strip().split('\t')
    sample_ids = header_parts[1:]
    new_sample_ids = []

    matched_indices = []
    for idx, sid in enumerate(sample_ids):
        if sid in metadata_dict:
            new_id = f"{sid}(" + ", ".join(metadata_dict[sid]) + ")"
            new_sample_ids.append(new_id)
            matched_indices.append(idx + 1)  # +1 because index 0 is gene ID
        else:
            new_sample_ids.append(f"{sid}(no metadata)")

    # Step 3: Write new output files ---
    # metadata_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix")
    # metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")

    with open(metadata_expression_file, 'w') as out_all, open(metadata_only_expression_file, 'w') as out_meta:
        # Write headers
        out_all.write(header_parts[0] + '\t' + '\t'.join(new_sample_ids) + '\n')
        out_meta.write(header_parts[0] + '\t' + '\t'.join([new_sample_ids[i - 1] for i in matched_indices]) + '\n')

        # Write expression values
        for line in lines[1:]:
            parts = line.strip().split('\t')
            gene_id = parts[0]
            expr_values = parts[1:]
            out_all.write(gene_id + '\t' + '\t'.join(expr_values) + '\n')
            meta_values = [expr_values[i - 1] for i in matched_indices]
            out_meta.write(gene_id + '\t' + '\t'.join(meta_values) + '\n')

    #print(f"Metadata-annotated expression matrix written to: {metadata_expression_file}")
    #print(f"Metadata-only expression matrix written to: {metadata_only_expression_file}")

def write_highest_expression_info(filtered_expression_file, highest_expression_file, metadata_expression_file=None, metadata_only_expression_file=None):
    """
    Write information about the highest expressed gene for each sample.
    """
    highest_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_highest_expression")

    # Step 1: Read original expression values
    with open(filtered_expression_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().split('\t')
    sample_ids = header[1:]
    gene_rows = [line.strip().split('\t') for line in lines[1:]]

    # Step 2: Extract metadata from annotated expression matrix (if available)
    metadata_info = {}
    metadata_available = metadata_expression_file and os.path.isfile(metadata_expression_file)

    if metadata_available:
        with open(metadata_expression_file, 'r') as f:
            meta_header = f.readline().strip().split('\t')[1:]
            for idx, sid in enumerate(sample_ids):
                if idx < len(meta_header):
                    match = re.search(r"\((.*?)\)", meta_header[idx])
                    metadata_info[sid] = match.group(1) if match else "no metadata"
                else:
                    metadata_info[sid] = "no metadata"
    else:
        metadata_info = {sid: "no metadata" for sid in sample_ids}

    # Step 3: If metadata_only_expression_file is provided, parse it for AltSample + TPM info
    meta_expr_data = {}
    if metadata_only_expression_file and os.path.isfile(metadata_only_expression_file):
        with open(metadata_only_expression_file, 'r') as f:
            meta_lines = f.readlines()
        meta_header_ids = meta_lines[0].strip().split('\t')[1:]
        for line in meta_lines[1:]:
            parts = line.strip().split('\t')
            gene_id = parts[0]
            try:
                values = list(map(float, parts[1:]))
            except ValueError:
                continue
            if values:
                max_val = max(values)
                max_idx = values.index(max_val)
                alt_sample = meta_header_ids[max_idx]
                match = re.search(r"(.+?)\((.*?)\)", alt_sample)
                if match:
                    alt_id = match.group(1)
                    alt_meta = match.group(2)
                    meta_expr_data[gene_id] = (f"{alt_id}({alt_meta})", max_val)
                else:
                    meta_expr_data[gene_id] = (alt_sample, max_val)
    else:
        meta_expr_data = {}

    # Step 4: Write output file
    with open(highest_expression_file, 'w') as out:
        out.write("GeneID\tBestSample\tBestTPM\tAltSample\tAltTPM\n")
        for row in gene_rows:
            gene_id = row[0]
            try:
                expr_vals = list(map(float, row[1:]))
            except ValueError:
                continue

            max_tpm = max(expr_vals)
            max_idx = expr_vals.index(max_tpm)
            sample_id = sample_ids[max_idx]
            condition = metadata_info.get(sample_id, "no metadata")
            sample_with_meta = f"{sample_id}({condition})"

            alt_sample = "-"
            alt_tpm = "-"

            if condition == "no metadata" and gene_id in meta_expr_data:
                alt_sample, alt_tpm_val = meta_expr_data[gene_id]
                alt_tpm = f"{alt_tpm_val:.3f}"

            out.write(f"{gene_id}\t{sample_with_meta}\t{max_tpm:.3f}\t{alt_sample}\t{alt_tpm}\n")

    #print(f"Highest expression summary written to: {output_file}")

    return highest_expression_file

#endregion

# region Paralogs
"""Region containing functions to search for paralog groups in phylogenetic trees, taking into account expression data if provided."""
def establish_paralog_groups( tree_file, member_candidates, dist_cutoff_factor ):
    """
    Establish paralog groups based on phylogenetic trees.
    """

    candidate_mapping_table = {}
    for gene in member_candidates:    #candidate genes of new species
        candidate_mapping_table.update( { gene: None } )
    
    # --- find node objects of reference genes --- #
    tree = dendropy.Tree.get_from_path( tree_file, "newick" )
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
    my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
    
    new_node_objects = {}    #get new family candidate node objects
    for node in tree.taxon_namespace:
        if node.label is not None:
            node.label = node.label.replace(" ", "_")
        try:
            candidate_mapping_table[ node.label ]
            new_node_objects.update( { node.label: node } )
        except KeyError:
            pass
    
    candidate_gene_nodes = [] # nodes of new candidates
    for gene in member_candidates:
        candidate_gene_nodes.append( new_node_objects[ gene ] )
    
    black_list = {}
    paralog_collection = []
    for i, t1 in enumerate( candidate_gene_nodes ):
        try:
            black_list[ t1.label ]
        except KeyError:
            paralogs = [ t1.label ]
            edge_distances = []
            patr_distances = {}
            for t2 in tree.taxon_namespace:    #calculate distance to all other sequences in tree
                try:
                    black_list[ t2.label ]
                except KeyError:
                    if t1.label != t2.label:
                        edge_distances.append( { 'id': t2.label, 'dist': pdm.path_edge_count( t1, t2) } )
                        patr_distances.update( { t2.label: pdm.patristic_distance( t1, t2 ) } )
            for each in list( sorted( edge_distances, key=itemgetter('dist') ) ):
                try:
                    candidate_mapping_table[ each['id'] ]
                    if patr_distances[ each['id'] ] < ( my_mean_nearest_taxon_distance*dist_cutoff_factor):
                        paralogs.append( each['id'] )
                        black_list.update( { each['id']: None } )
                except KeyError:
                    break    #next neighbour is not a new candidate => break extension of paralog group
            paralog_collection.append( paralogs )
            black_list.update( { t1.label: None } )

    return paralog_collection

def establish_paralog_groups_with_expression_filter(
    tree_file, member_candidates, dist_cutoff_factor,
    expression_matrix_path, min_paralog_tpm):
    """Create paralog groups while considering tissue-/condition-specific expression."""

    # Step 1: Read expression matrix and metadata
    expression_data = {}
    sample_conditions = {}

    with open(expression_matrix_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        header = lines[0].split('\t')
        sample_ids = header[1:]

        for idx, sid in enumerate(sample_ids):
            match = re.search(r"\((.*?)\)", sid)
            sample_conditions[sid] = match.group(1) if match else "NA"

        for line in lines[1:]:
            parts = line.split('\t')
            gene_id = parts[0]
            values = list(map(float, parts[1:]))
            expression_data[gene_id] = dict(zip(sample_ids, values))

    # Step 2: Load tree and identify candidate nodes
    candidate_mapping_table = {gene: None for gene in member_candidates}
    tree = dendropy.Tree.get_from_path(tree_file, "newick")
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
    mean_dist = pdm.mean_nearest_taxon_distance()

    new_node_objects = {node.label: node for node in tree.taxon_namespace if node.label in candidate_mapping_table}
    candidate_gene_nodes = [new_node_objects[gene] for gene in member_candidates]

    # Step 3: Grouping with expression-based separation
    black_list = {}
    paralog_collection = []

    for i, t1 in enumerate(candidate_gene_nodes):
        if t1.label in black_list:
            continue

        paralogs = [t1.label]
        edge_distances = []
        patr_distances = {}

        for t2 in tree.taxon_namespace:
            if t2.label in black_list or t2.label == t1.label:
                continue
            edge_distances.append({'id': t2.label, 'dist': pdm.path_edge_count(t1, t2)})
            patr_distances[t2.label] = pdm.patristic_distance(t1, t2)

        for entry in sorted(edge_distances, key=itemgetter('dist')):
            candidate = entry['id']
            if candidate not in candidate_mapping_table:
                break
            if patr_distances[candidate] < mean_dist * dist_cutoff_factor:
                # --- Expression-based filter --- #
                keep_separate = False
                if t1.label in expression_data and candidate in expression_data:
                    conds_t1 = {}
                    conds_t2 = {}
                    for sid in sample_ids:
                        expr1 = expression_data[t1.label].get(sid, 0.0)
                        expr2 = expression_data[candidate].get(sid, 0.0)
                        cond = sample_conditions.get(sid, "NA")

                        if expr1 >= min_paralog_tpm:
                            conds_t1[cond] = expr1
                        if expr2 >= min_paralog_tpm:
                            conds_t2[cond] = expr2

                    # If there are conditions in which only one of the two genes is expressed
                    unique_conds = set(conds_t1.keys()).symmetric_difference(set(conds_t2.keys()))
                    if unique_conds:
                        keep_separate = True

                if not keep_separate:
                    paralogs.append(candidate)
                    black_list[candidate] = None

        paralog_collection.append(paralogs)
        black_list[t1.label] = None

    return paralog_collection

def get_represenative_paralog_per_group(paralog_groups, clean_members, repr_clean_file, metadata_expression_path=None):
    """Selects the longest sequence per group and optionally appends tissue-/condition-specific metadata."""
    expression_annotation = {}

    if metadata_expression_path:
        with open(metadata_expression_path, 'r') as f:
            header = f.readline().strip().split('\t')
            sample_ids = header[1:]
            for sid in sample_ids:
                match = re.search(r"\((.*?)\)", sid)
                expression_annotation[sid] = match.group(1) if match else "NA"

    paralog_representatives = {}
    annotation_dict = {}  # for metadata output

    with open(repr_clean_file, "w") as out:
        for group in paralog_groups:
            seqs = [
                {'id': each, 'len': len(clean_members[each]), 'seq': clean_members[each]}
                for each in group
            ]
            representative = sorted(seqs, key=itemgetter('len', 'id'))[-1]
            rep_id = representative['id']
            out.write(f">{rep_id}\n{representative['seq']}\n")
            paralog_representatives[rep_id] = representative['seq']

            if metadata_expression_path:
                annots = []
                for gene in group:
                    if gene == rep_id:
                        continue
                    annots.append(gene)
                annotation_dict[rep_id] = annots

    # Extension of 08_representative_paralogs.txt
    if metadata_expression_path:
        meta_output_file = repr_clean_file.replace(".fasta", ".txt")
        with open(meta_output_file, "w") as out:
            out.write("RepresentativeSeqID\tMembersOfParalogGroup\n")
            for rep, group in paralog_representatives.items():
                meta_info = f"({','.join(annotation_dict.get(rep, []))})" if rep in annotation_dict else ""
                out.write(f"{rep}{meta_info}\t" + ";".join(group for group in paralog_groups if rep in group[0])[0] + "\n")

    return paralog_representatives

def build_paralog_groups_for_itol(
    baits_fasta_path: str,
    final_members: Dict[str, str],
    paralog_groups: List[List[str]],
    representatives: Dict[str, str],
) -> Dict[str, set[str]]:
    """
    Creates three groups for Step 8 (representative paralog groups):
      - group1_baits: all bait sequence IDs (from FASTA headers)
      - group2_representatives: the representative paralogs (keys from 'representatives')
      - group3_collapsed_paralogs: all paralog members from 'paralog_groups'
        that are not in 'representatives' (collapsed, not used in the tree)
    final_members: dict of sequences that are part of the final selection (IDs -> sequence)
    """
    # Group 1: Baits
    bait_records = read_fasta(baits_fasta_path)
    group1_baits = set(bait_records.keys())

    # Group 2: Representative paralogs (IDs only)
    rep_ids = set(representatives.keys())

    # Group 3: Collapsed paralogs = all members from groups except representatives
    collapsed = set()
    for grp in paralog_groups:
        for gid in grp:
            if gid not in rep_ids:
                collapsed.add(gid)

    # Optional: restrict to known IDs (if desired)
    # Only consider candidates present in final_members
    cand_ids = set(final_members.keys())
    group2_representatives = rep_ids & cand_ids if cand_ids else rep_ids
    group3_collapsed_paralogs = collapsed & cand_ids if cand_ids else collapsed

    return {
        "group1_baits": group1_baits,
        "group2_representatives": group2_representatives,
        "group3_collapsed_paralogs": group3_collapsed_paralogs,
    }

#endregion

# region Summary
"""Region to collect and summarize data in a summary table."""
def collect_summary_data(
        ref_mapping_file,
        baits_info,
        best_domain_match_file=None,
        motif_check_file_summary=None,
        highest_expression_file=None,
        paralog_group_file=None
        ):
    """"
    Collects data for the summary table.
    Returns a dictionary with summary data.
    """
    
    summary_data = {}

    # Step 1: Read ref_mapping_file → Column 2 (seq_id) and 3 (gene_id)
    with open(ref_mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            seq_id = parts[1]
            gene_id = parts[2]
            summary_data[seq_id] = [seq_id, gene_id]

    # Add Family/Subfamily info from baits_info
    for seq_id in summary_data:
        gene_id = summary_data[seq_id][1]  
        baits_info_row = baits_info.get(gene_id, {})
        col4 = baits_info_row.get("Family", "-")
        col5 = baits_info_row.get("Subfamily", "-")
        summary_data[seq_id].extend([col4, col5])

    # Step 2: Add best domain matches → Column 2 by ID in Column 1
    best_domain_dict = {}
    if best_domain_match_file:
        with open(best_domain_match_file, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    best_domain_dict[parts[0]] = parts[1]
    for seq_id in summary_data:
        summary_data[seq_id].append(best_domain_dict.get(seq_id, "-"))

    # Step 3: Add motif presence → Columns 2–7 by ID in Column 1
    motif_dict = {}
    if motif_check_file_summary:
        with open(motif_check_file_summary, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    motif_dict[parts[0]] = parts[1:7]
    for seq_id in summary_data:
        summary_data[seq_id].extend(motif_dict.get(seq_id, ["-"] * 6))

    # Step 4: Add highest expression info → Columns 2–5 by ID in Column 1
    expression_dict = {}
    if highest_expression_file:
        with open(highest_expression_file, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    expression_dict[parts[0]] = parts[1:5]
    for seq_id in summary_data:
        summary_data[seq_id].extend(expression_dict.get(seq_id, ["-"] * 4))

    # Step 5: Determine representative paralogs → ID before ":" = representative
    repr_ids = set()
    if paralog_group_file:
        with open(paralog_group_file, 'r') as f:
            next(f)
            for line in f:
                if "\t" in line:
                    parts = line.strip().split("\t")
                    repr_id = parts[0].strip()
                    repr_ids.add(repr_id)
    for seq_id in summary_data:
        summary_data[seq_id].append("YES" if seq_id in repr_ids else "NO")

    return summary_data

#endregion

#region Output
"""
Region providing functions to create various output files, including HTML, SVG and PDF/PNG
"""
def baits_info_to_txt(baits_info: Dict[str, Dict[str, str]], headers: List[str], output_path: str) -> None:
    """
    Write the baits_info dictionary to a tab-separated .txt file using the original column order.
    """
    if not baits_info or not headers:
        print(f"No data to write to {output_path}.")
        return

    try:
        with open(output_path, "w") as out:
            out.write("\t".join(headers) + "\n")
            for seq_id, subdict in sorted(baits_info.items()):
                row = [seq_id]  # ID is separate (first column)
                for h in headers[1:]:  # Skip 'ID', already used
                    row.append(subdict.get(h, ""))
                out.write("\t".join(row) + "\n")
    except Exception as e:
        raise IOError(f"Could not write baits info to {output_path}: {e}")

def group_by_leading_number(files):
    """
    Group files by leading number.
    """
    grouped = defaultdict(list)
    for html_file, orig_file in files:
        match = re.match(r"(\d+)_", orig_file)
        group = match.group(1) if match else "misc"
        grouped[group].append((html_file, orig_file))
    return grouped

def create_output_html_index(txt_folder, html_folder, link_overview_file,
                             tree_folder, fasta_folder, heatmap_folder=None):
    """
    Create an index HTML file for the output files.
    """
    os.makedirs(html_folder, exist_ok=True)

    # TXT-Files
    txt_files = [f for f in os.listdir(txt_folder) if f.endswith(".txt")]
    txt_grouped = defaultdict(list)

    for txt_file in txt_files:
        match = re.match(r"(\d+)_", txt_file)
        group = match.group(1) if match else "misc"
        txt_grouped[group].append(txt_file)

        txt_path = os.path.join(txt_folder, txt_file)
        html_file = txt_file.replace(".txt", ".html")
        html_path = os.path.join(html_folder, html_file)

        with open(txt_path, "r") as f_in, open(html_path, "w") as f_out:
            lines = [line.rstrip('\n') for line in f_in if line.strip() != ""]
            f_out.write("<html><body>\n")
            f_out.write(f"<h2>{txt_file}</h2>\n")
            rel_link = os.path.relpath(link_overview_file, html_folder)
            f_out.write(f"<p><a href='{os.path.basename(rel_link)}'>Overview</a></p>\n")
            f_out.write("<table border='1' cellspacing='0' cellpadding='3'>\n")
            for i, line in enumerate(lines):
                cols = line.split("\t")
                f_out.write("<tr>")
                for col in cols:
                    tag = "th" if i == 0 else "td"
                    safe_col = col.strip() if col.strip() != "" else "&nbsp;"
                    f_out.write(f"<{tag}>{safe_col}</{tag}>")
                f_out.write("</tr>\n")
            f_out.write("</table>\n</body></html>")

    # Heatmaps
    heatmap_htmls = []
    if os.path.isdir(heatmap_folder):
        heatmap_htmls = []
        for f in os.listdir(heatmap_folder):
            if f.endswith((".svg")):
                src_path = os.path.join(heatmap_folder, f)
                html_name = os.path.splitext(f)[0] + ".html"
                target_path = os.path.join(html_folder, html_name)
                convert_file_to_html(src_path, target_path, f, link_overview_file)
                heatmap_htmls.append((html_name, f))

    # Trees
    tree_htmls = []
    for f in os.listdir(tree_folder):
        if f.endswith((".svg")):
            src_path = os.path.join(tree_folder, f)
            html_name = os.path.splitext(f)[0] + ".html"
            target_path = os.path.join(html_folder, html_name)
            convert_file_to_html(src_path, target_path, f, link_overview_file)
            tree_htmls.append((html_name, f))

    # FASTA
    fasta_htmls = []
    for f in os.listdir(fasta_folder):
        if f.endswith(".fasta"):
            src_path = os.path.join(fasta_folder, f)
            html_name = os.path.splitext(f)[0] + ".html"
            target_path = os.path.join(html_folder, html_name)
            convert_file_to_html(src_path, target_path, f, link_overview_file)
            fasta_htmls.append((html_name, f))

    # Group by leading number
    txt_htmls = [(f.replace(".txt", ".html"), f) for group in txt_grouped.values() for f in group]
    if heatmap_htmls == []:
        other_grouped = group_by_leading_number(tree_htmls + fasta_htmls)
    else:
        other_grouped = group_by_leading_number(heatmap_htmls + tree_htmls + fasta_htmls)

    # Merged grouping
    merged_grouped = defaultdict(list)
    for group, files in txt_grouped.items():
        for f in files:
            merged_grouped[group].append((f.replace(".txt", ".html"), f))
    for group, files in other_grouped.items():
        merged_grouped[group].extend(files)

    # Write HTML overview
    with open(link_overview_file, "w") as fout:
        fout.write("<html><body>\n<h1>CYP Output Files</h1>\n")

        for group in sorted(merged_grouped.keys(), key=lambda x: int(x) if x.isdigit() else 999):
            fout.write(f"<h2>Section {group}</h2><ul>\n")
            for html_name, orig_file in sorted(merged_grouped[group]):
                fout.write(f"<li><a href='{html_name}'>{orig_file}</a></li>\n")
            fout.write("</ul><br>\n")

        fout.write("</body></html>\n")

def convert_file_to_html(src_path, target_path, filename, link_overview_file):
    """
    Creates a HTML file with embedded content from a file.
    Supports .txt, .png, .pdf, .fa/.fasta (as text display).
    Inserts a link to the overview page at the beginning.
    """
    ext = os.path.splitext(filename)[1].lower()
    os.makedirs(os.path.dirname(target_path), exist_ok=True)

    with open(target_path, "w") as f_out:
        f_out.write("<html><body>\n")
        f_out.write(f"<h2>{filename}</h2>\n")

        # Link to HTML overview
        rel_link = os.path.relpath(link_overview_file, os.path.dirname(target_path))
        f_out.write(f"<p><a href='{os.path.basename(rel_link)}'>Overview</a></p>\n")

        # Content deprending on file type
        if ext in [".svg"]:
            f_out.write(f"<object data='../{os.path.basename(os.path.dirname(src_path))}/{filename}' type='image/svg+xml' width='100%' height='800px'>\n")

        elif ext in [".png"]:
            f_out.write(f"<img src='../{os.path.basename(os.path.dirname(src_path))}/{filename}' style='max-width:100%;'/>\n")

        elif ext in [".pdf"]:
            f_out.write(f"<object data='../{os.path.basename(os.path.dirname(src_path))}/{filename}' type='application/pdf' width='100%' height='800px'>\n")
            f_out.write("PDF not shown. <a href='../{0}/{1}'>Download PDF</a>\n".format(os.path.basename(os.path.dirname(src_path)), filename))
            f_out.write("</object>\n")

        elif ext in [".fa", ".fasta", ".txt"]:
            with open(src_path, "r") as f_in:
                content = f_in.read()
            f_out.write("<pre style='font-size:10px; font-family: monospace;'>\n")
            f_out.write(content)
            f_out.write("</pre>\n")

        else:
            f_out.write(f"<p>Format {ext} is not supported.</p>\n")

        f_out.write("</body></html>\n")

def parse_newick(newick):
    """
    Simple Newick parser.
    Returns a nested tuple representation of the tree.
    """
    newick = newick.strip().rstrip(';')
    stack = []
    current = []
    label = ''

    for char in newick:
        if char == '(':
            stack.append(current)
            current = []
        elif char == ')':
            if label:
                current.append(label.strip())
                label = ''
            parent = stack.pop()
            parent.append(current)
            current = parent
        elif char == ',':
            if label:
                current.append(label.strip())
                label = ''
        else:
            label += char

    if label:
        current.append(label.strip())

    return current[0] if current else []

def count_leaves(tree):
    """
    Count the number of leaves in a Newick tree.
    """
    if isinstance(tree, str):
        return 1
    return sum(count_leaves(child) for child in tree)

def draw_tree_rectangular(tree, ax, x=0, y=0, dx=1.5, dy=1.5, fontsize=10):
    """
    Draw a Newick tree in a rectangular layout.
    """
    if isinstance(tree, str):
        ax.text(x, y, tree, verticalalignment='center', fontsize=fontsize)
        return x, y

    child_coords = []
    for i, child in enumerate(tree):
        cx, cy = draw_tree_rectangular(child, ax, x + dx, y - i * dy, dx, dy, fontsize)
        child_coords.append((cx, cy))

    avg_y = sum(cy for _, cy in child_coords) / len(child_coords)
    ax.plot([x + dx, x], [avg_y, avg_y], color='black')
    for cx, cy in child_coords:
        ax.plot([x, cx], [avg_y, cy], color='black')

    return x, avg_y

def draw_tree_polar(tree, ax, radius=1.0, angle=0.0, spread=2 * math.pi):
    """
    Draw a Newick tree in a polar layout.
    """
    def draw_node(node, r, theta, spread):
        if isinstance(node, str):
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            ax.text(x, y, node, ha='center', va='center', fontsize=8)
            return [(x, y)]

        n = len(node)
        coords = []
        for i, child in enumerate(node):
            angle_offset = (i + 0.5) * (spread / n) - (spread / 2)
            child_theta = theta + angle_offset
            child_coords = draw_node(child, r + radius, child_theta, spread / n)
            coords.extend(child_coords)

        parent_x = r * math.cos(theta)
        parent_y = r * math.sin(theta)
        for x, y in coords:
            ax.plot([parent_x, x], [parent_y, y], color='black')

        return [(parent_x, parent_y)]

    ax.set_aspect('equal')
    ax.axis('off')
    draw_node(tree, radius, angle, spread)


def draw_tree_radial(tree, ax, radius=1.0):
    """
    Draw a Newick tree in a radial layout.
    """
    def layout(node, depth=0, angle_start=0.0, angle_end=2 * math.pi):
        if isinstance(node, str):
            angle = (angle_start + angle_end) / 2
            x = depth * radius * math.cos(angle)
            y = depth * radius * math.sin(angle)
            ax.text(x, y, node, ha='center', va='center', fontsize=8)
            return [(x, y)]

        step = (angle_end - angle_start) / len(node)
        children_coords = []
        for i, child in enumerate(node):
            child_coords = layout(child, depth + 1, angle_start + i * step, angle_start + (i + 1) * step)
            children_coords.extend(child_coords)

        parent_angle = (angle_start + angle_end) / 2
        px = depth * radius * math.cos(parent_angle)
        py = depth * radius * math.sin(parent_angle)
        for cx, cy in children_coords:
            ax.plot([px, cx], [py, cy], color='black')

        return [(px, py)]

    ax.set_aspect('equal')
    ax.axis('off')
    layout(tree)

def plot_all_trees(result_folder: str, output_folder: str, bait_file: str, candidate_file: str,
                   format: str = None, layout: str = "Polar", 
                   cand_color: str = "green", bait_color: str = "red"):
    """
    Searches for Newick trees in the tree order and plots them as SVG and PDF/PNG files.
    """

    os.makedirs(output_folder, exist_ok=True)
    tree_files = [f for f in os.listdir(result_folder) if f.endswith(".tre")]

    bait_ids = set(read_fasta(bait_file).keys())
    cand_ids = set(read_fasta(candidate_file).keys())
    
    if not bait_ids:
        print(f"Warning: No bait IDs found from {bait_file}. Bait coloring will not be applied.")
    if not cand_ids:
        print(f"Warning: No candidate IDs found from {candidate_file}. Candidate coloring will not be applied.")

    for tree_file in tree_files:
        tree_path = os.path.join(result_folder, tree_file)
        try:
            with open(tree_path, 'r') as f:
                newick = f.read().strip()

            tree = Tree(newick, format=1)
            ts = TreeStyle()
            ts.show_leaf_name = False

            layout_mode = layout.lower()
            if layout_mode == "rectangular":
                ts.mode = "r"
            elif layout_mode == "radial":
                ts.mode = "c"
            else:  # Default to polar
                ts.mode = "c"

            leaf_count = len(tree.get_leaves())
            width = max(1200, int(leaf_count * 25))
            dpi = 300

            for node in tree.traverse():
                if node.is_leaf():
                    nstyle = NodeStyle()
                    nstyle["size"] = 0
                    node.set_style(nstyle)
                    
                    color = "black"
                    if node.name in bait_ids:
                        color = bait_color
                    elif node.name in cand_ids:
                        color = cand_color
                    
                    face = TextFace(node.name, fsize=10, fgcolor=color)
                    node.add_face(face, column=0, position="branch-right")

            base_name = os.path.splitext(tree_file)[0]

            # Always save SVG file
            svg_path = os.path.join(output_folder, f"{base_name}.svg")
            tree.render(svg_path, w=width, units="px", tree_style=ts)
            #print(f"Tree SVG saved to: {svg_path} ({leaf_count} leaves)")

            # Optional saving as PNG or PDF
            if format in ["png", "pdf"]:
                output_path = os.path.join(output_folder, f"{base_name}.{format}")
                tree.render(output_path, w=width, units="px", dpi=dpi, tree_style=ts)
                #print(f"Tree {format.upper()} saved to: {output_path}")

        except Exception as e:
            print(f"Error processing {tree_file}: {e}")

def generate_heatmaps_from_expression_data(txt_folder, heatmap_folder, format = None):
    """
    Creates heatmaps from expression data and saves them as SVG and PDF/PNG files.
    """

    filenames = [
        "07_filtered_expression_matrix.txt",
        "07_metadata_expression_matrix.txt",
        "07_metadata_only_expression_matrix.txt"
    ]

    for filename in filenames:
        file_path = os.path.join(txt_folder, filename)
        if os.path.isfile(file_path):
            try:
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                if df.empty:
                    print(f"File {filename} is empty. Heatmap creation will be skipped.")
                    continue

                n_rows, n_cols = df.shape
                base_fontsize = 10
                x_fontsize = max(4, min(12, base_fontsize * 40 / max(n_cols, 1)))
                y_fontsize = max(4, min(12, base_fontsize * 40 / max(n_rows, 1)))

                fig_width = min(150, max(10, n_cols * 0.3))
                fig_height = min(30, max(6, n_rows * 0.3))

                plt.figure(figsize=(fig_width, fig_height))
                ax = sns.heatmap(df, cmap="viridis", cbar=True)

                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=x_fontsize)
                ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=y_fontsize)
                plt.title(f"Heatmap: {filename}", fontsize=14)
                plt.tight_layout()

                base_name = os.path.splitext(filename)[0]

                # Always save SVG
                svg_path = os.path.join(heatmap_folder, f"{base_name}.svg")
                plt.savefig(svg_path, format="svg", dpi=300)
                print(f"Heatmap SVG saved at: {svg_path}")

                # Optional save as PNG or PDF
                if format in ["png", "pdf"]:
                    output_path = os.path.join(heatmap_folder, f"{base_name}.{format}")
                    plt.savefig(output_path, format=format, dpi=300)
                    #print(f"Heatmap {format.upper()} saved at: {output_path}")

                plt.close()
            except Exception as e:
                print(f"Error processing {filename}: {e}")
        else:
            print(f"File not found: {filename} – will be skipped.")

#endregion

# region Main
def main():
    print("\n")
    print("CYP_Annotator started." + "\n")
    print("Step 0: Preprocessing input data.")
    
    # Parse and validate command-line arguments
    args = parser.parse_args()
    args = validate_args(args)

    # Check required tools
    tools = check_required_tools()
    # use_blast = tools['BLAST']
    # use_hmmer = tools['HMMER']
    mafft_available = tools['MAFFT']
    muscle_available = tools['MUSCLE']

    # Handle file collection CSV if provided
    if args.file_collection:
        csv_paths = read_file_collection_csv(args.file_collection)
        if csv_paths.get('Data'):
            args.data = csv_paths['Data']
        else:
            # Assign paths from CSV if Data path not provided
            args.baits = csv_paths.get('Baits')
            args.baits_info = csv_paths.get('Baits_Info')
            args.hmm_domains = csv_paths.get('HMM_Domains')
            args.hmm_motifs = csv_paths.get('HMM_Motifs')
            args.protein_motifs = csv_paths.get('Protein_Motifs')
            args.expression = csv_paths.get('Expression')
            args.metadata = csv_paths.get('Metadata')

    # A SINGLE CALL to collect and process all files
    processed_paths = prepare_input_files(args)

    # Extract all paths from processed_paths dictionary
    baits_with_outgroup_path = processed_paths['baits_with_outgroup_path']
    baits_path = processed_paths['baits_no_outgroup_path']
    baits_info = processed_paths['baits_info']
    info_headers = processed_paths['baits_info_headers']
    transcript_baits_path = processed_paths['transcript_baits_path']
    protein_baits_path = processed_paths['protein_baits_path']
    literature_baits_path = processed_paths['literature_baits_path']
    outgroup_path = processed_paths['outgroup_path']
    subjects_paths = processed_paths['subject_files_paths']
    protein_motifs_path = processed_paths.get('protein_motifs_path')
    hmm_motifs_path = processed_paths.get('hmm_motifs_path')
    hmm_domains_path = processed_paths.get('hmm_domains_path')
    #baits_info_path = processed_paths.get('baits_info_path')
    expression_matrix_path = processed_paths.get('expression_path')
    metadata_path = processed_paths.get('metadata_path')
    
    # Set output folder
    output_folder = args.output_folder if args.output_folder is not None else "Output"
    processed_input_folder = args.processed_input_folder if args.processed_input_folder is not None else os.path.join(output_folder, "Processed_Input/")
    
    print("Step 0 completed. " + "\n")

    # --- separated analyses for each subject --- #             
    for jidx, subject in enumerate(subjects_paths):
        
        job_ID = Path(subject).stem
        job_output_folder = Path(output_folder)
        
        if len(subjects_paths) > 1:
            job_output_folder = job_output_folder / f"{jidx:05d}_{job_ID}"
        
        job_output_folder.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing job {jidx+1}")

        # Create folders
        result_folder = os.path.join(job_output_folder, "RESULTS/")
        supplement_folder = os.path.join(job_output_folder, "SUPPLEMENTS/")
        
        # Create subfolders
        txt_folder = os.path.join(result_folder, "TXT/")
        fasta_folder = os.path.join(result_folder, "FASTA/")
        tree_folder = os.path.join(result_folder, "TREES/")
        html_folder = os.path.join(result_folder, "HTML/")
        heatmap_folder = os.path.join(result_folder, "HEATMAPS/")

        # Make folders
        os.makedirs(result_folder, exist_ok=True)
        os.makedirs(supplement_folder, exist_ok=True)
        os.makedirs(txt_folder, exist_ok=True)
        os.makedirs(fasta_folder, exist_ok=True)
        os.makedirs(tree_folder, exist_ok=True)
        os.makedirs(html_folder, exist_ok=True)
        if expression_matrix_path:
            os.makedirs(heatmap_folder, exist_ok=True)

        tree_counter = 0
        
        # --- 00 Write baits_info to TXT --- #
        bait_groups = {
            seq_id: subdict.get("Evidence", "-")
            for seq_id, subdict in baits_info.items()
            }
        baits_info_file = os.path.join(txt_folder, "00_baits_info.txt")
        baits_info_to_txt(baits_info, info_headers, baits_info_file)

        # baits_info_file = txt_folder / "00_baits_info.txt"
        # baits_info_to_txt(baits_info, info_headers, str(baits_info_file))

        # --- 00 Generate documentation file (simplified call) --- #
        doc_file = txt_folder + args.name + "00_documentation.txt"
        generate_documentation_file(doc_file, baits_with_outgroup_path, baits_path, baits_info_file, hmm_domains_path, hmm_motifs_path, protein_motifs_path, 
                                    job_output_folder, subjects_paths, args.use_hmmer, args.mode_aln ,args.mode_tree, args.blastp, args.makeblastdb, args.hmmsearch, 
                                    args.cpu_max, args.mafft, args.muscle, args.raxml, args.fasttree,
                                    args.bitcutp, args.simcutp, args.poscutp, args.lencutp,
                                    args.minscore, args.numneighbours, args.neighbourdist, args.minneighbours, args.paralogdist,
                                    args.parallel, args.num_process_candidates, args.name, args.trim_names, args.collapse)
        
        start_time = datetime.datetime.now()

        # --- 01 Search initial candidates --- #
        print("Step 1: Search initial candidates")
        seq_search_result_file = os.path.join(supplement_folder, "01_seq_search_results.txt")
        blast_analyze_folder = os.path.join(supplement_folder, "01_blast_seq_search_analysis/")
        os.makedirs(blast_analyze_folder, exist_ok=True)

        hmm_file = args.hmm if args.hmm else os.path.join(blast_analyze_folder, "bait.hmm")
        # Create HMM file if not provided
        if args.hmm is None:
            bait_source = literature_baits_path if literature_baits_path else baits_path
            if not os.path.isfile(hmm_file):
                try:
                    hmm_build(bait_source, hmm_file)
                except Exception as e:
                    print(f"Error creating HMM file: {str(e)}")
                    raise SystemExit("Skript was canceled, as no HMM file was found.")

        search_files_exist = all ([
            os.path.isfile(seq_search_result_file) and os.path.getsize(seq_search_result_file) > 0,
            os.path.isdir(blast_analyze_folder)
        ])

        # Skip search if result file already exists
        if not search_files_exist:
            seq_search_results = {}
            use_hmmer = args.use_hmmer.upper() in ("Y", "YES")

            # --- HMMER-Path ---
            if use_hmmer:
                print("Running HMMER search for CYP candidates")
                hmm_result_file = os.path.join(supplement_folder, "hmm_results.txt")
                if not os.path.isfile(hmm_result_file):
                    run_hmmsearch(
                        hmm_file=hmm_file,
                        subject_file=subject,
                        output_file=hmm_result_file
                    )
                seq_search_results = load_hmmsearch_results(hmm_result_file)

            # --- BLAST-Path ---
            else:
                print("Running BLAST search for CYP candidates")
                blast_result_file = os.path.join(supplement_folder, "blast_results.txt")
                if not os.path.isfile(blast_result_file):
                    run_blast_search(
                        query_file=baits_path,
                        subject_file=subject,
                        output_file=blast_result_file,
                        blast_db_folder=blast_analyze_folder,
                        cpu_max=args.cpu_max
                    )
                seq_search_results = load_blast_results(
                    blast_result_file,
                    args.simcutp,
                    args.poscutp,
                    args.lencutp,
                    args.bitcutp
                )

            # Write out results
            with open(seq_search_result_file, 'w') as f:
                for seq_id in seq_search_results.keys():
                    f.write(f"{seq_id}\n")

        else:
            #print(f"Load results from {seq_search_result_file}")
            with open(seq_search_result_file) as f:
                seq_search_results = {line.strip(): None for line in f if line.strip()}

        # Load sequences
        subject_sequences = read_fasta(subject)
        subject_name_mapping_table = load_name_mapping_table(Path(processed_input_folder) / f"{job_ID}_mapping.txt")

        # Save candidates as FASTA
        candidate_file = os.path.join(fasta_folder, f"{args.name}01_initial_candidates.fasta")
        candidate_sequences = {seq_id: subject_sequences[seq_id] for seq_id in seq_search_results.keys() if seq_id in subject_sequences}
        cand_count = 0
        with open(candidate_file, "w") as out:
            for seq_id in seq_search_results.keys():
                if seq_id in subject_sequences:
                    cand_count += 1
                    out.write(f">{seq_id}\n{subject_sequences[seq_id]}\n")

        print(f"Step 1 completed. {cand_count} candidate sequences found.")

        # --- 02 construct phylogenetic tree and analyze tree file --- #
        print("Step 2: Assigning candidates to in and out groups")
        inout_output_folder = os.path.join(supplement_folder, f"{args.name}02_in_out_CYP_analysis_trees/")
        os.makedirs(inout_output_folder, exist_ok=True)


        in_list, out_list = create_in_out_anno(baits_path, outgroup_path)
        aln_candidate_file = "alignment_candidates.fasta"
        aln_input_file = "alignment_input.fasta"
        aln_file = "alignment_input.fasta.aln"
        cln_aln_file = "alignment_input.fasta.aln.cln"
        
        # HMM motif check
        tree_hmmsearch_results_file = os.path.join(inout_output_folder, "02_hmmsearch_results.txt")
        tree_hmmsearch_results = {}

        if not (os.path.isfile(tree_hmmsearch_results_file) and os.path.getsize(tree_hmmsearch_results_file) > 0):
            run_hmmsearch(
                hmm_file=hmm_file,
                subject_file=candidate_file,
                output_file=tree_hmmsearch_results_file
            )
        tree_hmmsearch_results = load_hmmsearch_results(tree_hmmsearch_results_file)
        
        # --- 02 first classification --- #
        filtered_members_file_f = os.path.join(inout_output_folder, f"{args.name}first_ingroup_CYPs.fasta")
        tmp_result_table_f = os.path.join(inout_output_folder, f"{args.name}first_in_out_CYP_analysis_results.txt")

        first_classification_files_exist = all ([
            os.path.isfile(filtered_members_file_f) and os.path.getsize(filtered_members_file_f) > 0,
            os.path.isfile(tmp_result_table_f) and os.path.getsize(tmp_result_table_f) > 0
        ])

        if not first_classification_files_exist:
            cyp_classification = {}
            # sys.stdout.write(f"Number of ingroup CYP baits: {len(in_list)}\n")
            # sys.stdout.write(f"Number of outgroup CYP baits: {len(out_list)}\n")
            sys.stdout.flush()

            # construct phylogenetic trees
            tree_files = parallel_tree_constructor(
                args.num_process_candidates, inout_output_folder, aln_candidate_file, aln_input_file,
                aln_file, cln_aln_file, baits_with_outgroup_path, candidate_file, "first_", "",
                args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree, args.iqtree,
                args.cpu_max, args.parallel
            )

            for tree_file in tree_files:
                # Replace spaces in leaf names before processing
                modified_tree_file = replace_spaces_in_tree(tree_file)

                classification = split_into_ingroup_and_outgroup(
                    modified_tree_file, in_list, out_list, args.numneighbours,
                    args.neighbourdist, args.minneighbours, tree_hmmsearch_results
                )
                cyp_classification.update(classification)

            # Prepare data for writing
            result_lines = []
            fasta_content = {}
            header = "ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n"
            
            candidate_order = sorted(cyp_classification.keys())
            for candidate in candidate_order:
                result_lines.append("\t".join(map(str, [
                    cyp_classification[candidate]["tree"],
                    subject_name_mapping_table[candidate],
                    candidate,
                    cyp_classification[candidate]["hmm"],
                    cyp_classification[candidate]["score"],
                    cyp_classification[candidate]["in"],
                    cyp_classification[candidate]["out"],
                ])))
                
                if cyp_classification[candidate]["score"] > args.minscore:
                    fasta_content[candidate] = subject_sequences[candidate]

            # Write output files
            with open(filtered_members_file_f, "w") as out, open(tmp_result_table_f, "w") as out2:
                out2.write("ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n")
                candidate_order = sorted(cyp_classification.keys())
                for candidate in candidate_order:
                    out2.write("\t".join(map(str, [
                        cyp_classification[candidate]["tree"],
                        subject_name_mapping_table[ candidate ],
                        candidate,
                        cyp_classification[candidate]["hmm"],
                        cyp_classification[candidate]["score"],
                        cyp_classification[candidate]["in"],
                        cyp_classification[candidate]["out"],
                    ])) + "\n")

                    if cyp_classification[ candidate ]["score"] > args.minscore and ( candidate in tree_hmmsearch_results or not args.filterdomain.upper() in ["YES", "Y"]):                            
                            out.write( ">" + candidate + "\n" + subject_sequences[ candidate ] + "\n" )

        # --- 02 second classification --- #
        filtered_members_file_s = os.path.join(inout_output_folder, f"{args.name}first_ingroup_CYPs.fasta")
        tmp_result_table_s = os.path.join(inout_output_folder, f"{args.name}first_in_out_CYP_analysis_results.txt")
        filtered_members_file = os.path.join(fasta_folder, f"{args.name}02_ingroup_CYPs.fasta")
        tmp_result_table = os.path.join(txt_folder, f"{args.name}02_in_out_CYP_analysis_results.txt")

        second_classification_files_exist = all ([
            os.path.isfile(filtered_members_file_s) and os.path.getsize(filtered_members_file_s) > 0,
            os.path.isfile(tmp_result_table_s) and os.path.getsize(tmp_result_table_s) > 0,
            os.path.isfile(filtered_members_file) and os.path.getsize(filtered_members_file) > 0,
            os.path.isfile(tmp_result_table) and os.path.getsize(tmp_result_table) > 0
        ])

        if not second_classification_files_exist:
            cyp_classification = {}
            #sys.stdout.write(f"Number of ingroup CYP baits: {len(in_list)}\n")
            #sys.stdout.write(f"Number of outgroup CYP baits: {len(out_list)}\n")
            sys.stdout.flush()

            # construct phylogenetic trees
            tree_files = parallel_tree_constructor(
                args.num_process_candidates, inout_output_folder, aln_candidate_file, aln_input_file,
                aln_file, cln_aln_file, baits_with_outgroup_path, candidate_file, "second_", "",
                args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree, args.iqtree,
                args.cpu_max, args.parallel
            )

            for tree_file in tree_files:
                # Replace spaces in leaf names before processing
                modified_tree_file = replace_spaces_in_tree(tree_file)
                
                classification = split_into_ingroup_and_outgroup(
                    modified_tree_file, in_list, out_list, args.numneighbours,
                    args.neighbourdist, args.minneighbours, tree_hmmsearch_results
                )
                cyp_classification.update(classification)

            # Prepare data for writing
            result_lines = []
            fasta_content = {}
            header = "ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n"
            
            candidate_order = sorted(cyp_classification.keys())
            for candidate in candidate_order:
                result_lines.append("\t".join(map(str, [
                    cyp_classification[candidate]["tree"],
                    subject_name_mapping_table[candidate],
                    candidate,
                    cyp_classification[candidate]["hmm"],
                    cyp_classification[candidate]["score"],
                    cyp_classification[candidate]["in"],
                    cyp_classification[candidate]["out"],
                ])))
                
                if cyp_classification[candidate]["score"] > args.minscore:
                    fasta_content[candidate] = subject_sequences[candidate]

            with open(filtered_members_file_s, "w") as out, open(tmp_result_table_s, "w") as out2:
                out2.write("ClassificationTree\tOriginalID\tCleanID\tHMMmotif\tScore\tIngroupMatches\tOutgroupMatches\n")
                candidate_order = sorted(cyp_classification.keys())
                for candidate in candidate_order:
                    out2.write("\t".join(map(str, [
                        cyp_classification[candidate]["tree"],
                        subject_name_mapping_table[ candidate ],
                        candidate,
                        cyp_classification[candidate]["hmm"],
                        cyp_classification[candidate]["score"],
                        cyp_classification[candidate]["in"],
                        cyp_classification[candidate]["out"],
                    ])) + "\n")

                    if cyp_classification[candidate]["score"] > args.minscore and (candidate in tree_hmmsearch_results or not args.filterdomain.upper() in ["YES", "Y"]):
                        out.write(f">{candidate}\n{subject_sequences[candidate]}\n")

            # Copy final files to results folder
            shutil.copyfile(filtered_members_file_s, filtered_members_file)
            shutil.copyfile(tmp_result_table_s, tmp_result_table)

        else:
            cyp_classification = load_in_out_classification_file(tmp_result_table)

        filtered_members = read_fasta(filtered_members_file)
        if len(filtered_members) < 1:
            sys.exit(f"ERROR: no CYPs detected.")

        print(f"Step 2 completed. {len(filtered_members)} ingroup candidates selected.")

        # --- 03 construct a filtered tree and find closest reference for baits --- #
        print("Steps 3 and 4: Construct phylogenetic trees and assign orthologs to candidates.")
        prefix=args.ortholog_prefix
        first_prefix=f"03_{args.ortholog_prefix}"
        second_prefix=f"04_{args.ortholog_prefix}"
        filtered_group_around_ref_file = os.path.join(txt_folder, f"{second_prefix}_candidates_group_around_bait.txt")
        filtered_ref_mapping_file = os.path.join(txt_folder, f"{second_prefix}_candidate_2_bait_mapping_file.txt")
        final_members_file = Path(fasta_folder) / f"{second_prefix}_candidate.fasta"
        linneage_specific_file = Path(fasta_folder) / f"{second_prefix}_linneage_specific_sequences.fasta"

        ortholog_files_exist = all ([
            os.path.isfile(filtered_group_around_ref_file) and os.path.getsize(filtered_group_around_ref_file) > 0,
            os.path.isfile(filtered_ref_mapping_file) and os.path.getsize(filtered_ref_mapping_file) > 0,
            os.path.isfile(final_members_file) and os.path.getsize(final_members_file) > 0
        ])

        if not (os.path.isdir(tree_folder) and any(f.startswith("03") for f in os.listdir(tree_folder))) \
        or not ortholog_files_exist:

            filtered_tree_path, final_tree_path = perform_candidate_tree_analysis(
                prefix=prefix,
                first_prefix=first_prefix,
                second_prefix=second_prefix,
                bait_fasta_path=baits_path,
                clean_members=filtered_members,
                clean_members_file=filtered_members_file,
                baits_info=baits_info,
                subject_name_mapping_table=subject_name_mapping_table,
                bait_groups=bait_groups,
                args=args,
                supplement_folder=supplement_folder,
                tree_folder=tree_folder,
                group_around_ref_file=filtered_group_around_ref_file,
                ref_mapping_file=filtered_ref_mapping_file,
                filtered_fasta_file=final_members_file,
                linneage_specific_fasta = linneage_specific_file
            )
        else:
            filtered_tree_path = os.path.join(tree_folder, f"{args.name}{first_prefix}_first_0_FastTree_tree.tre")
            final_tree_path = os.path.join(tree_folder, f"{args.name}{second_prefix}_final_0_FastTree_tree.tre")

        final_members = read_fasta(final_members_file)
        tree_counter += 2

        # --- 04 Optional: Individual analysis with selected bait sequences --- #
        if args.individual_tree.upper() in ["YES", "Y"]:
            if args.individual_ortholog_prefix is None:
                args.individual_ortholog_prefix = args.bait_keyword.split()[0]
            prefix=args.individual_ortholog_prefix
            first_prefix=f"03.1_{args.individual_ortholog_prefix}"
            second_prefix=f"04.1_{args.individual_ortholog_prefix}"
            ind_group_around_ref_file = os.path.join(txt_folder, f"{second_prefix}_candidates_group_around_bait.txt")
            ind_ref_mapping_file = os.path.join(txt_folder, f"{second_prefix}_candidate_2_bait_mapping_file.txt")
            ind_filtered_fasta_file = Path(fasta_folder) / f"{second_prefix}_candidate.fasta"

            individual_ortholog_files_exist = all ([
                os.path.isfile(ind_group_around_ref_file) and os.path.getsize(ind_group_around_ref_file) > 0,
                os.path.isfile(ind_ref_mapping_file) and os.path.getsize(ind_ref_mapping_file) > 0,
                os.path.isfile(ind_filtered_fasta_file) and os.path.getsize(ind_filtered_fasta_file) > 0
            ])

            if not (os.path.isdir(tree_folder) and any(f.startswith("04") for f in os.listdir(tree_folder))) \
            or not individual_ortholog_files_exist:
                filtered_bait_ids = filter_baits_by_info(baits_info, args.bait_column, args.bait_keyword)
                reduced_bait_fasta = os.path.join(fasta_folder, f"{args.name}04_filtered_baits.fasta")
                bait_sequences = read_fasta(baits_path)
                selected_seqs = {k: v for k, v in bait_sequences.items() if k in filtered_bait_ids}
                write_fasta(selected_seqs, Path(reduced_bait_fasta))

                ind_filtered_tree_path, ind_final_tree_path = perform_candidate_tree_analysis(
                    prefix=prefix,
                    first_prefix=first_prefix,
                    second_prefix=second_prefix,
                    bait_fasta_path=reduced_bait_fasta,
                    clean_members=filtered_members,
                    clean_members_file=filtered_members_file,
                    baits_info=baits_info,
                    subject_name_mapping_table=subject_name_mapping_table,
                    bait_groups=bait_groups,
                    args=args,
                    supplement_folder=supplement_folder,
                    tree_folder=tree_folder,
                    group_around_ref_file=ind_group_around_ref_file,
                    ref_mapping_file=ind_ref_mapping_file,
                    filtered_fasta_file=ind_filtered_fasta_file,
                    linneage_specific_fasta = ""
                )
            else:
                ind_filtered_tree_path = os.path.join(tree_folder, f"{args.name}{first_prefix}_first_0_FastTree_tree.tre")
                ind_final_tree_path = os.path.join(tree_folder, f"{args.name}{second_prefix}_final_0_FastTree_tree.tre")

            tree_counter += 2

        print("Steps 3 and 4 completed.")

        # --- 05 check for protein domains --- #
        print("Step 5: Check for protein domains.")
        domain_check_file_summary = os.path.join(txt_folder, f"{args.name}05_domain_check.txt")
        best_domain_match_file = os.path.join(txt_folder, f"{args.name}05_best_domain_match.txt")
        domain_check_hmm_result = os.path.join(supplement_folder, "05_domain_hmmsearch.txt")
        domain_check_hmm_output = os.path.join(supplement_folder, "05_domain_hmmsearch_output.txt")

        domain_files_exist_and_not_empty = all([
            os.path.isfile(domain_check_file_summary) and os.path.getsize(domain_check_file_summary) > 0,
            os.path.isfile(best_domain_match_file) and os.path.getsize(best_domain_match_file) > 0,
            os.path.isfile(domain_check_hmm_result) and os.path.getsize(domain_check_hmm_result) > 0,
            os.path.isfile(domain_check_hmm_output) and os.path.getsize(domain_check_hmm_output) > 0,
        ])

        if not domain_files_exist_and_not_empty:
            if hmm_domains_path is not None:
                domain_check_results, best_domain_match, domain_names = domain_check(
                    final_members_file,
                    hmm_domains_path,
                    domain_check_hmm_result,
                    domain_check_hmm_output,
                    args.hmmsearch,
                    args.domain_Score
                )

                ref_family = load_ortholog_family(filtered_ref_mapping_file)
                with open(domain_check_file_summary, "w") as out1:
                    out1.write("RefMember\t" + "\t".join(domain_names) + "\n")
                    candidates = list(sorted(final_members.keys()))
                    for candidate in candidates:
                        new_line = [subject_name_mapping_table[candidate]]
                        for dom in domain_names:
                            new_line.append("YES" if domain_check_results[candidate].get(dom, "") else "NO")
                        out1.write("\t".join(new_line) + "\n")

                with open(best_domain_match_file, "w") as out2:
                    out2.write("RefMember\tBestDomain\tOrtholog_Family\n")
                    for candidate in sorted(final_members.keys()):
                        best_dom = best_domain_match.get(candidate, "-")
                        candidate_name = subject_name_mapping_table[candidate]
                        ortholog_family = ref_family.get(candidate_name, "None")
                        out2.write(f"{candidate_name}\t{best_dom}\t{ortholog_family}\n")
                        
        print("Step 5 completed.")

        # --- 06 check for common protein motifs --- #
        print("Step 6: Check for common protein motifs.")
        static_summary_file = os.path.join(txt_folder, f"{args.name}06_static_motif_check_summary.txt")
        static_results_file = os.path.join(txt_folder, f"{args.name}06_static_motif_check_results.txt")
        hmm_summary_file = os.path.join(txt_folder, f"{args.name}06_hmm_motif_check_summary.txt")
        hmm_results_file = os.path.join(txt_folder, f"{args.name}06_hmm_motif_check_results.txt")

        motif_check_hmm_result = os.path.join(supplement_folder, "06_motif_hmmsearch.txt")
        motif_check_hmm_output = os.path.join(supplement_folder, "06_motif_hmmsearch_output.txt")

        static_motif_files_exist_and_not_empty = all([
            os.path.isfile(static_summary_file) and os.path.getsize(static_summary_file) > 0,
            os.path.isfile(static_results_file) and os.path.getsize(static_results_file) > 0,
        ])

        if not static_motif_files_exist_and_not_empty:
            # --- Static motif check ---
            static_motif_results = static_motif_check(final_members_file, protein_motifs_path)
            static_data, static_motif_names = static_motif_results

            # Claculate raw motifs lenghts
            raw_static_lengths = {}
            with open(protein_motifs_path, 'r') as f:
                for line in f.readlines()[1:]:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            name, pattern = parts[0], parts[1]
                            mlen = len(pattern.replace("(", "").replace(")", "").replace("/", ""))
                            raw_static_lengths[name] = mlen
            static_flanks = {motif: (0, 0) for motif in static_motif_names}

            write_motif_output(static_summary_file, static_results_file, static_data, static_motif_names,
                               raw_static_lengths, static_flanks, final_members, subject_name_mapping_table)

        hmm_motif_files_exist_and_not_empty = all([
            os.path.isfile(hmm_summary_file) and os.path.getsize(hmm_summary_file) > 0,
            os.path.isfile(hmm_results_file) and os.path.getsize(hmm_results_file) > 0,
            os.path.isfile(motif_check_hmm_result) and os.path.getsize(motif_check_hmm_result) > 0
        ])

        if not hmm_motif_files_exist_and_not_empty:
        # --- HMM motif check ---
            if hmm_motifs_path is None and protein_motifs_path is not None:
                CYP_source = literature_baits_path if literature_baits_path else baits_path
                hmm_motifs_path, raw_hmm_lengths = create_motif_from_baits(
                    protein_motifs_path, CYP_source, supplement_folder, 
                    mafft_available, muscle_available, flanklength=0
                )
            else:
                raw_hmm_lengths = {}

            hmm_motif_names, hmm_lengths = parse_motif_names_lengths(hmm_motifs_path)
            target_hmm_len = max(hmm_lengths.values())
            target_hmm_len += 4 if target_hmm_len % 2 == 0 else 5

            hmm_flanks = {}
            for motif, length in hmm_lengths.items():
                pad_total = target_hmm_len - length
                flank_left = pad_total // 2
                flank_right = pad_total - flank_left
                hmm_flanks[motif] = (flank_left, flank_right)

            hmm_motif_results = motif_check(
                final_members_file, hmm_motifs_path,
                motif_check_hmm_result, motif_check_hmm_output,
                args.hmmsearch, args.motif_cEvalue
            )
            hmm_data = hmm_motif_results[0]

            write_motif_output(hmm_summary_file, hmm_results_file, hmm_data, hmm_motif_names,
                               hmm_lengths, hmm_flanks, final_members, subject_name_mapping_table)

        # --- Create tree filtered by motifs ---
        prefix="Motifs"
        first_prefix=f"07.1_{args.ortholog_prefix}"
        second_prefix=f"07.2_{args.ortholog_prefix}"
        motifs_group_around_ref_file = os.path.join(txt_folder, f"{second_prefix}_candidates_group_around_bait.txt")
        motifs_ref_mapping_file = os.path.join(txt_folder, f"{second_prefix}_candidate_2_bait_mapping_file.txt")
        motifs_members_file = Path(fasta_folder) / f"{second_prefix}_candidate.fasta"

        motif_ortholog_files_exist = all ([
            os.path.isfile(motifs_group_around_ref_file) and os.path.getsize(motifs_group_around_ref_file) > 0,
            os.path.isfile(motifs_ref_mapping_file) and os.path.getsize(motifs_ref_mapping_file) > 0,
            os.path.isfile(motifs_members_file) and os.path.getsize(motifs_members_file) > 0
        ])

        if not (os.path.isdir(tree_folder) and any(f.startswith("07") for f in os.listdir(tree_folder))) \
        or not motif_ortholog_files_exist:
        
            seq_id_list = []
            motif_members = {}

            try:
                with open(static_summary_file, mode='r', newline='') as infile:
                    # CSV-Reader für eine Tab-separierte Datei erstellen
                    reader = csv.reader(infile, delimiter='\t')

                    # 2. Kopfzeile überspringen
                    header = next(reader, None)

                    # 3. Über die Datenzeilen iterieren
                    for row in reader:
                        # Sicherheitsabfrage, ob die Zeile genügend Spalten hat
                        if len(row) > 6:
                            seq_id = row[0]  # Erste Spalte
                            heme_motif = row[6] # Siebte Spalte
                            if heme_motif == 'YES':
                                seq_id_list.append(seq_id)
                                if seq_id in final_members.keys():
                                    motif_members[seq_id] = final_members[seq_id]


            except FileNotFoundError:
                print(f"Error: '{static_summary_file}' not found.")

            motif_fasta_file = os.path.join(supplement_folder, f"motifs_to_tree.fasta")
            write_fasta(motif_members, Path(motif_fasta_file))

            motifs_filtered_tree_path, motifs_final_tree_path = perform_candidate_tree_analysis(
                prefix=prefix,
                first_prefix=first_prefix,
                second_prefix=second_prefix,
                bait_fasta_path=baits_path,
                clean_members=motif_members,
                clean_members_file=motif_fasta_file,
                baits_info=baits_info,
                subject_name_mapping_table=subject_name_mapping_table,
                bait_groups=bait_groups,
                args=args,
                supplement_folder=supplement_folder,
                tree_folder=tree_folder,
                group_around_ref_file=motifs_group_around_ref_file,
                ref_mapping_file=motifs_ref_mapping_file,
                filtered_fasta_file=motifs_members_file,
                linneage_specific_fasta = ""
            )
        else:
            motifs_filtered_tree_path = os.path.join(tree_folder, f"{prefix}{first_prefix}_first_0_FastTree_tree.tre")
            motifs_final_tree_path = os.path.join(tree_folder, f"{prefix}{second_prefix}_final_0_FastTree_tree.tre")

        tree_counter += 2

        print("Step 6 completed.")

        # --- 07 filter expression matrix for candidates --- # 
        if expression_matrix_path is None:
            print("Step 7: No expression data provided. Skipping procession of expression data.")
            # highest_expression_file = os.path.join(txt_folder, f"{args.name}07_highest_expression_placeholder.txt")
            # with open(highest_expression_file, 'w') as f:
            #     f.write("\n")
            highest_expression_file = None
        else:
            print("Step 7: Process expression data.")
            filtered_expression_file = os.path.join(txt_folder, f"{args.name}07_filtered_expression_matrix.txt")

            if not (os.path.isfile(filtered_expression_file) and os.path.getsize(filtered_expression_file) > 0):
                if expression_matrix_path is not None and os.path.isfile(expression_matrix_path):
                    #print(f"Processing expression matrix: {expression_matrix_path}")

                    filtered_lines, expressed, potential_pseudogenes = Candidate_Expression(
                        expression_matrix_path,
                        final_members,
                        args.min_avg_tpm,
                        args.min_single_tpm
                    )

                    # Write filtered expression matrix
                    with open(filtered_expression_file, 'w') as out:
                        for line in filtered_lines:
                            out.write(line + '\n')

                    print(f"{len(expressed)} candidates classified as expressed.")
                    #print(f"{len(potential_pseudogenes)} candidates classified as potential pseudogenes.")
                else:
                    print("No expression matrix provided or file not found - skipping expression filtering step")
            #else:
                #print(f"Filtered expression matrix already exists: {filtered_expression_file}")
            metadata_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix")
            metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")
            metadata_files_exist_and_not_empty = all([
                    os.path.isfile(metadata_expression_file) and os.path.getsize(metadata_expression_file) > 0,
                    os.path.isfile(metadata_only_expression_file) and os.path.getsize(metadata_only_expression_file) > 0
                ])
            if metadata_path is not None and os.path.isfile(metadata_path) and not metadata_files_exist_and_not_empty:
                #print(f"Using metadata from: {metadata_path}")
                use_metadata(filtered_expression_file, metadata_path, metadata_expression_file, metadata_only_expression_file)
            # else:
            #     print("No metadata file provided or file not found - skipping metadata integration")

            highest_expression_file = write_highest_expression_info(
                filtered_expression_file,
                metadata_expression_file=filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix"),
                metadata_only_expression_file=filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")
            )
            print("Step 7 completed.")

        # --- 08 find in species-specific paralogs (in-paralogs) --- #
        print("Step 8: Find in-species paralogs.")
        if args.collapse.upper() in ['YES', 'Y']:
            repr_clean_file = os.path.join(fasta_folder, f"{args.name}08_representative_paralogs_sequences.fasta")
            paralog_group_file = os.path.join(txt_folder, f"{args.name}08_representative_paralogs_summary.txt")
            
            paralog_files_exist_and_not_empty = all([
                os.path.isfile(repr_clean_file) and os.path.getsize(repr_clean_file) > 0,
                os.path.isfile(paralog_group_file) and os.path.getsize(paralog_group_file) > 0
            ])

            if not paralog_files_exist_and_not_empty:
                if expression_matrix_path is not None and metadata_path is not None:
                    metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")
                    paralog_groups = establish_paralog_groups_with_expression_filter(final_tree_path, final_members.keys(), args.paralogdist, 
                                                                                    metadata_only_expression_file, args.min_paralog_tpm)
                else:
                    paralog_groups = establish_paralog_groups(final_tree_path, final_members.keys(), args.paralogdist)
                
                rep_per_group = get_represenative_paralog_per_group(paralog_groups, final_members, repr_clean_file)

                label_map_paralogs = {
                    "group1_baits": "Baits",
                    "group2_representatives": "Representative",
                    "group3_collapsed_paralogs": "CollapsedParalog",
                }
                palette_paralogs = {
                    "group1_baits": "#1f77b4",
                    "group2_representatives": "#2ca02c",
                    "group3_collapsed_paralogs": "#ff7f0e",
                }

                itol_groups_paralogs = build_paralog_groups_for_itol(
                    baits_fasta_path=baits_path,
                    final_members=final_members,
                    paralog_groups=paralog_groups,
                    representatives=rep_per_group,
                )

                write_itol_annotation(
                    groups=itol_groups_paralogs,
                    output_folder=tree_folder,
                    dataset_label=f"{args.name}_Paralog_Representatives",
                    filename_suffix="Paralog_Representatives",
                    label_map=label_map_paralogs,
                    palette=palette_paralogs,
                )  # [1][2]

                with open(paralog_group_file, "w") as out:
                    out.write("RepresentativeSeqID\tMembersOfParalogGroup\n")
                    for gene in rep_per_group.keys():
                        for group in paralog_groups:
                            if gene in group:
                                out.write(gene + "\t" + ";".join(group) + "\n")

            if not (os.path.isdir(tree_folder) and any(f.startswith("08") for f in os.listdir(tree_folder))):
                repr_tree_output_folder = os.path.join(supplement_folder, f"{args.name}08_repr_tree/")
                repr_fin_aln_candidate_file = "08_repr_fin_alignment_candidates.fasta"
                repr_fin_aln_input_file = "08_repr_fin_alignment_input.fasta"
                repr_fin_aln_file = "08_repr_fin_alignment_input.fasta.aln"
                repr_fin_cln_aln_file = "08_repr_fin_alignment_input.fasta.aln.cln"

                repr_tree_file = parallel_tree_constructor(
                    len(final_members),
                    repr_tree_output_folder,
                    repr_fin_aln_candidate_file,
                    repr_fin_aln_input_file,
                    repr_fin_aln_file,
                    repr_fin_cln_aln_file,
                    baits_path,
                    repr_clean_file,
                    "08_representatives_baits_",
                    "",
                    args.mode_aln,
                    args.mode_tree,
                    args.mafft,
                    args.muscle,
                    args.raxml,
                    args.fasttree,
                    args.iqtree,
                    args.cpu_max,
                    "N"
                )

                if repr_tree_file:
                    final_repr_tree_path = replace_spaces_in_tree(repr_tree_file[0])
                    shutil.copy(final_repr_tree_path, tree_folder)
                else:
                    print("Warning: No tree file was created during representative paralog processing.")

            tree_counter += 1
        
        print("Step 8 completed.")

        # --- 09 create summary file --- #
        print("Step 9: Create summary file.")
        summary_file = os.path.join(txt_folder, f"{args.name}09_summary.txt")

        if not (os.path.isfile(summary_file) and os.path.getsize(summary_file) > 0):
            #print("Generating summary file...")

            summary_data = collect_summary_data(
                ref_mapping_file=filtered_ref_mapping_file,
                baits_info=baits_info,
                best_domain_match_file=best_domain_match_file,
                motif_check_file_summary=static_summary_file,
                highest_expression_file=highest_expression_file,
                paralog_group_file=paralog_group_file
            )

            header = [
                "ID", "Bait_Reference_Gene", "", "",
                "Domain_Check",
                "Motifs_Check", "", "", "", "", "",
                "Expression", "", "", "",
                "Paralogs"
            ]

            subheader = [
                "ID", "Ref_Gene", "Ref_Family", "Ref_Subfamily",
                "Best_Domain",
                "Proline-enriched_region", "C-terminal_helix", "I-helix", "EXXR", "PERF", "Heme-binding_region",
                "Highest_Run", "TPM", "Highest_Run_with_Metadata", "TPM_with_Metadata",
                "repr_paralog"
            ]

            with open(summary_file, "w") as out:
                out.write("\t".join(header) + "\n")
                out.write("\t".join(subheader) + "\n")
                for seq_id in sorted(summary_data.keys()):
                    out.write("\t".join(summary_data[seq_id]) + "\n")

            #print(f"Summary file written to: {summary_file}")
        else:
            print(f"Summary file already exists: {summary_file}")

        print("Step 9 completed.")

        # --- 10 Export --- #
        print("Step 10: Export.")
        # --- Tree plotting --- #
        ete3_available = importlib.util.find_spec("ete3") is not None
        if ete3_available and (not os.path.exists(tree_folder) or not any(f.endswith((".png", ".pdf", ".svg")) for f in os.listdir(tree_folder)) or len(os.listdir(tree_folder)) <  (tree_counter * 3)):
            #print("Searching for .tre files and plotting trees...")
            plot_all_trees(tree_folder, tree_folder, baits_path, filtered_members_file, format="png", layout="Polar", cand_color="green", bait_color="red")
        elif not ete3_available:
            print("ete3 is not installed. If you wish to convert the treefiles, please install ete3 using 'pip install ete3'")
        else:
            print(f"Tree visualizations already exist in: {tree_folder}")

        # --- Create Heatmaps --- #
        if os.path.isdir(heatmap_folder):
            seaborn_available = importlib.util.find_spec("seaborn") is not None

            # Check if seaborn is available and if heatmaps already exist
            heatmaps_exist = any(f.endswith((".png", ".pdf", ".svg")) for f in os.listdir(heatmap_folder))

            if seaborn_available and (not heatmaps_exist or len(os.listdir(heatmap_folder)) < 6):
                #print("Create Heatmaps from expression data...")
                generate_heatmaps_from_expression_data(txt_folder, heatmap_folder, format="png")
            elif not seaborn_available:
                print("seaborn is not available – skipping heatmap creation.")
            else:
                print(f"Heatmaps already exist in: {heatmap_folder}")

        # --- documentation of execution time --- #
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        seconds = duration.total_seconds()
        with open(doc_file, "a") as out:
            out.write("\n" + "Execution time: " + f"{int(seconds // 3600)}h:{int((seconds % 3600) // 60)}min:{int(seconds % 60)}sec")               

        # --- Export + Overview --- #
        html_link_file = os.path.join(html_folder, "CYP_output_links.html")

        if not (os.path.isfile(html_link_file) and os.path.getsize(html_link_file) > 0):
            #print("Creating all HTML files and overview index...")
            create_output_html_index(txt_folder, html_folder, html_link_file,
                                     tree_folder, fasta_folder, heatmap_folder)
            #print(f"Output HTML overview written to: {html_link_file}")
        else:
            print(f"HTML overview already exists: {html_link_file}")     

        print("Step 10 completed.")

        sys.stdout.write("\n" + "Successfully finished execution. " + "\n")
        print("\n")

# endregion

if __name__ == "__main__":
    main()