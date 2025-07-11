# region Import
import argparse
import csv
import dendropy
import datetime
import importlib.util
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
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
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from operator import itemgetter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

# endregion

# region Constants
FASTA_EXTENSIONS = ('.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')
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
class Input:
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
        "--subject",
        type=str,
        nargs="+",  # Accepts one or more paths
        default=None,
        help="Path to subject (FASTA file(s) or folder(s)). Required if --data is not provided."
    )

    # Optional arguments
    parser.add_argument(
        "--outgroup",
        type=str,
        default=None,
        help="Path to outgroup FASTA file (optional)"
    )
    

    # Optional files
    parser.add_argument(
        "--protein_motifs",
        type=str,
        default=None,
        help="Path to protein motifs file (optional)"
    )
    parser.add_argument(
        "--hmm_motifs",
        type=str,
        default=None,
        help="Path to hmm motifs file (optional)"
    )
    parser.add_argument(
        "--hmm_domains",
        type=str,
        default=None,
        help="Path to hmm domains file (optional)"
    )
    parser.add_argument(
        "--baits_info",
        type=str,
        default=None,
        help="Path to baits info file (optional)"
    )

    # Motif Check
    parser.add_argument(
        "--static_motifs",
        type=str,
        default="n",
        choices=['y', 'n', 'yes', 'no'],
        help="Use static motif search instead of hmm (y/n, default: n)"
    )
    parser.add_argument(
        "--motif_cEvalue",
        type=float,
        default=0.001,
        help="c-Evalue for hmm motif integration (default: 0.001)"
    )

    # Domain Check
    parser.add_argument(
        "--domain_Score",
        type=float,
        default=100,
        help="c-Evalue for hmm motif integration (default: 100)"
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

    # Folders
    parser.add_argument(
        "--output_folder",
        type=str,
        default=None,
        help="Name of the output folder (optional)"
    )
    parser.add_argument(
        "--processed_input_folder",
        type=str,
        default=None,
        help="Name of the processed input folder (optional)"
    )


    # parser.add_argument(
    #     "--info",
    #     type=str,
    #     default=None,
    #     help="Path to info CSV file (optional)"
    # )
    parser.add_argument(
        "--trim_names",
        type=str,
        default="n",
        choices=['y', 'n', 'yes', 'no'],
        help="Trim sequence names at first space or tab (y/n, default: n)"
    )
    parser.add_argument(
        "--name",
        type=str,
        default="",
        help="STRING_USED_AS_PREFIX_IN_FILENAMES"
    )
    parser.add_argument(
        "--hmm",
        type=str,
        default=None,
        help="Path to bait HMM file"
    )
    parser.add_argument(
        "--use_hmm",
        type=str,
        default="n",
        choices=['y', 'n', 'yes', 'no'],
        help="Reduce BLAST candidates using HMMER (y/n, default: n)"
    )
    parser.add_argument(
        "--cpub",
        type=int,
        default=4,
        help="CPUs_TO_USE_FOR_BLASTp (default = 4)"
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
        default=80,
        help="BLASTP_MIN_LENGTH_CUTOFF (default = 80)"
    )
    parser.add_argument(
        "--bitcutp",
        type=int,
        default=60,
        help="BLASTP_BITSCORE_CUTOFF (default = 60)"
    )
    parser.add_argument(
        "--numprocesscandidates",
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
        choices=['fasttree', 'raxml'],
        help="Tool used for tree construction (fasttree/raxml, default: fasttree)"
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
        "--parallel",
        type=str,
        default="y",
        choices=['y', 'n', 'yes', 'no'],
        help="Run script in parallel mode (y/n, default: y)"
    )
    parser.add_argument(
        "--cpumax",
        type=int,
        default=4,
        help="MAX_CPUs_FOR_IN_OUT_CLASSIFICATION"
    )
    parser.add_argument(
        "--cpur",
        type=int,
        default=4,
        help="CPUs_TO_USE_FOR_RAxML"
    )
    parser.add_argument(
        "--numneighbours",
        type=int,
        default=10,
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
        default="y",
        choices=['y', 'n', 'yes', 'no'],
        help="DOMAIN_FILTER_FOR_CLASSIFICATION (y/n, default: y)"
    )

    # Orthologs
    parser.add_argument(
        "--ortholog_threshold",
        type=float,
        default=1.0,
        help="Theshold for patristic distance considering orthologs (default: 1.0)"
    )
    parser.add_argument(
        "--neighbour_threshold",
        type=float,
        default=2.0,
        help="Theshold for patristic considering further neighbours (default: 2.0)"
    )
    parser.add_argument(
        "--reroot",
        type=str,
        default="XP-020521354.1",
        help="Sequence ID to reroot final tre ('_' need to be replaced by '-')"
    )

    # create seperate trees for different bait qualities
    # Wahrscheinlich zu löschen
    parser.add_argument(
        "--transcript",
        type=str,
        default="y",
        choices=['y', 'n', 'yes', 'no'],
        help="Create separate tree with transcript level evidence baits as lowest quality? (y/n, default: y)"
    )
    parser.add_argument(
        "--protein",
        type=str,
        default="y",
        choices=['y', 'n', 'yes', 'no'],
        help="Create separate tree with protein level evidence baits as lowest quality? (y/n, default: y)"
    )
    parser.add_argument(
        "--literature",
        type=str,
        default="y",
        choices=['y', 'n', 'yes', 'no'],
        help="Create separate tree with only literature based baits? (y/n, default: y)"
    )

    @staticmethod
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

    # Check essential tool groups
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

def collect_files(
    data_folder: Optional[str] = None,
    subject_paths: Optional[List[str]] = None,
    bait_path: Optional[str] = None,
    protein_motifs_path: Optional[str] = None,
    hmm_motifs_path: Optional[str] = None,
    hmm_domains_path: Optional[str] = None,
    baits_info_path: Optional[str] = None,
    expression_path: Optional[str] = None,
    metadata_path: Optional[str] = None,
    output_dir: str = "Processed_Input",
    trim_names: str = 'n'
) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """
    Collect all input files from data folder and/or explicit paths.
    Returns tuple of (file_paths, processed_paths) dictionaries.
    """
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(exist_ok=True)
    
    file_paths = {
        'subject_files': [],
        'bait_file': None,
        'baits_info_file': None,
        'hmm_domains_file': None,
        'protein_motifs_file': None,
        'hmm_motifs_file': None,
        'expression_file': None,
        'metadata_file': None
    }
    
    # Collect from data folder if provided
    if data_folder:
        data_path = Path(data_folder)
        if not data_path.is_dir():
            raise ValueError(f"Data folder does not exist: {data_folder}")
        
        # Subject files
        subject_file = data_path / "subject.fasta"
        if subject_file.exists():
            file_paths['subject_files'].append(str(subject_file))

        subjects_dir = data_path / "subjects"
        if subjects_dir.exists():
            file_paths['subject_files'].extend(collect_fasta_files([str(subjects_dir)]))
        
        # Bait and info files
        bait_file = data_path / "baits.fasta"
        if bait_file.exists():
            file_paths['bait_file'] = str(bait_file)

        baits_info_file = data_path / "baits_info.csv"
        if baits_info_file.exists():
            file_paths['baits_info_file'] = str(baits_info_file)

        # Optional files
        file_paths['protein_motifs_file'] = find_case_insensitive_file(data_path, "protein_motifs.txt")
        file_paths['hmm_motifs_file'] = find_case_insensitive_file(data_path, "hmm_motifs.hmm")
        file_paths['hmm_domains_file'] = find_case_insensitive_file(data_path, "hmm_domains.hmm")
        file_paths['expression_file'] = find_case_insensitive_file(data_path, "expression.txt")
        file_paths['metadata_file'] = find_case_insensitive_file(data_path, "metadata.csv")

    # Collect from explicit paths (overrides data folder files)
    if subject_paths:
        valid_subject_paths = [p for p in subject_paths if p is not None]
        if valid_subject_paths:
            file_paths['subject_files'].extend(collect_fasta_files(valid_subject_paths))
    
    if bait_path and os.path.isfile(bait_path):
        file_paths['bait_file'] = bait_path
    
    # Process optional files
    file_paths['protein_motifs_file'] = protein_motifs_path if (protein_motifs_path and os.path.isfile(protein_motifs_path)) else file_paths['protein_motifs_file']
    file_paths['hmm_motifs_file'] = hmm_motifs_path if (hmm_motifs_path and os.path.isfile(hmm_motifs_path)) else file_paths['hmm_motifs_file']
    file_paths['hmm_domains_file'] = hmm_domains_path if (hmm_domains_path and os.path.isfile(hmm_domains_path)) else file_paths['hmm_domains_file']
    file_paths['baits_info_file'] = baits_info_path if (baits_info_path and os.path.isfile(baits_info_path)) else file_paths['baits_info_file']
    file_paths['expression_file'] = expression_path if (expression_path and os.path.isfile(expression_path)) else file_paths['expression_file']
    file_paths['metadata_file'] = metadata_path if (metadata_path and os.path.isfile(metadata_path)) else file_paths['metadata_file']

    # Validate required files
    if not file_paths['subject_files']:
        raise ValueError("No valid subject files found")
    if not file_paths['bait_file']:
        raise ValueError("No valid bait file found")
    if not file_paths['baits_info_file']:
        raise ValueError("No valid baits info file found")

    # Process all files
    processed_paths = process_all_files(file_paths, output_dir, trim_names)
    
    return file_paths, processed_paths

def process_all_files(file_paths: Dict[str, List[str]], output_dir: str = "Processed_Input", trim_names: str = 'n') -> Dict[str, str]:
    """
    Process all collected files through the pipeline.
    Returns dictionary of processed file paths.
    """
    output_dir_path = Path(output_dir)
    processed_paths = {}

    # Process bait file and create group-specific files
    bait_results = process_baits_file(
        file_paths['bait_file'],
        file_paths['baits_info_file'],
        output_dir,
        trim_names
    )
    processed_paths.update({
        'baits_path': bait_results[0],
        'baits_no_outgroup_path': bait_results[1],
        'transcript_baits_path': bait_results[2],
        'protein_baits_path': bait_results[3],
        'literature_baits_path': bait_results[4],
        'outgroup_path': bait_results[5]
    })

    # Process subject files
    subject_files_paths = []
    for subject_file in file_paths['subject_files']:
        subject_files_paths.append(process_single_file(subject_file, Path(subject_file).stem, output_dir, trim_names))
    processed_paths['subject_files_paths'] = subject_files_paths

    # Process other files
    if file_paths['protein_motifs_file']:
        processed_paths['protein_motifs_path'] = copy_file(
            file_paths['protein_motifs_file'],
            "protein_motifs",
            output_dir
        )

    if file_paths['hmm_motifs_file']:
        processed_paths['hmm_motifs_path'] = copy_file(
            file_paths['hmm_motifs_file'],
            "hmm_motifs",
            output_dir
        )

    if file_paths['hmm_domains_file']:
        processed_paths['hmm_domains_path'] = copy_file(
            file_paths['hmm_domains_file'],
            "hmm_domains",
            output_dir
        )

    if file_paths['baits_info_file']:
        processed_paths['baits_info_path'] = copy_file(
            file_paths['baits_info_file'],
            "baits_info",
            output_dir
        )

    if file_paths['expression_file']:
        processed_paths['expression_path'] = process_expression_file(
            file_paths['expression_file'],
            "expression",
            output_dir,
            trim_names
        )

    if file_paths['metadata_file']:
        processed_paths['metadata_path'] = copy_file(
            file_paths['metadata_file'],
            "metadata",
            output_dir
        )

    return processed_paths

def process_baits_file(
    bait_file: str,
    baits_info_file: str,
    output_dir: str = "Processed_Input",
    trim_names: str = 'n'
) -> Tuple[str, str, str, str, str, str]:
    """
    Process bait file and create group-specific files based on baits_info.csv.
    Returns tuple of paths: (all_baits, no_outgroup, transcript, protein, literature, outgroup)
    """
    output_dir_path = Path(output_dir)
    
    # Process main bait file
    all_baits_path = process_single_file(bait_file, "baits_with_outgroup", output_dir, trim_names)
    sequences = read_fasta(all_baits_path)

    # Process baits info file    
    baits_info, info_headers = read_baits_info(baits_info_file)
    bait_groups = {
        seq_id: subdict.get("Type", "-")
        for seq_id, subdict in baits_info.items()
        }

    # Initialize output files
    output_files = {
        'no_outgroup': output_dir_path / "clean_baits.fasta",
        'transcript': output_dir_path / "transcript_baits.fasta",
        'protein': output_dir_path / "protein_baits.fasta",
        'literature': output_dir_path / "literature_baits.fasta",
        'outgroup': output_dir_path / "outgroup.fasta"
    }

    # Check if files already exist
    all_exist = all(f.exists() for f in output_files.values())
    if all_exist:
        print("Group-specific bait files already exist, skipping creation.")
        return (
            all_baits_path,
            str(output_files['no_outgroup']),
            str(output_files['transcript']),
            str(output_files['protein']),
            str(output_files['literature']),
            str(output_files['outgroup'])
        )

    # Sort sequences into groups
    grouped_sequences = {
        'no_outgroup': {},
        'transcript': {},
        'protein': {},
        'literature': {},
        'outgroup': {}
    }

    for seq_id, seq in sequences.items():
        group = bait_groups.get(seq_id, "")
        group = group.lower()
        
        if "outgroup" in group:
            grouped_sequences['outgroup'][seq_id] = seq
        else:
            grouped_sequences['no_outgroup'][seq_id] = seq
            if "transcript" in group:
                grouped_sequences['transcript'][seq_id] = seq
            if "protein" in group:
                grouped_sequences['protein'][seq_id] = seq
            if "literature" in group:
                grouped_sequences['literature'][seq_id] = seq

    # Write output files only if they don't exist or are empty
    for group_name, seq_dict in grouped_sequences.items():
        if seq_dict:  # Only write if we have sequences
            output_path = output_files[group_name]
            if not output_path.exists() or output_path.stat().st_size == 0:
                write_fasta(seq_dict, output_path)

    return (
        all_baits_path,
        str(output_files['no_outgroup']),
        str(output_files['transcript']),
        str(output_files['protein']),
        str(output_files['literature']),
        str(output_files['outgroup'])
    )

def read_baits_info(baits_info_file: str) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    Read baits info CSV file and return:
    - a dictionary of sequence ID to subdictionary (column -> value)
    - the original column headers in order
    """
    baits_info = {}
    headers = []

    try:
        with open(baits_info_file, 'r', newline='') as f:
            reader = csv.reader(f)
            headers = next(reader)  # Read header line
            for row in reader:
                if len(row) < 1:
                    continue
                seq_id = row[0].strip()
                subdict = {
                    headers[i].strip(): row[i].strip()
                    for i in range(1, min(len(headers), len(row)))
                }
                baits_info[seq_id] = subdict
    except Exception as e:
        raise ValueError(f"Invalid baits info file {baits_info_file}: {str(e)}")

    return baits_info, headers


def process_single_file(
    input_file: str,
    output_name: str,
    output_dir: str = "Processed_Input",
    trim_names: str = 'n'
) -> str:
    """Process a single FASTA file through the complete pipeline."""
    output_dir_path = Path(output_dir)
    output_fasta = output_dir_path / f"{output_name}.fasta"
    mapping_file = output_dir_path / f"{output_name}_mapping.txt"

    # Skip processing if output files already exist
    if output_fasta.exists() and mapping_file.exists():
        print(f"{input_file} already processed, skipping.")
        return str(output_fasta)

    # Load and process sequences
    sequences = load_and_process_sequences(input_file, trim_names)
    
    # Clean sequence IDs and create mapping
    cleaned_sequences = {}
    sequence_mapping = {}
    
    for original_id, sequence in sequences.items():
        cleaned_id = clean_sequence_id(original_id, trim_names)
        cleaned_sequences[cleaned_id] = sequence
        sequence_mapping[original_id] = cleaned_id

    # Write outputs
    write_fasta(cleaned_sequences, output_fasta)
    write_mapping_table(sequence_mapping, mapping_file)
    
    print(f"Successfully processed {input_file} to: {output_fasta}")
    return str(output_fasta)

def load_and_process_sequences(file_path: str, trim_names: str = 'n') -> Dict[str, str]:
    """Load FASTA sequences and translate if nucleotide."""
    sequences = read_fasta(file_path)
    
    # Determine sequence type and translate if needed
    if sequences and is_nucleotide_sequence(next(iter(sequences.values()))):
        sequences = translate_sequences(sequences)
    
    # Clean sequence IDs if requested
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
    nucleotide_chars = set("ACGTN")
    return all(c in nucleotide_chars for c in clean_seq)

def translate_sequences(sequences: Dict[str, str]) -> Dict[str, str]:
    """Translate nucleotide sequences to peptide sequences."""
    translated = {}
    
    for seq_id, seq in sequences.items():
        clean_seq = re.sub(r'[\- ]', '', seq)
        peptide = []
        
        for i in range(len(clean_seq) // 3):
            codon = clean_seq[i*3:i*3+3]
            peptide.append(GENETIC_CODE.get(codon, "*"))
        
        translated[seq_id] = "".join(peptide)
    
    return translated

def clean_sequence_id(seq_id: str, trim_names: str = 'n') -> str:
    """Clean sequence ID by removing forbidden characters and optionally trimming."""
    # Trim at first space or tab if requested
    if trim_names.upper() in ['Y', 'YES']:
        seq_id = seq_id.split()[0] if seq_id.split() else seq_id
    
    # Replace forbidden characters
    for char in FORBIDDEN_CHARS:
        seq_id = seq_id.replace(char, "-")
    
    # Encode to ASCII to remove non-ASCII characters
    seq_id = seq_id.encode("ascii", "ignore").decode()
    
    return seq_id

def write_fasta(sequences: Dict[str, str], output_path: Path) -> None:
    """Write sequences to FASTA file."""
    # Skip if file already exists and is not empty
    if output_path.exists() and output_path.stat().st_size > 0:
        return
    
    with open(output_path, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

def write_mapping_table(mapping: Dict[str, str], output_path: Path) -> None:
    """Write sequence ID mapping table."""
    # Skip if file already exists and is not empty
    if output_path.exists() and output_path.stat().st_size > 0:
        return
    
    with open(output_path, 'w') as f:
        f.write("InitialID\tCleanID\n")
        for original_id, cleaned_id in mapping.items():
            f.write(f"{original_id}\t{cleaned_id}\n")

def process_expression_file(
    input_file: str,
    output_name: str,
    output_dir: str = "Processed_Input",
    trim_names: str = 'n'
) -> str:
    """Process expression file, cleaning IDs in the first column."""
    output_path = Path(output_dir) / f"{output_name}.txt"
    
    # Skip if file already exists and is not empty
    if output_path.exists() and output_path.stat().st_size > 0:
        return str(output_path)
    
    try:
        with open(input_file, 'r') as infile, open(output_path, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)
            
            for line in infile:
                parts = line.strip().split('\t')
                if not parts:
                    continue
                
                # Clean the ID in first column
                parts[0] = clean_sequence_id(parts[0], trim_names)
                outfile.write('\t'.join(parts) + '\n')
        
        return str(output_path)
    except Exception as e:
        raise ValueError(f"Failed to process expression file {input_file}: {str(e)}")

def copy_file(
    input_file: str,
    output_name: str,
    output_dir: str = "Processed_Input"
) -> str:
    """Copy file to output directory with new name."""
    output_dir_path = Path(output_dir)
    if "hmm" in output_name:
        output_path = output_dir_path / f"{output_name}.hmm"
    else:
        output_path = output_dir_path / f"{output_name}.txt"
    
    # Skip if file already exists and is not empty
    if output_path.exists() and output_path.stat().st_size > 0:
        return str(output_path)
    
    try:
        subprocess.run(["cp", input_file, output_path], check=True)
        return str(output_path)
    except Exception as e:
        raise ValueError(f"Failed to copy {input_file} to {output_path}: {str(e)}")

def collect_fasta_files(paths: List[str]) -> List[str]:
    """Collect FASTA files from given paths (files or directories)."""
    files = []
    for path_str in paths:
        path = Path(path_str)
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS:
            files.append(str(path))
        elif path.is_dir():
            for f in path.glob("*.fa*"):
                if f.suffix.lower() in FASTA_EXTENSIONS:
                    files.append(str(f))
        else:
            print(f"Warning: Path not found or invalid: {path}")
    return files

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
    """Load subject name mapping table from a tab-delimited file.
    
    Args:
        mapping_table_file: Path to the mapping table file
        
    Returns:
        Dictionary mapping clean IDs to original IDs
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
#endregion

# region Search
def run_blast_search(query_file, subject_file, output_file, blast_db_folder, cpub=1, makeblastdb="makeblastdb", blastp="blastp"):
    """Run BLAST search and return results file path"""
    blast_db = os.path.join(blast_db_folder, "blastdb")
    
    # Create BLAST database
    p = subprocess.Popen(
        args=f"{makeblastdb} -in {subject_file} -out {blast_db} -dbtype prot",
        shell=True
    )
    p.communicate()
    
    # Run BLAST search
    p = subprocess.Popen(
        args=f"{blastp} -query {query_file} -db {blast_db} -out {output_file} "
                f"-outfmt 6 -evalue 0.001 -num_threads {cpub}",
        shell=True
    )
    p.communicate()
    
    return output_file
    
def run_hmmsearch(hmm_file: str, subject_file: str, output_file: str, hmmsearch: str = "hmmsearch", domtblout: bool = False) -> str:
    """Run HMMER search and return results file path
    
    Args:
        hmm_file: Path to HMM file
        subject_file: Path to subject sequences file
        output_file: Path for output results
        hmmsearch: Path to hmmsearch executable
        domtblout: Whether to use domtblout format
        
    Returns:
        Path to results file
    """
    waste_file = os.path.join(os.path.dirname(output_file), "hmmsearch_waste.txt")
    
    # Build command based on output format
    if domtblout:
        cmd = f"{hmmsearch} --domtblout {output_file} {hmm_file} {subject_file} > {waste_file}"
    else:
        cmd = f"{hmmsearch} --tblout {output_file} {hmm_file} {subject_file} > {waste_file}"
    
    # Execute command
    p = subprocess.Popen(args=cmd, shell=True)
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
def create_alignments(num_par: int, tree_output_folder: str, name: str, 
                        number: int, aln_input_file: str, aln_file: str, 
                        mode_aln: str, mafft: str, muscle: str, cpu_par: int, 
                        cpu_max: int) -> None:
    """
    Creates alignments for all candidate chunks in parallel.
    
    Args:
        num_par: Number of parallel chunks
        tree_output_folder: Output directory path
        Xname: Base name for files
        Xnumber: Number identifier
        aln_input_file: Input alignment file suffix
        aln_file: Output alignment file suffix
        mode_aln: Alignment mode ("mafft" or "muscle")
        mafft: Path to MAFFT executable
        muscle: Path to MUSCLE executable
        cpu_par: CPUs per process
        cpu_max: Maximum CPUs available
    """
    alignment_processes = []
    cpu_use = 0
    
    for i in range(num_par):
        aln_input = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_input_file}")
        aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_file}")
        
        # Wait for CPU resources if needed
        while cpu_use + cpu_par > cpu_max:
            time.sleep(2)  # Shorter sleep interval for better responsiveness
            # Clean up finished processes and recalculate CPU usage
            active_processes = []
            for p in alignment_processes:
                if p.poll() is None:  # Process still running
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
                print(f"Started alignment process {i+1}")
                
            except Exception as e:
                print(f"Error starting alignment {i+1}: {str(e)}")
    
    # Wait for all alignment processes to complete
    for p in alignment_processes:
        returncode = p.wait()
        if returncode != 0:
            print(f"Warning: Alignment process failed with return code {returncode}")
    
    print("All alignment processes completed")

def alignment_trimming(aln_file: str, cln_aln_file: str, occupancy: float = 0.1) -> None:
    """
    Trims alignment by removing columns with low occupancy.
    
    Args:
        aln_file: Input alignment file path
        cln_aln_file: Output cleaned alignment file path
        occupancy: Minimum occupancy threshold (default: 0.1)
    """
    alignment = load_alignment(aln_file, {})
    
    # --- if there is an alignment (expected case) 
    if len(list(alignment.keys())) > 0:
        # --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
        valid_index = []
        for idx, aa in enumerate(list(alignment.values())[0]):
            counter = 0
            for key in list(alignment.keys()):
                if alignment[key][idx] != "-":
                    counter += 1
            if counter / float(len(list(alignment.keys()))) > occupancy:
                valid_index.append(idx)
            
        # --- generate new sequences --- #
        with open(cln_aln_file, "w") as out:
            for key in list(alignment.keys()):
                seq = alignment[key]
                new_seq = []
                for idx in valid_index:
                    new_seq.append(seq[idx])
                out.write(">" + key + '\n' + "".join(new_seq) + '\n')
    # --- just in case the alignment file is empty (is this possible?) ---#
    else:
        with open(cln_aln_file, "w") as out:
            out.write("")

def load_alignment(aln_file: str, tmp_mapping: Dict[str, str]) -> Dict[str, str]:
    """
    Loads alignment from FASTA file.
    
    Args:
        aln_file: Path to alignment file
        tmp_mapping: Dictionary for header mapping
        
    Returns:
        Dictionary with sequence headers as keys and sequences as values
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

# def parallel_tree_constructor(num_process_candidates: int, tree_output_folder: str, 
#                             aln_candidate_file: str, aln_input_file: str, 
#                             aln_file: str, cln_aln_file: str, bait_file: str, 
#                             candidate_file: str, name: str, number: int, 
#                             mode_aln: str, mode_tree: str, mafft: str, muscle: str, 
#                             raxml: str, fasttree: str, cpu_max: int, cpur: int, parallel_mode: str,
#                             cand_color: str = "green", bait_color: str = "red") -> List[str]:

#     """
#     Handles parallel construction of alignments and phylogenetic trees.
    
#     Returns a list of paths to the generated tree files.
#     """
#     # Input validation
#     if not os.path.isfile(bait_file):
#         raise FileNotFoundError(f"Bait sequence file not found: {bait_file}")
#     if not os.path.isfile(candidate_file):
#         raise FileNotFoundError(f"Candidate file not found: {candidate_file}")
    
#     try:
#         candidates = read_fasta(candidate_file)
#         baits = read_fasta(bait_file)
#         cand_ids = candidates.keys()
#         bait_ids = baits.keys()
#         if not candidates:
#             raise ValueError("No candidates found in input file")
        
#         # Calculate parallelization parameters
#         candidates_number = min(num_process_candidates, len(candidates))
#         num_par = math.ceil(len(candidates) / candidates_number)
#         cpu_par = math.floor(cpu_max / num_par) if parallel_mode.upper() in ["YES", "Y"] else cpu_max
        
#         # Adjust CPU allocation
#         if parallel_mode.upper() in ["YES", "Y"] and cpu_par < 8:
#             cpu_par = max(cpur, 1)
        
#         if parallel_mode.upper() in ["NO", "N"]:
#             num_par = 1
#             candidates_number = len(candidates)
        
#         # Create output directory
#         Path(tree_output_folder).mkdir(parents=True, exist_ok=True)
        
#         print(f"Processing {len(candidates)} candidates in {num_par} parallel chunks")
        
#         # Step 1: Create candidate files and combine with bait sequences
#         candidate_keys = list(candidates.keys())
#         # Shuffle candidate files
#         random.shuffle(candidate_keys)

#         for i in range(num_par):
#             start = i * candidates_number
#             end = min((i + 1) * candidates_number, len(candidates))
            
#             cand_file = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_candidate_file}")
#             aln_input = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_input_file}")
            
#             if not os.path.isfile(aln_input):
#                 # Create candidate file
#                 with open(cand_file, "w") as out:
#                     for c in candidate_keys[start:end]:
#                         out.write(f'>{c}\n{candidates[c]}\n')
                
#                 # Combine bait and candidate files safely
#                 with open(aln_input, 'w') as outfile:
#                     with open(bait_file, 'r') as infile1:
#                         outfile.write(infile1.read())
#                     with open(cand_file, 'r') as infile2:
#                         outfile.write(infile2.read())
        
#         # Step 2: Create alignments using the Alignment class
#         create_alignments(
#             num_par=num_par,
#             tree_output_folder=tree_output_folder,
#             name=name,
#             number=number,
#             aln_input_file=aln_input_file,
#             aln_file=aln_file,
#             mode_aln=mode_aln,
#             mafft=mafft,
#             muscle=muscle,
#             cpu_par=cpu_par,
#             cpu_max=cpu_max
#         )
        
#         # Step 3: Create phylogenetic trees with improved process management
#         tree_processes = []
#         trees = []
#         cpu_use = 0
        
#         for i in range(num_par):
#             aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_file}")
#             cln_aln = os.path.join(tree_output_folder, f"{name}{number}{i}_{cln_aln_file}")
            
#             # Wait for CPU resources if needed
#             while cpu_use + cpu_par > cpu_max:
#                 time.sleep(2)
#                 # Clean up finished processes and recalculate CPU usage
#                 active_processes = []
#                 for p in tree_processes:
#                     if p.poll() is None:  # Process still running
#                         active_processes.append(p)
#                     else:
#                         # Check for errors in completed processes
#                         if p.returncode != 0:
#                             print(f"Warning: Tree construction process failed with return code {p.returncode}")
#                 tree_processes = active_processes
#                 cpu_use = len(tree_processes) * cpu_par
            
#             # Perform alignment trimming using the Alignment class
#             if not os.path.isfile(cln_aln):
#                 try:
#                     alignment_trimming(aln, cln_aln, occupancy=0.01)
#                 except Exception as e:
#                     print(f"Alignment trimming failed for chunk {i+1}: {str(e)}")
#                     continue
            
#             # Tree construction
#             try:
#                 if mode_tree == "raxml":
#                     prefix = os.path.join(tree_output_folder, f"{name}{number}{i}_RAxML_tree")
#                     tree_file = f"{prefix}.raxml.bestTree"
#                     trees.append(tree_file)
                    
#                     if not os.path.isfile(tree_file):
#                         cmd = [raxml, "--all", "--threads", str(cpu_par), "--model", "LG+G8+F", 
#                                 "--msa", cln_aln, "--prefix", prefix]
#                         p = subprocess.Popen(cmd, stderr=subprocess.PIPE)
#                         tree_processes.append(p)
#                         cpu_use += cpu_par
#                         print(f"Started RAxML tree construction {i+1}")
#                 else:
#                     tree_file = os.path.join(tree_output_folder, f"{name}{number}{i}_FastTree_tree.tre")
#                     trees.append(tree_file)
                    
#                     if not os.path.isfile(tree_file):
#                         cmd = [fasttree, "-wag", "-nopr"]
#                         with open(cln_aln, 'r') as infile, open(tree_file, 'w') as outfile:
#                             p = subprocess.Popen(cmd, stdin=infile, stdout=outfile, stderr=subprocess.PIPE)
#                         tree_processes.append(p)
#                         cpu_use += cpu_par
#                         print(f"Started FastTree construction {i+1}")
                        
#             except Exception as e:
#                 print(f"Error starting tree construction {i+1}: {str(e)}")
        
#         # Wait for all tree construction processes to complete
#         for p in tree_processes:
#             returncode = p.wait()
#             if returncode != 0:
#                 print(f"Warning: Tree construction process failed with return code {returncode}")
        
#         print(f"Successfully created {len(trees)} phylogenetic trees")
        
#         # Filter out trees that weren't actually created due to errors
#         existing_trees = [tree for tree in trees if os.path.isfile(tree)]
        
#         if "ete3" in sys.modules:
#             for tree_file in existing_trees:
#                 try:
#                     t = Tree(tree_file)
#                     for node in t.traverse():
#                         if node.is_leaf():
#                             nstyle = NodeStyle()
#                             nstyle["size"] = 0  # No circle
#                             node.set_style(nstyle)
                            
#                             # Color assignment via text label
#                             if node.name in bait_ids:
#                                 face = TextFace(node.name, fsize=10, fgcolor=bait_color)
#                             elif node.name in cand_ids:
#                                 face = TextFace(node.name, fsize=10, fgcolor=cand_color)
#                             else:
#                                 face = TextFace(node.name, fsize=10, fgcolor="black")

#                             node.add_face(face, column=0, position="branch-right")
                    
#                     # Tree display style
#                     ts = TreeStyle()
#                     ts.show_leaf_name = False  # We add leaf names manually as faces
#                     ts.show_branch_support = True
#                     ts.scale = 120

#                     img_out = tree_file + ".pdf"
#                     t.render(img_out, tree_style=ts, w=800)
#                     print(f"Colored tree saved to {img_out}")
#                 except Exception as e:
#                     print(f"Failed to render tree {tree_file}: {str(e)}")
#         else:
#             print("ETE3 not installed – skipping tree colorization.")

#         return existing_trees
        
#     except Exception as e:
#         print(f"Error in parallel_tree_constructor: {str(e)}")
#         raise

def parallel_tree_constructor(num_process_candidates: int, tree_output_folder: str, 
                            aln_candidate_file: str, aln_input_file: str, 
                            aln_file: str, cln_aln_file: str, bait_file: str, 
                            candidate_file: str, name: str, number: int, 
                            mode_aln: str, mode_tree: str, mafft: str, muscle: str, 
                            raxml: str, fasttree: str, cpu_max: int, cpur: int, parallel_mode: str,
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
        cand_ids = set(candidates.keys()) # Verwende Sets für schnellere Lookups
        bait_ids = set(baits.keys())     # Verwende Sets für schnellere Lookups
        if not candidates:
            raise ValueError("No candidates found in input file")
        
        # Calculate parallelization parameters
        candidates_number = min(num_process_candidates, len(candidates))
        num_par = math.ceil(len(candidates) / candidates_number)
        cpu_par = math.floor(cpu_max / num_par) if parallel_mode.upper() in ["YES", "Y"] else cpu_max
        
        # Adjust CPU allocation
        if parallel_mode.upper() in ["YES", "Y"] and cpu_par < 8:
            cpu_par = max(cpur, 1)
        
        if parallel_mode.upper() in ["NO", "N"]:
            num_par = 1
            candidates_number = len(candidates)
        
        # Create output directory
        Path(tree_output_folder).mkdir(parents=True, exist_ok=True)
        
        print(f"Processing {len(candidates)} candidates in {num_par} parallel chunks")
        
        # Step 1: Create candidate files and combine with bait sequences
        candidate_keys = list(candidates.keys())
        # Shuffle candidate files
        random.shuffle(candidate_keys)

        for i in range(num_par):
            start = i * candidates_number
            end = min((i + 1) * candidates_number, len(candidates))
            
            cand_file = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_candidate_file}")
            aln_input = os.path.join(tree_output_folder, f"{name}{number}{i}_{aln_input_file}")
            
            # Überprüfen, ob die aln_input Datei bereits existiert und nicht leer ist
            # Dies verhindert unnötiges Überschreiben bei wiederholten Läufen
            if not os.path.isfile(aln_input) or os.path.getsize(aln_input) == 0:
                # Create candidate file
                with open(cand_file, "w") as out:
                    for c in candidate_keys[start:end]:
                        out.write(f'>{c}\n{candidates[c]}\n')
                
                # Combine bait and candidate files safely
                with open(aln_input, 'w') as outfile:
                    with open(bait_file, 'r') as infile1:
                        outfile.write(infile1.read())
                    with open(cand_file, 'r') as infile2:
                        outfile.write(infile2.read())
                print(f"Combined bait and candidate file created: {aln_input}")
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
                        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE) # Capture stdout as well
                        tree_processes.append(p)
                        cpu_use += cpu_par
                        print(f"Started RAxML tree construction {i+1} with PID {p.pid}")
                    else:
                        print(f"Skipping RAxML for {tree_file}, file already exists.")

                elif mode_tree == "fasttree": # Explizit "fasttree" verwenden
                    tree_file = os.path.join(tree_output_folder, f"{name}{number}{i}_FastTree_tree.tre")
                    
                    if not os.path.isfile(tree_file):
                        cmd = [fasttree, "-wag", "-nopr", cln_aln] # FastTree kann direkt aus einer Datei lesen
                        print(f"Executing FastTree command: {' '.join(cmd)}")
                        with open(tree_file, 'w') as outfile:
                            p = subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE)
                        tree_processes.append(p)
                        cpu_use += cpu_par
                        print(f"Started FastTree construction {i+1} with PID {p.pid}")
                    else:
                        print(f"Skipping FastTree for {tree_file}, file already exists.")
                else:
                    print(f"Unsupported tree mode: {mode_tree}. Skipping tree construction for chunk {i+1}.")
                    continue # Springe zur nächsten Iteration

                if tree_file: # Füge nur hinzu, wenn ein tree_file Name generiert wurde
                    trees.append(tree_file)
                        
            except Exception as e:
                print(f"Error starting tree construction {i+1}: {str(e)}")
                import traceback
                traceback.print_exc() # Für detailliertere Fehlermeldungen

        # Wait for all tree construction processes to complete
        for p in tree_processes:
            try:
                stdout, stderr = p.communicate(timeout=3600) # Wait with timeout (e.g., 1 hour)
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
        
        print(f"Successfully started {len(trees)} phylogenetic tree processes (some may have failed).")
        
        # Filter out trees that weren't actually created due to errors
        existing_trees = [tree for tree in trees if os.path.isfile(tree) and os.path.getsize(tree) > 0]
        print(f"Found {len(existing_trees)} actually existing tree files.")

        return existing_trees
        
    except Exception as e:
        print(f"Error in parallel_tree_constructor: {str(e)}")
        import traceback
        traceback.print_exc() # Für detailliertere Fehlermeldungen
        raise

# def replace_spaces_in_tree(tree_file):
#     """
#     Replace spaces with underscores in leaf names of a Newick tree file without using BioPython.
    
#     Args:
#         tree_file: Path to the input tree file
        
#     Returns:
#         Path to the modified tree file (same as input)
#     """
#     with open(tree_file, 'r') as f:
#         tree_str = f.read()
    
#     def replace_spaces(match):
#         leaf_name = match.group(1)
#         branch_length = match.group(2) if match.group(2) else ''
#         return leaf_name.replace(' ', '_') + branch_length
    
#     # Pattern explanation:
#     # ([^()\[\]:,]+) - captures leaf names (anything that's not tree structure characters)
#     # (:\d*\.?\d*)? - optionally captures branch lengths that might follow
#     modified_tree = re.sub(r'([^()\[\]:,]+)(:\d*\.?\d*)?', replace_spaces, tree_str)
    
#     with open(tree_file, 'w') as f:
#         f.write(modified_tree)
    
#     return tree_file

def replace_spaces_in_tree(tree_file):
    """
    Replace spaces with underscores in leaf names of a Newick tree file without using BioPython.
    
    Args:
        tree_file: Path to the input tree file
        
    Returns:
        Path to the modified tree file (same as input)
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
def create_in_out_anno(baits_path: str, outgroup_path: str) -> Tuple[List[str], List[str]]:
    """Create in_list and out_list from bait and outgroup FASTA files.
    
    Args:
        baits_path: Path to the bait sequences FASTA file
        outgroup_path: Optional path to the outgroup sequences FASTA file
        
    Returns:
        Tuple containing:
        - in_list: List of bait sequence IDs
        - out_list: List of outgroup sequence IDs (empty list if no outgroup)
        
    Raises:
        ValueError: If baits file cannot be read or is empty
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
    
    # --- preparation of data structure --- #
    groups_around_ref_gene = {}
    for gene in (in_list + out_list):
        groups_around_ref_gene.update({gene: []})
    
    # --- find node objects of reference genes --- #
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
                    elif seq_label in present_out_list:  # Changed from 'if' to 'elif'
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
    
    print(f"Successfully classified {len(results)} sequences from tree {os.path.basename(tree_file)}")
    return results

# def candidates_baits_assignment(baits_path, tree_file, member_candidates, bait_groups, reroot_id = None, patristic_threshold=1.0):
#     """Assign new candidates to reference members with group hierarchy"""
#     ref_members = read_fasta(baits_path)

#     # Initialize data structures
#     new2ref_mapping_table = {gene: [] for gene in member_candidates}
#     new_per_ref_mem = {gene: [] for gene in ref_members.keys()}
    
#     # Load tree and reroot FIRST
#     tree = dendropy.Tree.get_from_path(tree_file, "newick")
    
#     # Rerooting BEFORE creating distance matrix
#     # if reroot_id is not None:
#     #     node_to_reroot = tree.find_node_with_taxon_label(reroot_id)
#     #     if node_to_reroot is not None:
#     #         tree.reroot_at_node(node_to_reroot)
#     #     else:
#     #         print(f"Warning: Reroot node '{reroot_id}' not found, using midpoint rooting")
#     #         tree.reroot_at_midpoint()
#     # else:
#     tree.reroot_at_midpoint()

#     # Create distance matrix AFTER rerooting
#     pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
    
#     # Get node objects AFTER rerooting
#     ref_node_objects = {}
#     new_node_objects = {}
    
#     for node in tree.taxon_namespace:
#         if node.label in ref_members:
#             ref_node_objects[node.label] = node
#         elif node.label in member_candidates:
#             new_node_objects[node.label] = node
    
#     # Prepare node lists
#     ref_nodes = [ref_node_objects.get(gene) for gene in ref_members.keys() if gene in ref_node_objects]
#     candidate_nodes = [new_node_objects[gene] for gene in member_candidates]
    
#     # Assign candidates to baits with group hierarchy
#     for candidate_node in candidate_nodes:
#         # Berechne Distanzen zu allen Referenzknoten
#         distances = []
#         for ref_node in ref_nodes:
#             try:
#                 edge_dist = pdm.path_edge_count(candidate_node, ref_node)
#                 patr_dist = pdm.patristic_distance(candidate_node, ref_node)
#                 distances.append({
#                     'ref_id': ref_node.label,
#                     'group': bait_groups.get(ref_node.label, 'unknown'),
#                     'edges': edge_dist,
#                     'patr': patr_dist
#                 })
#             except KeyError as e:
#                 print(f"Warning: Could not calculate distance between {candidate_node.label} and {ref_node.label}: {e}")
#                 continue

#         # Sortiere nach Kantenlänge
#         distances.sort(key=lambda x: x['edges'])

#         # Suche gültige Zuordnungen unterhalb des Thresholds
#         assignments = []
#         groups_found = set()

#         for dist in distances:
#             if dist['patr'] <= patristic_threshold:
#                 ref_id = dist['ref_id']
#                 group = dist['group']

#                 if group.lower() == 'homology' and 'homology' not in groups_found:
#                     assignments.append(dist)
#                     groups_found.add('homology')
#                 elif group.lower() == 'transcript' and 'transcript' not in groups_found and 'homology' in groups_found:
#                     assignments.append(dist)
#                     groups_found.add('transcript')
#                 elif group.lower() == 'protein' and 'protein' not in groups_found and 'transcript' in groups_found:
#                     assignments.append(dist)
#                     groups_found.add('protein')
#                 elif group.lower() == 'literature' and 'literature' not in groups_found:
#                     assignments.append(dist)
#                     groups_found.add('literature')
#                     break

#             if 'literature' in groups_found:
#                 break

#         # Fallback: Keine gültige Referenz unter Threshold
#         if not assignments:
#             assignments = [{
#                 'ref_id': 'No match in reference sequences below patristic distance threshold',
#                 'group': 'None',
#                 'edges': 'None',
#                 'patr': 'None'
#             }]

#         # Speichere Zuordnung für Kandidaten
#         new2ref_mapping_table[candidate_node.label] = assignments

#         # Rückrichtung nur für echte Zuordnungen unterhalb des Thresholds
#         for assignment in assignments:
#             if isinstance(assignment.get('patr'), (int, float)) and assignment['patr'] <= patristic_threshold:
#                 new_per_ref_mem[assignment['ref_id']].append({
#                     'candidate_id': candidate_node.label,
#                     'edges': assignment['edges'],
#                     'patr': assignment['patr']
#                 })


    
#     return new2ref_mapping_table, new_per_ref_mem
            
def candidates_baits_assignment(baits_path, tree_file, member_candidates, bait_groups, reroot_id=None, patristic_threshold=1.0, candidate_threshold=5.0):
    """Assign new candidates to reference members with group hierarchy"""
    ref_members = read_fasta(baits_path)

    new2ref_mapping_table = {gene: [] for gene in member_candidates}
    new_per_ref_mem = {gene: [] for gene in ref_members.keys()}
    neighbour_members = {}

    tree = dendropy.Tree.get_from_path(tree_file, "newick")
    # Rerooting BEFORE creating distance matrix
    # if reroot_id is not None:
    #     node_to_reroot = tree.find_node_with_taxon_label(reroot_id)
    #     if node_to_reroot is not None:
    #         tree.reroot_at_node(node_to_reroot)
    #     else:
    #         print(f"Warning: Reroot node '{reroot_id}' not found, using midpoint rooting")
    #         tree.reroot_at_midpoint()
    # else:
    tree.reroot_at_midpoint()
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)

    ref_node_objects = {}
    new_node_objects = {}

    for node in tree.taxon_namespace:
        if node.label in ref_members:
            ref_node_objects[node.label] = node
        elif node.label in member_candidates:
            new_node_objects[node.label] = node

    ref_nodes = [ref_node_objects.get(gene) for gene in ref_members.keys() if gene in ref_node_objects]
    candidate_nodes = [new_node_objects[gene] for gene in member_candidates]

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

        distances.sort(key=lambda x: x['edges'])

        assignments = []
        groups_found = set()

        for dist in distances:
            if dist['patr'] <= patristic_threshold:
                group = dist['group'].lower()
                if group == 'homology' and 'homology' not in groups_found:
                    assignments.append(dist)
                    groups_found.add('homology')
                elif group == 'transcript' and 'transcript' not in groups_found and 'homology' in groups_found:
                    assignments.append(dist)
                    groups_found.add('transcript')
                elif group == 'protein' and 'protein' not in groups_found and 'transcript' in groups_found:
                    assignments.append(dist)
                    groups_found.add('protein')
                elif group == 'literature' and 'literature' not in groups_found:
                    assignments.append(dist)
                    groups_found.add('literature')
                    break
            if 'literature' in groups_found:
                break

        if not assignments:
            assignments = [{
                'ref_id': 'No match in reference sequences below patristic distance threshold',
                'group': 'None',
                'edges': 'None',
                'patr': 'None'
            }]

        new2ref_mapping_table[candidate_node.label] = assignments

        for assignment in assignments:
            if isinstance(assignment.get('patr'), (int, float)) and assignment['patr'] <= patristic_threshold:
                new_per_ref_mem[assignment['ref_id']].append({
                    'candidate_id': candidate_node.label,
                    'edges': assignment['edges'],
                    'patr': assignment['patr']
                })

        for dist in distances:
            if dist['patr'] <= candidate_threshold:
                ref_id = dist['ref_id']
                neighbour_members.setdefault(ref_id, []).append({
                    'candidate_id': candidate_node.label,
                    'edges': dist['edges'],
                    'patr': dist['patr']
                })

    return new2ref_mapping_table, new_per_ref_mem, neighbour_members

#endregion

# region Domain Check
def domain_check(
    sequences_file: str, 
    domains_file: str, 
    hmmresult: str, 
    hmmoutput: str, 
    hmmsearch: str = "hmmsearch", 
    hmm_score: float = 100.0
) -> Tuple[Dict, Dict, List]:
    """Screen sequences for domain hits using HMMER
    
    Args:
        sequences_file: Path to input sequences FASTA file
        domains_file: Path to HMM domains file
        hmmresult: Path for domain table output
        hmmoutput: Path for HMMER output
        hmmsearch: Path to hmmsearch executable
        hmm_score: Minimum score threshold
        
    Returns:
        Tuple of (results_dict, best_domain_matches, domain_names)
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



# def load_ortholog_family(ref_mapping_file):
#     candidate_to_family = {}
#     if not os.path.isfile(ref_mapping_file):
#         return candidate_to_family

#     with open(ref_mapping_file) as f:
#         header = f.readline().strip().split("\t")

#         # Spaltenindizes ermitteln
#         idx_new_member = header.index("NewMember")
#         idx_ortho_family = header.index("Ortholog_Family") if "Ortholog_Family" in header else None
#         idx_neigh_family = header.index("Neighbour_Family") if "Neighbour_Family" in header else None

#         current_new_member = None  # Zum Zwischenspeichern, falls Zeile leer ist

#         for line in f:
#             cols = line.strip().split("\t")

#             # Leere Zeilen überspringen
#             if not any(cols):
#                 continue

#             # Kandidaten-ID extrahieren (auch bei Wiederholungen mit Leerzeichen)
#             new_member_raw = cols[idx_new_member].strip()
#             if new_member_raw:
#                 current_new_member = new_member_raw

#             if not current_new_member:
#                 continue  # kann nicht zugeordnet werden

#             # Zuerst versuchen: Ortholog-Family
#             family = ""
#             if idx_ortho_family is not None and idx_ortho_family < len(cols):
#                 family = cols[idx_ortho_family].strip()

#             # Falls keine Ortholog-Family vorhanden, versuche Nachbar-Family
#             if not family and idx_neigh_family is not None and idx_neigh_family < len(cols):
#                 family = cols[idx_neigh_family].strip()

#             if family:  # Nur speichern, wenn Family existiert
#                 candidate_to_family[current_new_member] = family

#     return candidate_to_family

def load_ortholog_family(ref_mapping_file):
    candidate_to_family = {}
    if not os.path.isfile(ref_mapping_file):
        return candidate_to_family

    with open(ref_mapping_file) as f:
        header = f.readline().strip().split("\t")

        # Spaltenindizes
        idx_new_member = header.index("NewMember")
        idx_ortho_family = header.index("Ortholog_Family") if "Ortholog_Family" in header else None
        idx_neigh_family = header.index("Neighbour_Family") if "Neighbour_Family" in header else None

        current_new_member = None

        for line in f:
            cols = line.strip().split("\t")

            if not any(cols):
                continue

            # Hole NewMember (ggf. erneut verwenden, wenn Zeile leer ist)
            new_member = cols[idx_new_member].strip() if idx_new_member < len(cols) else ""
            if new_member:
                current_new_member = new_member

            if not current_new_member:
                continue

            # Hole Ortholog_Family
            ortho_family = cols[idx_ortho_family].strip() if idx_ortho_family is not None and idx_ortho_family < len(cols) else ""
            neigh_family = cols[idx_neigh_family].strip() if idx_neigh_family is not None and idx_neigh_family < len(cols) else ""

            # Wenn Ortholog_Family ungültig ist (leer oder "None"), verwende Neighbour_Family
            family = ortho_family if (ortho_family and ortho_family != "None" and ortho_family != "-") else neigh_family

            # Nur speichern, wenn ein Family-Wert gefunden wurde
            if family:
                candidate_to_family[current_new_member] = family

    return candidate_to_family



#endregion

# region Motif Check
def static_motif_check(fasta_file: str, motifs_path: str) -> Dict[str, Dict[str, str]]:
    """
    Überprüft Proteinsequenzen auf das Vorhandensein spezifizierter Motive,
    wobei Reihenfolge und Positionsgrenzen berücksichtigt werden.
    Gibt ein Dictionary analog zur Struktur der Funktion `motif_check` zurück.
    """

    # --- Schritt 1: Sequenzen und Motiv-Definitionen einlesen ---
    sequences = read_fasta(fasta_file)

    motifs = []
    motif_names_list = []
    with open(motifs_path, 'r') as mfile:
        for line in mfile.readlines()[1:]:  # Header überspringen
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

    # --- Schritt 2: Verarbeitung pro Sequenz ---
    results = {}
    for seq_id, sequence in sequences.items():
        # Initialisiere Ergebnisstruktur für diese Sequenz
        seq_result = {motif['name']: "" for motif in motifs}

        # Finde alle erlaubten Treffer pro Motiv
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

        # Backtracking zur Motivabdeckung in richtiger Reihenfolge
        valid_paths = []

        def backtrack(current_path: List[Tuple[str]], last_end_pos: int, remaining_motifs: List[Dict]):
            if not remaining_motifs:
                valid_paths.append(current_path)
                return

            current_motif = remaining_motifs[0]
            motif_name = current_motif['name']

            # Option 1: Motiv überspringen
            backtrack(current_path + [(motif_name, "")], last_end_pos, remaining_motifs[1:])

            # Option 2: Treffer verwenden, wenn Position passt
            for start, end, match in all_matches[motif_name]:
                if start >= last_end_pos:
                    backtrack(current_path + [(motif_name, match)], end, remaining_motifs[1:])

        backtrack([], 0, motifs)

        # Wähle besten Pfad: Maximale Anzahl nicht-leerer Treffer
        if valid_paths:
            best_path = max(valid_paths, key=lambda p: sum(1 for _, m in p if m))
            for motif_name, match in best_path:
                if match:
                    seq_result[motif_name] = match

        results[seq_id] = seq_result

    return results, motif_names_list

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
    
    Args:
        sequences_file: Path to input sequences FASTA file
        motifs_file: Path to HMM motifs file
        hmmresult: Path for motif table output
        hmmoutput: Path for HMMER output
        hmmsearch: Path to hmmsearch executable
        cEvalue: Maximum E-value cutoff
        
    Returns:
        Tuple of (results_dict, motif_names)
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
) -> str:
    """
    Erstellt ein kombiniertes HMM-Profil aus Motiv-Definitionen und Bait-Sequenzen.

    Der Prozess umfasst:
    1.  Einlesen von Motiv-Definitionen (Name, Muster, Score-Threshold, Positionen).
    2.  Einlesen von Bait-Proteinsequenzen (FASTA).
    3.  Suche nach den Motiven in den Bait-Sequenzen mittels BLOSUM62-Scoring.
    4.  Extraktion der gefundenen Motiv-Sequenzen inklusive 3 flankierender Aminosäuren.
    5.  Erstellung eines Multiplen Sequenzalignments für jede Gruppe von Motiv-Sequenzen (MAFFT oder MUSCLE).
    6.  Erstellung eines HMM-Profils aus jedem Alignment (`hmmbuild`).
    7.  Zusammenführung aller HMM-Profile in eine einzige Datei.
    8.  Aufräumen aller temporären Dateien.

    Args:
        protein_motifs_path: Pfad zur tab-getrennten Motiv-Definitionsdatei.
        CYP_source: Pfad zur FASTA-Datei mit den Bait-Sequenzen.
        output_dir: Verzeichnis, in dem temporäre Dateien und die finale HMM-Datei gespeichert werden.
        mafft_available: Flag, ob MAFFT verfügbar ist.
        muscle_available: Flag, ob MUSCLE verfügbar ist.
        mafft_path: Pfad zur MAFFT-Executable.
        muscle_path: Pfad zur MUSCLE-Executable.
        hmmbuild_path: Pfad zur hmmbuild-Executable.

    Returns:
        Der Pfad zur finalen, kombinierten HMM-Datei.
        
    Raises:
        FileNotFoundError: Wenn eine benötigte Datei nicht gefunden wird.
        RuntimeError: Wenn weder MAFFT noch MUSCLE verfügbar ist oder ein `subprocess`-Aufruf fehlschlägt.
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
            if b == 'X': continue
            score += BLOSUM62.get((a, b), BLOSUM62.get((b, a), -4))
        return score

    # --- Hauptlogik ---
    
    # Temporäres Verzeichnis für alle Zwischenschritte
    temp_dir = os.path.join(output_dir, "temp_motif_build")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    try:
        # --- Schritt 1: Eingabedateien parsen ---
        print("Step 1: Parsing motif definitions and bait sequences...")
        motifs = []
        with open(protein_motifs_path, 'r') as mfile:
            for line in mfile.readlines()[1:]: # Header überspringen
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        name, pattern, thr, min_p, max_p = parts[:5]
                        motifs.append((name, pattern, float(thr), int(min_p), int(max_p)))

        bait_sequences = read_fasta(CYP_source)
        motif_records = {name: [] for name, _, _, _, _ in motifs}

        # --- Schritt 2: Motive in Bait-Sequenzen suchen und extrahieren ---
        print("Step 2: Searching for motifs in bait sequences...")
        for seq_id, seq in bait_sequences.items():
            seq = seq.upper()
            current_pos = 0
            for name, pattern, threshold, min_pos, max_pos in motifs:
                mlen = len(pattern.replace('(', '').replace(')', '').replace('/', ''))
                best_score, best_seq, best_start = -float('inf'), '', -1

                motif_variants = _generate_combinations(pattern)
                
                search_start = max(current_pos, min_pos - 1)
                search_end = min(len(seq) - mlen, max_pos - mlen + 1)

                for i in range(search_start, search_end):
                    window = seq[i:i+mlen]
                    for variant in motif_variants:
                        sc = _score_window(window, variant)
                        if sc >= threshold and sc > best_score:
                            best_score, best_seq, best_start = sc, window, i
                
                if best_start != -1:
                    start_flank = max(0, best_start - flanklength)
                    end_flank = min(len(seq), best_start + mlen + flanklength)
                    frag = seq[start_flank:end_flank]
                    header = f"{seq_id}_{name}_{best_start+1}-{best_start+mlen}"
                    motif_records[name].append((header, frag))
                    current_pos = best_start + mlen
        
        # --- Schritt 3-6: Für jedes Motiv FASTA -> Align -> HMM ---
        individual_hmm_files = []
        for name, records in motif_records.items():
            if len(records) < 2: # Alignment/HMM nicht sinnvoll für < 2 Sequenzen
                print(f"  Skipping HMM build for motif '{name}': only {len(records)} sequence(s) found.")
                continue

            print(f"Step 3-6: Processing motif '{name}'...")
            # Temporäre Dateipfade
            fasta_path = os.path.join(temp_dir, f"{name}.fasta")
            aln_path = os.path.join(temp_dir, f"{name}.aln")
            hmm_path = os.path.join(temp_dir, f"{name}.hmm")

            # Schritt 3: FASTA-Datei für das Motiv schreiben
            with open(fasta_path, 'w') as f:
                for header, seq in records:
                    f.write(f">{header}\n{seq}\n")

            # Schritt 4: Alignment durchführen
            print(f"  Aligning sequences for '{name}'...")
            if mafft_available:
                cmd = [mafft_path, "--auto", "--quiet", fasta_path]
                with open(aln_path, 'w') as out_file:
                    subprocess.run(cmd, stdout=out_file, check=True)
            elif muscle_available:
                cmd = [muscle_path, "-align", fasta_path, "-output", aln_path]
                subprocess.run(cmd, check=True, capture_output=True)
            else:
                raise RuntimeError("No alignment tool (MAFFT or MUSCLE) is available.")

            # Schritt 5: HMM erstellen
            print(f"  Building HMM for '{name}'...")
            cmd = [hmmbuild_path, hmm_path, aln_path]
            subprocess.run(cmd, check=True, capture_output=True)
            
            individual_hmm_files.append(hmm_path)

        # --- Schritt 7: Alle HMMs zusammenführen ---
        if not individual_hmm_files:
            raise RuntimeError("No HMMs could be built. Check motif definitions and bait sequences.")
            
        print("Step 7: Combining all created HMMs...")
        combined_hmm_path = os.path.join(output_dir, "motifs_from_baits.hmm")
        with open(combined_hmm_path, 'wb') as outfile:
            for hmm_file in individual_hmm_files:
                with open(hmm_file, 'rb') as infile:
                    outfile.write(infile.read())

        print(f"Successfully created combined HMM file: {combined_hmm_path}")
        return combined_hmm_path

    except subprocess.CalledProcessError as e:
        # Fehler von MAFFT, MUSCLE, oder hmmbuild abfangen
        print(f"Error during subprocess execution: {e.cmd}")
        print(f"Return Code: {e.returncode}")
        print(f"Stderr: {e.stderr.decode() if e.stderr else 'N/A'}")
        raise RuntimeError(f"A required tool failed to execute: {e.cmd}")
    finally:
        # --- Schritt 8: Temporäres Verzeichnis aufräumen ---
        if os.path.exists(temp_dir):
            print("Step 8: Cleaning up temporary files...")
            shutil.rmtree(temp_dir)

def parse_motif_names_lengths(hmm_file):
    names = []
    lengths = {}
    with open(hmm_file) as f:
        current_motif = None
        for line in f:
            if line.startswith("NAME"):
                current_motif = line.strip().split()[1]
                names.append(current_motif)
            elif line.startswith("LENG") and current_motif:
                lengths[current_motif] = int(line.strip().split()[1])
    return names, lengths
# endregion

# region ExpressionMatrix
def Candidate_Expression(expression_matrix, candidates, min_avg_tpm, min_single_tpm):
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

# def use_metadata(filtered_expression_file, metadata_path):
#     metadata_dict = {}
#     metadata_ids = set()
#     header_positions = {}

#     # Einlesen der Metadaten
#     with open(metadata_path, 'r', newline='') as meta_file:
#         reader = csv.reader(meta_file)
#         meta_header = next(reader)
        
#         # Suche nach 'tissue' und allen Spalten mit 'condition'
#         for idx, col_name in enumerate(meta_header):
#             if col_name == "tissue" or col_name.startswith("condition"):
#                 header_positions[col_name] = idx
        
#         # Restliche Metadaten einlesen
#         for row in reader:
#             sample_id = row[0]
#             metadata_ids.add(sample_id)
#             metadata_dict[sample_id] = [row[idx] for col, idx in header_positions.items()]
    
#     # Einlesen der gefilterten Expressionsmatrix
#     with open(filtered_expression_file, 'r') as f:
#         lines = f.readlines()

#     if not lines:
#         print("Warning: Filtered expression matrix is empty")
#         return

#     header_parts = lines[0].strip().split('\t')
#     sample_ids = header_parts[1:]  # ab zweiter Spalte
#     new_sample_ids = []

#     matched_indices = []
#     for idx, sid in enumerate(sample_ids):
#         if sid in metadata_dict:
#             new_id = f"{sid}(" + ", ".join(metadata_dict[sid]) + ")"
#             new_sample_ids.append(new_id)
#             matched_indices.append(idx + 1)  # +1 wegen Tab-Datei
#         else:
#             new_sample_ids.append(f"{sid}(no metadata)")

#     # Neue Dateien schreiben
#     metadata_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix")
#     metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")

#     with open(metadata_expression_file, 'w') as out_all, open(metadata_only_expression_file, 'w') as out_meta:
#         out_all.write(header_parts[0] + '\t' + '\t'.join(new_sample_ids) + '\n')
#         out_meta.write(header_parts[0] + '\t' + '\t'.join([new_sample_ids[i - 1] for i in matched_indices]) + '\n')

#         for line in lines[1:]:
#             parts = line.strip().split('\t')
#             gene_id = parts[0]
#             expr_values = parts[1:]
#             out_all.write(gene_id + '\t' + '\t'.join(expr_values) + '\n')
#             meta_values = [expr_values[i - 1] for i in matched_indices]
#             out_meta.write(gene_id + '\t' + '\t'.join(meta_values) + '\n')

#     print(f"Metadata-annotated expression matrix written to: {metadata_expression_file}")
#     print(f"Metadata-only expression matrix written to: {metadata_only_expression_file}")

def count_commas(line):
    in_quotes = False
    comma_count = 0
    for char in line:
        if char == '"':
            in_quotes = not in_quotes
        elif char == ',' and not in_quotes:
            comma_count += 1
    return comma_count

def use_metadata(filtered_expression_file, metadata_path):
    metadata_dict = {}
    metadata_ids = set()
    header_positions = {}

    unused_metadata_path = os.path.join(os.path.dirname(filtered_expression_file), "unused_metadata.txt")
    unused_lines = []

    # Einlesen der Metadaten
    with open(metadata_path, 'r', newline='') as meta_file:
        all_lines = meta_file.readlines()

    if len(all_lines) < 2:
        print("Metadata file is empty or has no data rows.")
        return

    meta_header = next(csv.reader([all_lines[0]]))
    
    # Suche nach 'tissue' und allen Spalten mit 'condition'
    for idx, col_name in enumerate(meta_header):
        if col_name == "tissue" or col_name.startswith("condition"):
            header_positions[col_name] = idx

    # Erwartete Anzahl an Kommata (erste Datenzeile ohne Header)
    expected_commas = count_commas(all_lines[1])

    # Datenzeilen verarbeiten
    for line in all_lines[1:]:
        if count_commas(line) != expected_commas:
            unused_lines.append(line.strip())
            continue
        row = next(csv.reader([line]))
        sample_id = row[0]
        metadata_ids.add(sample_id)
        metadata_dict[sample_id] = [row[idx] for col, idx in header_positions.items()]

    # Unbenutzte Metadaten schreiben
    if unused_lines:
        with open(unused_metadata_path, 'w') as out_unused:
            out_unused.write("These metadata lines had inconsistent comma counts (outside of quotes) and were ignored:\n")
            for l in unused_lines:
                out_unused.write(l + '\n')
        print(f"Skipped {len(unused_lines)} inconsistent metadata lines. See: {unused_metadata_path}")

    # Einlesen der gefilterten Expressionsmatrix
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
            matched_indices.append(idx + 1)
        else:
            new_sample_ids.append(f"{sid}(no metadata)")

    # Neue Dateien schreiben
    metadata_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix")
    metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")

    with open(metadata_expression_file, 'w') as out_all, open(metadata_only_expression_file, 'w') as out_meta:
        out_all.write(header_parts[0] + '\t' + '\t'.join(new_sample_ids) + '\n')
        out_meta.write(header_parts[0] + '\t' + '\t'.join([new_sample_ids[i - 1] for i in matched_indices]) + '\n')

        for line in lines[1:]:
            parts = line.strip().split('\t')
            gene_id = parts[0]
            expr_values = parts[1:]
            out_all.write(gene_id + '\t' + '\t'.join(expr_values) + '\n')
            meta_values = [expr_values[i - 1] for i in matched_indices]
            out_meta.write(gene_id + '\t' + '\t'.join(meta_values) + '\n')

    print(f"Metadata-annotated expression matrix written to: {metadata_expression_file}")
    print(f"Metadata-only expression matrix written to: {metadata_only_expression_file}")

def write_highest_expression_info(filtered_expression_file, metadata_expression_file=None, metadata_only_expression_file=None):
    import re

    output_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_highest_expression")

    # Schritt 1: Originale Expressionen lesen
    with open(filtered_expression_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().split('\t')
    sample_ids = header[1:]
    gene_rows = [line.strip().split('\t') for line in lines[1:]]

    # Schritt 2: Metadaten aus annotierter Datei holen (falls vorhanden)
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

    # Schritt 3: Falls metadata_only_expression_file vorhanden, lesen wir es für AltSample+TPM
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

    # Schritt 4: Ausgabe vorbereiten
    with open(output_file, 'w') as out:
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

    print(f"Highest expression summary written to: {output_file}")

    return output_file

#endregion

# region Paralogs
"""Both functions from bHLH-Annotator"""
def establish_paralog_groups( tree_file, member_candidates, dist_cutoff_factor ):
    """! @brief construct paralog groups """
    """From bHLH-Annotator"""
    
    candidate_mapping_table = {}
    for gene in member_candidates:    #candidate genes of new species
        candidate_mapping_table.update( { gene: None } )
    
    # --- find node objects of reference genes --- #
    tree = dendropy.Tree.get_from_path( tree_file, "newick" )
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
    my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
    
    new_node_objects = {}    #get new family candidate node objects
    for node in tree.taxon_namespace:
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
    """Paraloggruppen erstellen unter Berücksichtigung von tissue/condition-spezifischer Expression."""

    # (1) Lese Expressionen + Metadaten
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

    # (2) Original: Baum laden & Kandidaten identifizieren
    candidate_mapping_table = {gene: None for gene in member_candidates}
    tree = dendropy.Tree.get_from_path(tree_file, "newick")
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
    mean_dist = pdm.mean_nearest_taxon_distance()

    new_node_objects = {node.label: node for node in tree.taxon_namespace if node.label in candidate_mapping_table}
    candidate_gene_nodes = [new_node_objects[gene] for gene in member_candidates]

    # (3) Modifiziertes Gruppieren
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
                # --- Expression-basierter Filter --- #
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

                    # Wenn es Bedingungen gibt, in denen nur jeweils einer der beiden exprimiert ist
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
    """Wählt längste Sequenz pro Gruppe & ergänzt bei Bedarf tissue/condition-Metadaten."""

    import re
    expression_annotation = {}

    if metadata_expression_path:
        with open(metadata_expression_path, 'r') as f:
            header = f.readline().strip().split('\t')
            sample_ids = header[1:]
            for sid in sample_ids:
                match = re.search(r"\((.*?)\)", sid)
                expression_annotation[sid] = match.group(1) if match else "NA"

    paralog_representatives = {}
    annotation_dict = {}  # für Ausgabe

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

    # Erweiterung der 08_representative_paralogs.txt
    # if metadata_expression_path:
    #     meta_output_file = repr_clean_file.replace(".fasta", ".txt")
    #     with open(meta_output_file, "w") as out:
    #         out.write("RepresentativeSeqID\tMembersOfParalogGroup\n")
    #         for rep, group in paralog_representatives.items():
    #             meta_info = f"({','.join(annotation_dict.get(rep, []))})" if rep in annotation_dict else ""
    #             out.write(f"{rep}{meta_info}\t" + ";".join(group for group in paralog_groups if rep in group[0])[0] + "\n")

    return paralog_representatives
#endregion

# region Summary
def collect_summary_data(
        ref_mapping_file,
        baits_info,
        best_domain_match_file=None,
        motif_check_file_summary=None,
        highest_expression_file=None,
        paralog_group_file=None
        ):
    
    summary_data = {}

    # 1. ref_mapping_file: Column 2
    with open(ref_mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            seq_id = parts[1]
            gene_id = parts[2]
            summary_data[seq_id] = [seq_id, gene_id]  # Start mit ID und Gene-ID

    for seq_id in summary_data:
        gene_id = summary_data[seq_id][1]  
        baits_info_row = baits_info.get(gene_id, {})
        col4 = baits_info_row.get("Family", "-")
        col5 = baits_info_row.get("Subfamily", "-")
        summary_data[seq_id].extend([col4, col5])

    # 2. best_domain_match_file: Spalte 2 anhand ID aus Spalte 1
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

    # 3. motif_check_file_summary: Spalten 2-7 anhand ID aus Spalte 1
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

    # 4. highest_expression_file: Spalten 2-5 anhand ID aus Spalte 1
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

    # 5. paralog_group_file: ID vor dem ":" → repräsentative Sequenz
    repr_ids = set()
    if paralog_group_file:
        with open(paralog_group_file, 'r') as f:
            for line in f:
                if ":" in line:
                    parts = line.strip().split(":")
                    repr_id = parts[0].strip()
                    repr_ids.add(repr_id)
    for seq_id in summary_data:
        summary_data[seq_id].append("yes" if seq_id in repr_ids else "no")

    return summary_data

#endregion

#region Output
# def write_output_file(
#     file_path: str,
#     content: Union[str, List[str], Dict],
#     header: Optional[Union[str, List[str]]] = None,
#     overwrite: bool = False,
#     mode: str = 'w'
# ) -> None:
#     """
#     Unified function to write various output files with consistent handling.

#     Args:
#         file_path: Path to output file
#         content: Content to write (string, list of lines, or dictionary)
#         header: Optional header line(s)
#         overwrite: Whether to overwrite existing file (default False)
#         mode: Write mode ('w' for write, 'a' for append')

#     Raises:
#         ValueError: If file exists and overwrite=False
#         IOError: If writing fails
#     """
#     if os.path.exists(file_path) and not overwrite and mode == 'w':
#         return

#     try:
#         with open(file_path, mode) as f:
#             # Write header
#             if header is not None:
#                 if isinstance(header, list):
#                     f.write("\t".join(header))
#                 else:
#                     f.write(str(header))
#                 f.write("\n")  # genau ein Zeilenumbruch nach Header

#             # Write content
#             if isinstance(content, str):
#                 f.write(content.rstrip("\n") + "\n")
            
#             elif isinstance(content, list):
#                 for line in content:
#                     if isinstance(line, list):
#                         line_str = "\t".join(map(str, line))
#                     else:
#                         line_str = str(line)
#                     f.write(line_str.rstrip("\n") + "\n")
            
#             elif isinstance(content, dict):
#                 for key, value in sorted(content.items()):
#                     if isinstance(value, (list, tuple)):
#                         line = f"{key}\t" + "\t".join(map(str, value))
#                     else:
#                         line = f"{key}\t{value}"
#                     f.write(str(line).rstrip("\n") + "\n")

#     except Exception as e:
#         raise IOError(f"Failed to write {file_path}: {str(e)}")

    
# def write_fasta_file(file_path: str, sequences: Dict[str, str], overwrite: bool = False) -> None:
#     """Specialized function for writing FASTA files"""
#     if os.path.exists(file_path) and not overwrite:
#         return
        
#     with open(file_path, 'w') as f:
#         for seq_id, sequence in sequences.items():
#             f.write(f">{seq_id}\n")
#             for i in range(0, len(sequence), 80):
#                 f.write(sequence[i:i+80] + "\n")

# def convert_txt_to_html(txt_folder, html_folder, link_overview_file):
#     import os
#     import re
#     from collections import defaultdict

#     os.makedirs(html_folder, exist_ok=True)

#     txt_files = [f for f in os.listdir(txt_folder) if f.endswith(".txt")]
#     link_dict = defaultdict(list)

#     # Group by prefix
#     for txt_file in txt_files:
#         match = re.match(r"(\d+)_", txt_file)
#         group = match.group(1) if match else "misc"
#         link_dict[group].append(txt_file)

#     # Create HTML Files
#     for group, files in link_dict.items():
#         for txt_file in sorted(files):
#             txt_path = os.path.join(txt_folder, txt_file)
#             html_file = txt_file.replace(".txt", ".html")
#             html_path = os.path.join(html_folder, html_file)

#             with open(txt_path, "r") as f_in, open(html_path, "w") as f_out:
#                 lines = f_in.readlines()

#                 f_out.write("<html><body>\n")
#                 f_out.write(f"<h2>{txt_file}</h2>\n")
#                 f_out.write("<table border='1' cellspacing='0' cellpadding='3'>\n")

#                 for i, line in enumerate(lines):
#                     cols = line.strip().split("\t")
#                     f_out.write("<tr>")
#                     for col in cols:
#                         tag = "th" if i == 0 else "td"
#                         f_out.write(f"<{tag}>{col}</{tag}>")
#                     if i == 0:
#                         f_out.write("<th>Overview</th>")
#                     else:
#                         rel_link = os.path.relpath(link_overview_file, html_folder)
#                         f_out.write(f"<td><a href='{os.path.basename(rel_link)}'>Overview</a></td>")
#                     f_out.write("</tr>\n")

#                 f_out.write("</table>\n</body></html>")

#     # Create Overview File
#     with open(link_overview_file, "w") as f_link:
#         f_link.write("<html><body>\n<h1>CYP Output Files</h1>\n")

#         for group in sorted(link_dict.keys(), key=lambda x: int(x) if x.isdigit() else 999):
#             f_link.write(f"<h2>Section {group}</h2>\n<ul>\n")
#             for txt_file in sorted(link_dict[group]):
#                 html_file = txt_file.replace(".txt", ".html")
#                 f_link.write(f"<li><a href='{html_file}'>{html_file}</a></li>\n")
#             f_link.write("</ul><br>\n")

#         f_link.write("</body></html>")

def baits_info_to_txt(baits_info: Dict[str, Dict[str, str]], headers: List[str], output_path: str) -> None:
    """
    Write the baits_info dictionary to a tab-separated .txt file using the original column order.

    :param baits_info: Dictionary from read_baits_info function
    :param headers: Original headers from the CSV (including 'ID' as first column)
    :param output_path: Path to the output .txt file
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
        print(f"Written baits info to {output_path}")
    except Exception as e:
        raise IOError(f"Could not write baits info to {output_path}: {e}")


def convert_txt_to_html(txt_folder, html_folder, link_overview_file):
    import os
    import re
    from collections import defaultdict

    os.makedirs(html_folder, exist_ok=True)

    txt_files = [f for f in os.listdir(txt_folder) if f.endswith(".txt")]
    link_dict = defaultdict(list)

    # Gruppiere Dateien nach Präfix
    for txt_file in txt_files:
        match = re.match(r"(\d+)_", txt_file)
        group = match.group(1) if match else "misc"
        link_dict[group].append(txt_file)

    # Erzeuge HTML-Dateien
    for group, files in link_dict.items():
        for txt_file in sorted(files):
            txt_path = os.path.join(txt_folder, txt_file)
            html_file = txt_file.replace(".txt", ".html")
            html_path = os.path.join(html_folder, html_file)

            with open(txt_path, "r") as f_in, open(html_path, "w") as f_out:
                lines = [line.rstrip('\n') for line in f_in if line.strip() != ""]

                f_out.write("<html><body>\n")
                f_out.write(f"<h2>{txt_file}</h2>\n")

                # Overview-Link unter Überschrift
                rel_link = os.path.relpath(link_overview_file, html_folder)
                f_out.write(f"<p><a href='{os.path.basename(rel_link)}'>Overview</a></p>\n")

                f_out.write("<table border='1' cellspacing='0' cellpadding='3'>\n")

                for i, line in enumerate(lines):
                    cols = line.split("\t")
                    f_out.write("<tr>")
                    for col in cols:
                        tag = "th" if i == 0 else "td"
                        safe_col = col.strip() if col.strip() != "" else "&nbsp;"  # Leere Zelle sichtbar machen
                        f_out.write(f"<{tag}>{safe_col}</{tag}>")
                    f_out.write("</tr>\n")

                f_out.write("</table>\n</body></html>")

    # Erzeuge Übersicht
    with open(link_overview_file, "w") as f_link:
        f_link.write("<html><body>\n<h1>CYP Output Files</h1>\n")

        for group in sorted(link_dict.keys(), key=lambda x: int(x) if x.isdigit() else 999):
            f_link.write(f"<h2>Section {group}</h2>\n<ul>\n")
            for txt_file in sorted(link_dict[group]):
                html_file = txt_file.replace(".txt", ".html")
                f_link.write(f"<li><a href='{html_file}'>{html_file}</a></li>\n")
            f_link.write("</ul><br>\n")

        f_link.write("</body></html>")

def parse_newick(newick):
    """
    Very simple Newick parser (limited support).
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
    if isinstance(tree, str):
        return 1
    return sum(count_leaves(child) for child in tree)

def draw_tree_rectangular(tree, ax, x=0, y=0, dx=1.5, dy=1.5, fontsize=10):
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
                   format: str = "png", layout: str = "Polar", 
                   cand_color: str = "green", bait_color: str = "red"):
    """
    Findet alle .tre-Dateien im result_folder und speichert visualisierte Bäume
    im output_folder mit gewähltem Layout: Rectangular, Polar, Radial.
    Erhöht automatisch die Auflösung für bessere Lesbarkeit.
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
            #ts.show_branch_support = True
            # ts.scale = 120
            # ts.branch_vertical_margin = 0

            layout_mode = layout.lower()
            if layout_mode == "rectangular":
                ts.mode = "r"
            elif layout_mode == "radial":
                ts.mode = "c"
            else:  # Default to polar
                ts.mode = "c"

            # Dynamische Bildgröße basierend auf Blattanzahl
            leaf_count = len(tree.get_leaves())
            width = max(1200, int(leaf_count * 25))  # mehr Platz für Labels
            dpi = 300  # hohe Auflösung

            # Optional größere Schriftgröße (standardmäßig übernimmt ete3 die Skalierung)
            # ts.scale =  200  # Skalierung des gesamten Baums

            output_filename = os.path.splitext(tree_file)[0] + f".{format}"
            output_path = os.path.join(output_folder, output_filename)
                
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
                    else:
                        color = "black"
                    
                    face = TextFace(node.name, fsize=10, fgcolor=color)
                    node.add_face(face, column=0, position="branch-right")

            tree.render(output_path, w=width, units="px", dpi=dpi, tree_style=ts)
            print(f"Tree saved to: {output_path} ({leaf_count} leaves, width={width}px, dpi={dpi})")

        except Exception as e:
            print(f"Error processing {tree_file}: {e}")

def generate_heatmaps_from_expression_data(txt_folder, heatmap_folder, format = "png"):
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

                # Dynamische Schriftgrößen
                n_rows, n_cols = df.shape
                base_fontsize = 10
                x_fontsize = max(4, min(12, base_fontsize * 40 / max(n_cols, 1)))
                y_fontsize = max(4, min(12, base_fontsize * 40 / max(n_rows, 1)))

                # Dynamische Bildgröße
                #fig_width = min(30, max(8, n_cols * 0.3))
                fig_width = min(150, max(10, n_cols * 0.3))
                fig_height = min(30, max(6, n_rows * 0.3))

                plt.figure(figsize=(fig_width, fig_height))
                ax = sns.heatmap(df, cmap="viridis", cbar=True)

                # Achsenbeschriftungen mit dynamischer Schriftgröße
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=x_fontsize)
                ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=y_fontsize)
                plt.title(f"Heatmap: {filename}", fontsize=14)
                plt.tight_layout()

                base_name = os.path.splitext(filename)[0]
                output_path = os.path.join(heatmap_folder, f"{base_name}.{format}")
                plt.savefig(output_path, dpi=300)
                plt.close()
                print(f"Heatmap saved at: {output_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
        else:
            print(f"File not found: {filename} – will be skipped.")
#endregion

# region Main
def main():
    # Parse command-line arguments
    args = Input.parser.parse_args()
    args = Input.validate_args(args)

    # Check required tools
    tools = check_required_tools()
    use_blast = tools['BLAST']
    use_hmmer = tools['HMMER']
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

    # Collect and process all files
    file_paths, processed_paths = collect_files(
        data_folder=args.data,
        subject_paths=[args.subject],
        bait_path=args.baits,
        protein_motifs_path=args.protein_motifs,
        hmm_motifs_path=args.hmm_motifs,
        hmm_domains_path=args.hmm_domains,
        baits_info_path=args.baits_info,
        expression_path=args.expression,
        metadata_path=args.metadata,
        output_dir=args.processed_input_folder or "Processed_Input",
        trim_names=args.trim_names
    )

    # Extract all paths from processed_paths dictionary
    baits_with_outgroup_path = processed_paths['baits_path']
    baits_path = processed_paths['baits_no_outgroup_path']
    transcript_baits_path = processed_paths['transcript_baits_path']
    protein_baits_path = processed_paths['protein_baits_path']
    literature_baits_path = processed_paths['literature_baits_path']
    outgroup_path = processed_paths['outgroup_path']
    subjects_paths = processed_paths['subject_files_paths']
    protein_motifs_path = processed_paths.get('protein_motifs_path')
    hmm_motifs_path = processed_paths.get('hmm_motifs_path')
    hmm_domains_path = processed_paths.get('hmm_domains_path')
    baits_info_path = processed_paths.get('baits_info_path')
    expression_matrix_path = processed_paths.get('expression_path')
    metadata_path = processed_paths.get('metadata_path')

    # Set output folder
    output_folder = args.output_folder if args.output_folder is not None else "Output"
    processed_input_folder = args.processed_input_folder or "Processed_Input"

    # --- separated analyses for each subject --- #             
    for jidx, subject in enumerate(subjects_paths):    #use jidx to generate unique IDs for all jobs
        
        # define Job-ID
        job_ID = subject.split("/")[-1].split(".")[0]
        
        # --- prepare output folder for each job if there are multiple --- #   
        if len(subjects_paths) == 1:
            job_output_folder = output_folder 
        else:
            job_output_folder = output_folder + "/" + str(jidx).zfill(5) + "_" + job_ID + "/"
        
        if not job_output_folder.endswith("/"):
            job_output_folder = job_output_folder + "/" 
        
        if not os.path.exists(job_output_folder):
            os.makedirs(job_output_folder)
        
        # Neue Ordnerstruktur
        result_folder = os.path.join(job_output_folder, "RESULTS/")
        supplement_folder = os.path.join(job_output_folder, "SUPPLEMENTS/")
        
        # Unterordner in RESULTS
        txt_folder = os.path.join(result_folder, "TXT/")
        fasta_folder = os.path.join(result_folder, "FASTA/")
        tree_folder = os.path.join(result_folder, "TREES/")
        html_folder = os.path.join(result_folder, "HTML/")
        heatmap_folder = os.path.join(result_folder, "HEATMAPS/")
        
        # Ordner erstellen
        os.makedirs(result_folder, exist_ok=True)
        os.makedirs(supplement_folder, exist_ok=True)
        os.makedirs(txt_folder, exist_ok=True)
        os.makedirs(fasta_folder, exist_ok=True)
        os.makedirs(tree_folder, exist_ok=True)
        os.makedirs(html_folder, exist_ok=True)
        if expression_matrix_path:
            os.makedirs(heatmap_folder, exist_ok=True)
        
        # Read baits_info_path
        baits_info, info_headers = read_baits_info(baits_info_path)
        bait_groups = {
            seq_id: subdict.get("Type", "-")
            for seq_id, subdict in baits_info.items()
            }
        
        baits_info_file = os.path.join(txt_folder, "00_baits_info.txt")
        baits_info_to_txt(baits_info, info_headers, baits_info_file)

        start_time = datetime.datetime.now()

        # --- 01 Search initial candidates --- #
        seq_search_result_file = os.path.join(supplement_folder, "01_seq_search_results.txt")
        blast_analyze_folder = os.path.join(supplement_folder, "01_blast_seq_search_analysis/")
        os.makedirs(blast_analyze_folder, exist_ok=True)

        # Skip search if result file already exists
        if not os.path.isfile(seq_search_result_file):
            # Initialisierung
            seq_search_results = {}

            # Entscheidung: HMMER oder BLAST?
            use_hmmer = args.use_hmm.upper() in ("Y", "YES")
            use_blast = not use_hmmer  # exklusiv

            # --- HMMER-Weg ---
            if use_hmmer:
                # Pfad zur HMM-Datei
                hmm_file = args.hmm if args.hmm else os.path.join(job_output_folder, "bait.hmm")

                # Falls keine HMM-Datei übergeben wurde: erzeugen
                if not args.hmm:
                    bait_source = literature_baits_path if literature_baits_path else baits_path
                    if not os.path.isfile(hmm_file):  # nur wenn nicht vorhanden
                        try:
                            hmm_build(bait_source, hmm_file)
                        except Exception as e:
                            print(f"Fehler beim Erzeugen der HMM: {str(e)}")
                            raise SystemExit("Abbruch, da keine gültige HMM-Datei vorhanden ist.")

                # Suche ausführen
                hmm_result_file = os.path.join(job_output_folder, "hmm_results.txt")
                if not os.path.isfile(hmm_result_file):
                    run_hmmsearch(
                        hmm_file=hmm_file,
                        subject_file=subject,
                        output_file=hmm_result_file
                    )
                seq_search_results = load_hmmsearch_results(hmm_result_file)

            # --- BLAST-Weg ---
            elif use_blast:
                blast_result_file = os.path.join(supplement_folder, "blast_results.txt")
                if not os.path.isfile(blast_result_file):
                    run_blast_search(
                        query_file=baits_path,
                        subject_file=subject,
                        output_file=blast_result_file,
                        blast_db_folder=blast_analyze_folder,
                        cpub=args.cpub
                    )
                seq_search_results = load_blast_results(
                    blast_result_file,
                    args.simcutp,
                    args.poscutp,
                    args.lencutp,
                    args.bitcutp
                )

            # Ergebnisse speichern mit write_output_file
            with open(seq_search_result_file, 'w') as f:
                for seq_id in seq_search_results.keys():
                    f.write(f"{seq_id}\n")

        else:
            print(f"Load results from {seq_search_result_file}")
            with open(seq_search_result_file) as f:
                seq_search_results = {line.strip(): None for line in f if line.strip()}

        # Sequences laden
        subject_sequences = read_fasta(subject)
        subject_name_mapping_table = load_name_mapping_table(Path(processed_input_folder) / f"{job_ID}_mapping.txt")

        # Kandidaten als FASTA speichern mit write_fasta_file
        candidate_file = os.path.join(fasta_folder, f"{args.name}01_initial_candidates.fasta")
        candidate_sequences = {seq_id: subject_sequences[seq_id] for seq_id in seq_search_results.keys() if seq_id in subject_sequences}
        with open(candidate_file, "w") as out:
            for seq_id in seq_search_results.keys():
                if seq_id in subject_sequences:
                    out.write(f">{seq_id}\n{subject_sequences[seq_id]}\n")

        # --- 02 construct phylogenetic tree and analyze tree file --- #
        tree_output_folder = os.path.join(supplement_folder, f"{args.name}02_in_out_CYP_analysis_trees/")

        in_list, out_list = create_in_out_anno(baits_path, outgroup_path)
        aln_candidate_file = "alignment_candidates.fasta"
        aln_input_file = "alignment_input.fasta"
        aln_file = "alignment_input.fasta.aln"
        cln_aln_file = "alignment_input.fasta.aln.cln"
        
        # hmm motif
        hmmsearch_seq_file = os.path.join(job_output_folder, "02_hmmsearch_sequences.fasta")
        aln_hmmsearch_results_file = os.path.join(job_output_folder, "02_hmmsearch_results.txt")
        aln_hmmsearch_results = {}
        
        if os.path.isfile(aln_hmmsearch_results_file):
            aln_hmmsearch_results = load_hmmsearch_results(aln_hmmsearch_results_file)
        
        # --- 02 first classification --- #
        clean_members_file_f = os.path.join(tree_output_folder, f"{args.name}first_ingroup_CYPs.fasta")
        tmp_result_table_f = os.path.join(tree_output_folder, f"{args.name}first_in_out_CYP_analysis_results.txt")

        if not os.path.isfile(tmp_result_table_f):
            cyp_classification = {}
            sys.stdout.write(f"Number of ingroup CYP baits: {len(in_list)}\n")
            sys.stdout.write(f"Number of outgroup CYP baits: {len(out_list)}\n")
            sys.stdout.flush()

            # construct phylogenetic trees
            tree_files = parallel_tree_constructor(
                args.numprocesscandidates, tree_output_folder, aln_candidate_file, aln_input_file,
                aln_file, cln_aln_file, baits_with_outgroup_path, candidate_file, "first_", "",
                args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree,
                args.cpumax, args.cpur, args.parallel
            )

            for tree_file in tree_files:
                # Replace spaces in leaf names before processing
                modified_tree_file = replace_spaces_in_tree(tree_file)

                classification = split_into_ingroup_and_outgroup(
                    modified_tree_file, in_list, out_list, args.numneighbours,
                    args.neighbourdist, args.minneighbours, aln_hmmsearch_results
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

            # Write files using helper functions
            # write_output_file(
            #     file_path=tmp_result_table_f,
            #     content=result_lines,
            #     header=header,
            #     overwrite=True
            # )
            
            # write_fasta_file(
            #     file_path=clean_members_file_f,
            #     sequences=fasta_content,
            #     overwrite=True
            # )

            with open(clean_members_file_f, "w") as out, open(tmp_result_table_f, "w") as out2:
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

                    if cyp_classification[ candidate ]["score"] > args.minscore:
                            #and ( candidate in aln_hmmsearch_results or not args.filterdomain.upper() in ["YES", "Y"]):                            
                            out.write( ">" + candidate + "\n" + subject_sequences[ candidate ] + "\n" )

        # --- 02 second classification --- #
        clean_members_file_s = os.path.join(tree_output_folder, f"{args.name}first_ingroup_CYPs.fasta")
        tmp_result_table_s = os.path.join(tree_output_folder, f"{args.name}first_in_out_CYP_analysis_results.txt")
        clean_members_file = os.path.join(fasta_folder, f"{args.name}02_ingroup_CYPs.fasta")
        tmp_result_table = os.path.join(txt_folder, f"{args.name}02_in_out_CYP_analysis_results.txt")

        if not os.path.isfile(tmp_result_table):
            cyp_classification = {}
            sys.stdout.write(f"Number of ingroup CYP baits: {len(in_list)}\n")
            sys.stdout.write(f"Number of outgroup CYP baits: {len(out_list)}\n")
            sys.stdout.flush()

            # construct phylogenetic trees
            tree_files = parallel_tree_constructor(
                args.numprocesscandidates, tree_output_folder, aln_candidate_file, aln_input_file,
                aln_file, cln_aln_file, baits_with_outgroup_path, candidate_file, "second_", "",
                args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree,
                args.cpumax, args.cpur, args.parallel
            )

            for tree_file in tree_files:
                # Replace spaces in leaf names before processing
                modified_tree_file = replace_spaces_in_tree(tree_file)
                
                classification = split_into_ingroup_and_outgroup(
                    modified_tree_file, in_list, out_list, args.numneighbours,
                    args.neighbourdist, args.minneighbours, aln_hmmsearch_results
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

            with open(clean_members_file_s, "w") as out, open(tmp_result_table_s, "w") as out2:
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

                    if cyp_classification[candidate]["score"] > args.minscore:
                        #and (candidate in aln_hmmsearch_results or not args.filterdomain.upper() in ["YES", "Y"]):
                        out.write(f">{candidate}\n{subject_sequences[candidate]}\n")

            # Copy final files
            shutil.copyfile(clean_members_file_s, clean_members_file)
            shutil.copyfile(tmp_result_table_s, tmp_result_table)

        else:
            cyp_classification = load_in_out_classification_file(tmp_result_table)

        clean_members = read_fasta(clean_members_file)
        if len(clean_members) < 1:
            sys.exit(f"ERROR: no CYPs detected.")

        # --- 03 construct a final tree --- #
        # construct final tree with all baits
        if not (os.path.isdir(tree_folder) and any(f.startswith("03") for f in os.listdir(tree_folder))):
            tree_output_folder = os.path.join(supplement_folder, f"{args.name}03_final_tree/")
            fin_aln_candidate_file = f"{args.name}03_fin_alignment_candidates.fasta"
            fin_aln_input_file = f"{args.name}03_fin_alignment_input.fasta"
            fin_aln_file = f"{args.name}03_fin_alignment_input.fasta.aln"
            fin_cln_aln_file = f"{args.name}03_fin_alignment_input.fasta.aln.cln"

            tree_file = parallel_tree_constructor(
                    len(clean_members), tree_output_folder, fin_aln_candidate_file, fin_aln_input_file,
                    fin_aln_file, fin_cln_aln_file, baits_path, clean_members_file, "03_final_", "",
                    args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree,
                    args.cpumax, args.cpur, "N"
                )
            
            final_tree_path = replace_spaces_in_tree(tree_file[0])
            shutil.copy(final_tree_path, tree_folder)
        else:
            final_tree_path = os.path.join(tree_folder, f"{args.name}03_final_0_FastTree_tree.tre")

        # --- 04 find closest reference for baits --- #            
        # group_around_ref_file = os.path.join(txt_folder, f"{args.name}04_candidates_group_around_bait.txt")
        # ref_mapping_file = os.path.join(txt_folder, f"{args.name}04_candidate_2_bait_mapping_file.txt")

        # ref_mapping_table, per_ref_members = candidates_baits_assignment(
        #         baits_path, 
        #         final_tree_path, 
        #         clean_members.keys(),
        #         bait_groups,
        #         args.reroot,
        #         args.ortholog_threshold
        #     )

        # if not os.path.isfile(group_around_ref_file) or os.path.getsize(group_around_ref_file) == 0:
            
        #     # Write group around reference file
        #     with open(group_around_ref_file, "w") as out:
        #         out.write("Ref_Member\tRef_Group\tCandidate_Members\n")
        #         for ref_gene in sorted(per_ref_members.keys()):
        #             original_new_names = []
        #             for candidate in per_ref_members[ref_gene]:
        #                 original_new_names.append(subject_name_mapping_table[candidate['candidate_id']])
        #             group = bait_groups.get(ref_gene, 'unknown')
        #             out.write(f"{ref_gene}\t{group}\t{'; '.join(original_new_names)}\n")

        # if not os.path.isfile(ref_mapping_file) or os.path.getsize(ref_mapping_file) == 0:
            
        #     #baits_info = load_baits_info(baits_info_path)

        #     # Write new to reference mapping file
        #     with open(ref_mapping_file, "w") as out:
        #         if baits_info:
        #             out.write("NewMember\tOriginalID\tOrthologs\tOrtholog_Group\tOrtholog_Family\tOrtholog_Subfamily\tOrtholog_Edge_Distance\tOrtholog_Patristic_Distance\n")
        #             ref_name = ""
        #             for new_gene in sorted(ref_mapping_table.keys()):
        #                 for assignment in ref_mapping_table[new_gene]:
        #                     baits_info_row = baits_info.get(assignment['ref_id'], {})
        #                     if assignment['ref_id'] in baits_info.keys():
        #                         out.write("\t".join([
        #                             new_gene if new_gene != ref_name else " " * len(new_gene),
        #                             subject_name_mapping_table[new_gene] if new_gene != ref_name else " " * len(subject_name_mapping_table[new_gene]),
        #                             assignment['ref_id'],
        #                             assignment['group'],
        #                             # baits_info[assignment['ref_id']]['family'],
        #                             # baits_info[assignment['ref_id']]['subfamily'],
        #                             baits_info_row.get("Family", "-"),
        #                             baits_info_row.get("Subfamily", "-"),
        #                             str(assignment['edges']),
        #                             str(assignment['patr'])
        #                         ]) + "\n")
        #                     else:
        #                         out.write("\t".join([
        #                             new_gene if new_gene != ref_name else " " * len(new_gene),
        #                             subject_name_mapping_table[new_gene] if new_gene != ref_name else " " * len(subject_name_mapping_table[new_gene]),
        #                             assignment['ref_id'],
        #                             assignment['group'],
        #                             "None",
        #                             "None",
        #                             str(assignment['edges']),
        #                             str(assignment['patr'])
        #                         ]) + "\n")
        #                     ref_name = new_gene
        #         else:
        #             out.write("NewMember\tOriginalID\tRefMember\tRefGroup\tEdgeDistance\tPatristicDistance\n")
        #             ref_name = ""
        #             for new_gene in sorted(ref_mapping_table.keys()):
        #                 for assignment in ref_mapping_table[new_gene]:
        #                     out.write("\t".join([
        #                         new_gene if new_gene != ref_name else " " * len(new_gene),
        #                         subject_name_mapping_table[new_gene] if new_gene != ref_name else " " * len(subject_name_mapping_table[new_gene]),
        #                         assignment['ref_id'],
        #                         assignment['group'],
        #                         str(assignment['edges']),
        #                         str(assignment['patr'])
        #                     ]) + "\n")
        #                     ref_name = new_gene

        # --- 04 find closest reference for baits --- #            
        group_around_ref_file = os.path.join(txt_folder, f"{args.name}04_candidates_group_around_bait.txt")
        ref_mapping_file = os.path.join(txt_folder, f"{args.name}04_candidate_2_bait_mapping_file.txt")
        neighbour_fasta_file = Path(fasta_folder) / f"{args.name}04_candidate_neighbours.fasta"

        ref_mapping_table, per_ref_members, neighbour_members = candidates_baits_assignment(
            baits_path, 
            final_tree_path, 
            clean_members.keys(),
            bait_groups,
            args.reroot,
            args.ortholog_threshold,
            args.neighbour_threshold
        )

        if not os.path.isfile(group_around_ref_file) or os.path.getsize(group_around_ref_file) == 0:
            with open(group_around_ref_file, "w") as out:
                out.write("Ref_Member\tRef_Group\tClose_Candidate_Members\tCandidate_Members\n")
                for ref_gene in sorted(set(per_ref_members.keys()).union(neighbour_members.keys())):
                    close_members = [subject_name_mapping_table[c['candidate_id']] for c in per_ref_members.get(ref_gene, [])]
                    neighbour_names = [subject_name_mapping_table[c['candidate_id']] for c in neighbour_members.get(ref_gene, [])]
                    group = bait_groups.get(ref_gene, 'unknown')
                    out.write(f"{ref_gene}\t{group}\t{'; '.join(close_members)}\t{'; '.join(neighbour_names)}\n")

        if not os.path.isfile(ref_mapping_file) or os.path.getsize(ref_mapping_file) == 0:
            with open(ref_mapping_file, "w") as out:
                if baits_info:
                    out.write("NewMember\tOriginalID\tOrthologs\tOrtholog_Group\tOrtholog_Family\tOrtholog_Subfamily\tOrtholog_Edge_Distance\tOrtholog_Patristic_Distance\tNeighbour\tNeighbour_Group\tNeighbour_Family\tNeighbour_Subfamily\tNeighbour_Edge_Distance\tNeighbour_Patristic_Distance\n")
                    ref_name = ""
                    for new_gene in sorted(ref_mapping_table.keys()):
                        orthologs = ref_mapping_table[new_gene]
                        neighbours = []
                        for ref_id, candidates in neighbour_members.items():
                            for entry in candidates:
                                if entry['candidate_id'] == new_gene:
                                    neighbours.append({
                                        'ref_id': ref_id,
                                        'group': bait_groups.get(ref_id, 'unknown'),
                                        'edges': entry['edges'],
                                        'patr': entry['patr']
                                    })
                        max_len = max(len(orthologs), len(neighbours))
                        for i in range(max_len):
                            ortho = orthologs[i] if i < len(orthologs) else {}
                            neigh = neighbours[i] if i < len(neighbours) else {}
                            ortho_id = ortho.get('ref_id', '')
                            ortho_grp = ortho.get('group', '')
                            ortho_fam = baits_info.get(ortho_id, {}).get('Family', '-') if ortho_id in baits_info else '-'
                            ortho_subfam = baits_info.get(ortho_id, {}).get('Subfamily', '-') if ortho_id in baits_info else '-'
                            neigh_id = neigh.get('ref_id', '')
                            neigh_grp = neigh.get('group', '')
                            neigh_fam = baits_info.get(neigh_id, {}).get('Family', '-') if neigh_id in baits_info else '-'
                            neigh_subfam = baits_info.get(neigh_id, {}).get('Subfamily', '-') if neigh_id in baits_info else '-'
                            out.write("\t".join([
                                new_gene if new_gene != ref_name else " " * len(new_gene),
                                subject_name_mapping_table[new_gene] if i == 0 else " " * len(subject_name_mapping_table[new_gene]),
                                ortho_id,
                                ortho_grp,
                                ortho_fam,
                                ortho_subfam,
                                str(ortho.get('edges', '')),
                                str(ortho.get('patr', '')),
                                neigh_id,
                                neigh_grp,
                                neigh_fam,
                                neigh_subfam,
                                str(neigh.get('edges', '')),
                                str(neigh.get('patr', ''))
                            ]) + "\n")
                            ref_name = new_gene
                else:
                    out.write("NewMember\tOriginalID\tRefMember\tRefGroup\tEdgeDistance\tPatristicDistance\tNeighbour\tNeighbourGroup\tNeighbourEdgeDistance\tNeighbourPatristicDistance\n")
                    for new_gene in sorted(ref_mapping_table.keys()):
                        orthologs = ref_mapping_table[new_gene]
                        neighbours = []
                        for ref_id, candidates in neighbour_members.items():
                            for entry in candidates:
                                if entry['candidate_id'] == new_gene:
                                    neighbours.append({
                                        'ref_id': ref_id,
                                        'group': bait_groups.get(ref_id, 'unknown'),
                                        'edges': entry['edges'],
                                        'patr': entry['patr']
                                    })
                        max_len = max(len(orthologs), len(neighbours))
                        for i in range(max_len):
                            ortho = orthologs[i] if i < len(orthologs) else {}
                            neigh = neighbours[i] if i < len(neighbours) else {}
                            out.write("\t".join([
                                new_gene if i == 0 else " " * len(new_gene),
                                subject_name_mapping_table[new_gene] if i == 0 else " " * len(subject_name_mapping_table[new_gene]),
                                ortho.get('ref_id', ''),
                                ortho.get('group', ''),
                                str(ortho.get('edges', '')),
                                str(ortho.get('patr', '')),
                                neigh.get('ref_id', ''),
                                neigh.get('group', ''),
                                str(neigh.get('edges', '')),
                                str(neigh.get('patr', ''))
                            ]) + "\n")

        # Write neighbour FASTA
        neighbour_ids = set(entry['candidate_id'] for members in neighbour_members.values() for entry in members)
        neighbour_seqs = {seq_id: clean_members[seq_id] for seq_id in neighbour_ids if seq_id in clean_members}
        write_fasta(neighbour_seqs, neighbour_fasta_file)

        # construct final tree with all candidates with neighbours
        if not (os.path.isdir(tree_folder) and any(f.startswith("04") for f in os.listdir(tree_folder))):
            neigh_tree_output_folder = os.path.join(supplement_folder, f"{args.name}04_neighbouring_tree/")
            neigh_aln_candidate_file = f"{args.name}04_fin_alignment_candidates.fasta"
            neigh_aln_input_file = f"{args.name}04_fin_alignment_input.fasta"
            neigh_aln_file = f"{args.name}04_fin_alignment_input.fasta.aln"
            neigh_cln_aln_file = f"{args.name}04_fin_alignment_input.fasta.aln.cln"

            tree_file = parallel_tree_constructor(
                    len(clean_members), neigh_tree_output_folder, neigh_aln_candidate_file, neigh_aln_input_file,
                    neigh_aln_file, neigh_cln_aln_file, baits_path, neighbour_fasta_file, "04_neighbouring_", "",
                    args.mode_aln, args.mode_tree, args.mafft, args.muscle, args.raxml, args.fasttree,
                    args.cpumax, args.cpur, "N"
                )
            
            neigh_tree_path = replace_spaces_in_tree(tree_file[0])
            shutil.copy(neigh_tree_path, tree_folder)

        # --- 05 check for protein domains --- #
        domain_check_file_summary = os.path.join(txt_folder, f"{args.name}05_domain_check.txt")
        best_domain_match_file = os.path.join(txt_folder, f"{args.name}05_best_domain_match.txt")
        domain_check_hmm_result = os.path.join(supplement_folder, "05_domain_hmmsearch.txt")
        domain_check_hmm_output = os.path.join(supplement_folder, "05_domain_hmmsearch_output.txt")

        if not os.path.isfile(domain_check_file_summary):
            if hmm_domains_path is not None:
                domain_check_results, best_domain_match, domain_names = domain_check(
                    clean_members_file,
                    hmm_domains_path,
                    domain_check_hmm_result,
                    domain_check_hmm_output,
                    args.hmmsearch,
                    args.domain_Score
                )

                ref_family = load_ortholog_family(ref_mapping_file)
                with open(domain_check_file_summary, "w") as out1:
                    out1.write("RefMember\t" + "\t".join(domain_names) + "\n")
                    candidates = list(sorted(clean_members.keys()))
                    for candidate in candidates:
                        new_line = [subject_name_mapping_table[candidate]]
                        for dom in domain_names:
                            new_line.append("YES" if domain_check_results[candidate].get(dom, "") else "NO")
                        out1.write("\t".join(new_line) + "\n")

                with open(best_domain_match_file, "w") as out2:
                    out2.write("RefMember\tBestDomain\tOrtholog_Family\n")
                    for candidate in sorted(clean_members.keys()):
                        best_dom = best_domain_match.get(candidate, "-")
                        candidate_name = subject_name_mapping_table[candidate]
                        ortholog_family = ref_family.get(candidate_name, "None")
                        out2.write(f"{candidate_name}\t{best_dom}\t{ortholog_family}\n")

        # --- 06 check for common protein motifs --- #
        motif_check_file_summary = os.path.join(txt_folder, f"{args.name}06_motif_check.txt")
        motif_check_file_seqs = os.path.join(txt_folder, f"{args.name}06_motif_check_details.txt")
        motif_check_hmm_result = os.path.join(supplement_folder, "06_motif_hmmsearch.txt")
        motif_check_hmm_output = os.path.join(supplement_folder, "06_motif_hmmsearch_output.txt")
        
        if not os.path.isfile(motif_check_file_summary):
            if args.static_motifs.upper() in ['YES', 'Y']:
                if protein_motifs_path is not None:
                    motif_check_results, motif_names = static_motif_check(clean_members_file, protein_motifs_path)
                    flanklength = None
            else:
                if hmm_motifs_path is not None:
                    motif_check_results = motif_check(clean_members_file, hmm_motifs_path, motif_check_hmm_result, motif_check_hmm_output, args.hmmsearch, args.motif_cEvalue)
                elif protein_motifs_path is not None:
                    CYP_source = literature_baits_path if literature_baits_path else baits_path
                    flanklength = 3
                    hmm_motifs_path = create_motif_from_baits(protein_motifs_path, CYP_source, supplement_folder, mafft_available, muscle_available, flanklength)
                    motif_check_results = motif_check(clean_members_file, hmm_motifs_path, motif_check_hmm_result, motif_check_hmm_output, args.hmmsearch, args.motif_cEvalue)
                motif_names, motif_lengths = parse_motif_names_lengths(hmm_motifs_path)
            
            # # Prepare motif check data
            # motif_summary_lines = []
            # motif_details_lines = []
            
            # motif_summary_header = "RefMember\t" + "\t".join(motif_names) + "\n"
            # motif_details_header = "RefMember\t" + "\t".join(motif_names) + "\n"
            
            # motif_summary_lines.append(motif_summary_header)
            # motif_details_lines.append(motif_details_header)
            
            # candidates = list(sorted(clean_members.keys()))
            # for candidate in candidates:
            #     new_line_details = [subject_name_mapping_table[candidate]]
            #     new_line_summary = [subject_name_mapping_table[candidate]]
                
            #     for mot in motif_names:
            #         seq = motif_check_results[0][candidate].get(mot, "")
            #         if flanklength is not None:
            #             if seq:
            #                 trimmed_seq = seq[flanklength:-flanklength] if len(seq) > 2 * flanklength else seq
            #             else:
            #                 motif_len = motif_lengths.get(mot, 8)
            #                 trimmed_seq = "#" * motif_len
            #             new_line_details.append(trimmed_seq)
            #         else:
            #             new_line_details.append(seq)
            #         new_line_summary.append("YES" if seq else "NO")
                
            #     motif_summary_lines.append("\t".join(new_line_summary) + "\n")
            #     motif_details_lines.append("\t".join(new_line_details) + "\n")
            
            # # Write motif check files
            # write_output_file(
            #     file_path=motif_check_file_summary,
            #     content=motif_summary_lines,
            #     overwrite=True
            # )
            
            # write_output_file(
            #     file_path=motif_check_file_seqs,
            #     content=motif_details_lines,
            #     overwrite=True
            # )

            with open(motif_check_file_summary, "w") as out1:
                out1.write("RefMember\t" + "\t".join(motif_names) + "\n")
                with open(motif_check_file_seqs, "w") as out2:
                    out2.write("RefMember\t" + "\t".join(motif_names) + "\n")
                    candidates = list(sorted(clean_members.keys()))
                    for candidate in candidates:
                        new_line_details = [subject_name_mapping_table[candidate]]
                        new_line_summary = [subject_name_mapping_table[candidate]]
                        for mot in motif_names:
                            seq = motif_check_results[0][candidate].get(mot, "")
                            if flanklength is not None:
                                if seq:
                                    trimmed_seq = seq[flanklength:-flanklength] if len(seq) > 2 * flanklength else seq
                                else:
                                    motif_len = motif_lengths.get(mot, 8)
                                    trimmed_seq = "#" * motif_len
                                new_line_details.append(trimmed_seq)
                            else:
                                new_line_details.append(seq)
                            new_line_summary.append("YES" if seq else "NO")
                        out1.write("\t".join(new_line_summary) + "\n")
                        out2.write("\t".join(new_line_details) + "\n")
        
        # --- 07 filter expression matrix for candidates --- #
        filtered_expression_file = os.path.join(txt_folder, f"{args.name}07_filtered_expression_matrix.txt")

        if not os.path.isfile(filtered_expression_file):
            if expression_matrix_path is not None and os.path.isfile(expression_matrix_path):
                print(f"Processing expression matrix: {expression_matrix_path}")

                # Werte für Filterung (z.B. min_avg_tpm = 1.0, max_single_tpm = 5.0)
                filtered_lines, expressed, potential_pseudogenes = Candidate_Expression(
                    expression_matrix_path,
                    clean_members,
                    args.min_avg_tpm,
                    args.min_single_tpm
                )

                # Schreibe gefilterte Expressionsmatrix (nur "expressed")
                with open(filtered_expression_file, 'w') as out:
                    for line in filtered_lines:
                        out.write(line + '\n')

                print(f"{len(expressed)} candidates classified as expressed.")
                print(f"{len(potential_pseudogenes)} candidates classified as potential pseudogenes.")
            else:
                print("No expression matrix provided or file not found - skipping expression filtering step")
        else:
            print(f"Filtered expression matrix already exists: {filtered_expression_file}")

        if metadata_path is not None and os.path.isfile(metadata_path):
            print(f"Using metadata from: {metadata_path}")
            use_metadata(filtered_expression_file, metadata_path)
        else:
            print("No metadata file provided or file not found - skipping metadata integration")

        # --- 07.1 schreibe höchste Expression je Kandidat --- #
        highest_expression_file = write_highest_expression_info(
            filtered_expression_file,
            metadata_expression_file=filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_expression_matrix"),
            metadata_only_expression_file=filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")
        )

        # --- 08 find in species-specific paralogs (in-paralogs) --- #
        if args.collapse.upper() in ['YES', 'Y']:
            repr_clean_file = os.path.join(fasta_folder, f"{args.name}08_representative_paralogs.fasta")
            paralog_group_file = os.path.join(txt_folder, f"{args.name}08_representative_paralogs.txt")
            
            if not os.path.isfile(repr_clean_file) or not os.path.isfile(paralog_group_file):
                if expression_matrix_path is not None and metadata_path is not None:
                    metadata_only_expression_file = filtered_expression_file.replace("07_filtered_expression_matrix", "07_metadata_only_expression_matrix")
                    paralog_groups = establish_paralog_groups_with_expression_filter(final_tree_path, clean_members.keys(), args.paralogdist, 
                                                                                    metadata_only_expression_file, args.min_paralog_tpm)
                else:
                    paralog_groups = establish_paralog_groups(final_tree_path, clean_members.keys(), args.paralogdist)
                
                rep_per_group = get_represenative_paralog_per_group(paralog_groups, clean_members, repr_clean_file)
                
                # # Prepare paralog group data
                # paralog_lines = ["RepresentativeSeqID\tMembersOfParalogGroup\n"]
                # for rep_id, rep_seq in rep_per_group.items():
                #     for group in paralog_groups:
                #         if rep_id in group:
                #             group_members = "; ".join(group)
                #             paralog_lines.append(f"{rep_id}\t{group_members}\n")
                
                # # Write paralog group file
                # write_output_file(
                #     file_path=paralog_group_file,
                #     content=paralog_lines,
                #     overwrite=True
                # )

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

                # Parallel Tree Construction
                repr_tree_file = parallel_tree_constructor(
                    len(clean_members),
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
                    args.cpumax,
                    args.cpur,
                    "N"
                )

                # Final tree path = first tree file (if any)
                if repr_tree_file:
                    final_repr_tree_path = replace_spaces_in_tree(repr_tree_file[0])
                    shutil.copy(final_repr_tree_path, tree_folder)
                else:
                    print("Warning: No tree file was created during representative paralog processing.")

        # --- 09 create summary file --- #
        summary_file = os.path.join(txt_folder, f"{args.name}09_summary.txt")

        if not os.path.isfile(summary_file):
            print("Generating summary file...")

            summary_data = collect_summary_data(
                ref_mapping_file=ref_mapping_file,
                baits_info=baits_info,
                best_domain_match_file=best_domain_match_file,
                motif_check_file_summary=motif_check_file_summary,
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

            # # Prepare summary data
            # summary_lines = []
            # summary_lines.append("\t".join(header) + "\n")
            # summary_lines.append("\t".join(subheader) + "\n")
            # for seq_id in sorted(summary_data.keys()):
            #     summary_lines.append("\t".join(summary_data[seq_id]) + "\n")

            # # Write summary file
            # write_output_file(
            #     file_path=summary_file,
            #     content=summary_lines,
            #     overwrite=True
            # )

            with open(summary_file, "w") as out:
                out.write("\t".join(header) + "\n")
                out.write("\t".join(subheader) + "\n")
                for seq_id in sorted(summary_data.keys()):
                    out.write("\t".join(summary_data[seq_id]) + "\n")

            print(f"Summary file written to: {summary_file}")
        else:
            print(f"Summary file already exists: {summary_file}")

        # --- 10 Export --- #
        html_folder = os.path.join(result_folder, "HTML")
        html_link_file = os.path.join(html_folder, "CYP_output_links.html")

        if not os.path.isfile(html_link_file):
            print("Converting result .txt files to .html and creating output index...")
            convert_txt_to_html(txt_folder, html_folder, html_link_file)
            print(f"Output HTML overview written to: {html_link_file}")
        else:
            print(f"HTML overview already exists: {html_link_file}")

        # --- Tree plotting --- #
        ete3_available = importlib.util.find_spec("ete3") is not None
        tree_output_folder = os.path.join(result_folder, "TREES")
        if ete3_available and (not os.path.exists(tree_output_folder) or not any(f.endswith((".png", ".pdf")) for f in os.listdir(tree_output_folder))):
            print("Searching for .tre files and plotting trees...")
            plot_all_trees(tree_folder, tree_folder, baits_path, clean_members_file, format="png", layout="Polar", cand_color="green", bait_color="red")
        elif not ete3_available:
            print("ete3 is not installed. If you wish to convert the treefiles, please install ete3 using 'pip install ete3'")
        else:
            print(f"Tree visualizations already exist in: {tree_output_folder}")

        # --- Create Heatmaps --- #
        seaborn_available = importlib.util.find_spec("seaborn") is not None
        # heatmap_output_folder = os.path.join(result_folder, "HEATMAPS")
        # os.makedirs(heatmap_output_folder, exist_ok=True)

        # Prüfen, ob seaborn vorhanden ist UND ob bereits Heatmaps existieren
        heatmaps_exist = any(f.endswith((".png", ".pdf")) for f in os.listdir(heatmap_folder))

        if seaborn_available and not heatmaps_exist:
            print("Create Heatmaps from expression data...")
            generate_heatmaps_from_expression_data(txt_folder, heatmap_folder, format="png")
        elif not seaborn_available:
            print("seaborn is not available – skipping heatmap creation.")
        else:
            print(f"Heatmaps already exist in: {heatmap_folder}")

# endregion

if __name__ == "__main__":
    main()