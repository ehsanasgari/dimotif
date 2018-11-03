import sys
sys.path.append('../')
from utility.file_utility import FileUtility
import collections
import pandas as pd
import tqdm
import itertools
import numpy as np
from make_representations.cpe_efficient import train_cpe
from multiprocessing import Pool

#############################################################
# Simple script for learning segmentation steps from a fasta file
# Output: the file containing merging steps (i.e., "path_to_mergings"),
# can be used instead of Swiss-Prot merging steps
#############################################################

# Inputs
seq_dict=FileUtility.read_fasta_sequences_ids('sequences.fasta')
max_symbols=10000
min_freq_for_merging=10

# Output
path_to_mergings='ppe_mergings.txt'
path_to_merging_freqs='ppe_freq.txt'

#############################################################

SID=list(seq_dict.keys())
SID.sort()
seqs=[seq_dict[seqID][0] for seqID in SID]
train_cpe(seqs,path_to_mergings, max_symbols, path_to_merging_freqs, min_frequency=min_freq_for_merging)
