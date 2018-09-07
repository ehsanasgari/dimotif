#! /usr/bin/python

# -*- coding: utf-8 -*-
__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "DIMOTIF 2018"
__website__ = "llp.berkeley.edu/dimotif"


import sys
sys.path.append('../')
from utility.file_utility import FileUtility
from make_representations.cpe_apply import CPE
import tqdm
import numpy as np
from nltk import FreqDist

### Scrip to reproduce secondary structure mapping

def according_segmentation(segments,label):
    '''
    To segment sequences according to their secondary structures
    :param segments:
    :param label:
    :return:
    '''
    seg_idx=list(np.cumsum([len(x) for x in segments]))
    return [label[0 if idx==0 else seg_idx[idx-1]:pos] for idx, pos in  enumerate(seg_idx)]

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)



## script for segmentation of PDB secondary structures according to PPE units
## trained over Swiss Prot for different vocabulary sizes (here from "preferred numbers")
## but it can be also sampled from alpha distribution of SwissProt PPE lengths changes
sampled_lengths=[10000,20000,50000,100000,200000,500000,-1]
triples=dict()
for i in sampled_lengths:
    print(i)
    f=open('../data_config/swissprot_ppe','r')
    CPE_Applier=CPE(f,separator='', merge_size=i)
    sequences=FileUtility.read_fasta_sequences('../data_config/ss_N.txt')
    for pdb_idx, (x, y) in tqdm.tqdm(enumerate(pairwise(sequences))):
        segments=CPE_Applier.segment(x).split()
        label_segments=according_segmentation(segments,y)
        if i not in triples:
            triples[i]=[]
        triples[i]+=[(seg,label_segments[idx],pdb_idx) for idx,seg in enumerate(segments)]
for i in sampled_lengths:
    FileUtility.save_obj('../data_config/pdbsegments_'+str(i),triples[i])

## mapping of motifs to PDB ids
seq_ids=[x.strip() for x in FileUtility.load_list('../data_config/ss_N.txt') if x.strip()[0]=='>']
idx2pdb={idx:':'.join(val[1::].split(':')[0:2]) for idx, val in enumerate(seq_ids[::2])}
pdb2idx={':'.join(val[1::].split(':')[0:2]):idx for idx, val in enumerate(seq_ids[::2])}

structure_dict=dict()
for i in sampled_lengths:
    print(i)
    for seg,sec,idx in tqdm.tqdm(triples[i]):
        if seg not in structure_dict:
            structure_dict[seg]=dict()
        if sec not in structure_dict[seg]:
                structure_dict[seg][sec]=set()
        structure_dict[seg][sec].add(idx2pdb[idx])
FileUtility.save_obj('../data_config/motif2struct_pdbidx',structure_dict)


## mapping motifs to 10 most common secondary structures
import collections
structure_dict=collections.OrderedDict(structure_dict)
from nltk import  FreqDist
seg2structcount=dict()
for seg, struct_dict in tqdm.tqdm(structure_dict.items()):
    temp=[]
    for sec,pdbs in struct_dict.items():
        temp+=[sec]*len(pdbs)
    sumfreq=np.sum([y for x,y in FreqDist(temp).most_common(10)])
    seg2structcount[seg]=[(x,round(float(y/sumfreq),2)) for (x,y) in FreqDist(temp).most_common(10)]
FileUtility.save_obj('../data_config/seg2sec',seg2structcount)
