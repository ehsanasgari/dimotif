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
from multiprocessing import Pool
import tqdm
from sklearn.feature_extraction.text import TfidfVectorizer
import numpy as np
from nltk import FreqDist
from proteinseq_util.biophysical import ProtSeqProp

class MotifProperties(object):
    def __init__(self):
        self.seq2freqstructs = FileUtility.load_obj('../data_config/seg2sec.pickle')
        # color dictionary for secondary structures
        color_dict={'e':'yellow', 'g':'blue', 'h':'blue', 'n':'red', 's':'red', 't':'red'}
        #H = alpha helix
        #B = residue in isolated beta-bridge
        #E = extended strand, participates in beta ladder
        #G = 3-helix (3/10 helix)
        #I = 5 helix (pi helix)
        #T = hydrogen bonded turn
        #S = bend
        #N = loop or other irregular structure
    def getMotifStructure(self, motif):
        return self.seq2freqstructs[motif] if motif in self.seq2freqstructs else None

    def get_motifs_pss_biophys(self, motifs):
        '''
            Extract the most frequent structures of protein motifs
            according to PDB sequences
            +
            The motif average bioPhysical properties
        '''

        ## Secondary structure
        motif2class=dict()
        for motif in tqdm.tqdm(motifs):
            if motif in self.seq2freqstructs:
                possiblities=self.seq2freqstructs[motif]
                poss=0
                pattern=''
                count=0
                while poss<0.5:
                    poss+=possiblities[count][1]
                    pattern+=''.join([possiblities[count][0]]*int(round(possiblities[count][1]*100)+1))
                    count+=1
                motif2class[motif]=FreqDist(list(pattern)).most_common(1)[0][0]

        ## Motif vector of properties, where the order in the vector is as follows
        #'mean_molecular_weight'
        #'mean_flexibility'
        #'instability'
        #'mean_surface_accessibility'
        #'mean_kd_hydrophobicity'
        #'mean_hydrophilicity'
        prop_vec=dict()
        for motif in tqdm.tqdm(motifs):
            PSP=ProtSeqProp(motif.upper())
            if not PSP.extended:
                prop_vec[motif]=[PSP.prop[x] for x in ['mean_molecular_weight','mean_flexibility','instability','mean_surface_accessibility','mean_kd_hydrophobicity','mean_hydrophilicity']]
        return motif2class, prop_vec

