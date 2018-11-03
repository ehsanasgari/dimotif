
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
from multiprocessing import Pool
#import sentencepiece as spm
from sklearn.model_selection import train_test_split
from utility.math_utility import normalize_mat
import matplotlib
import matplotlib.pyplot as plt
from pylab import plot, show
import numpy as np
import numpy as np
import scipy.stats as st
import random

def sequences2segmented(sequences, vocab_sizes, nump):
    sequences=[(idx,x) for idx,x in enumerate(sequences)]
    segmented_seqs=dict()
    for vocab in tqdm.tqdm(vocab_sizes):
        f=open('../../protein_datasets/segmentations/swissprot_cpe','r')
        CPE_Applier=CPE(f,separator='', merge_size=vocab)
        pool = Pool(processes=nump)
        for idx,seg in pool.imap_unordered(CPE_Applier.segment_with_keys, sequences, chunksize=nump):
            if idx not in segmented_seqs:
                    segmented_seqs[idx]=[]
            segmented_seqs[idx].append(seg)
        pool.close()
    return [segmented_seqs[idx] for idx,x in enumerate(sequences)]


# read the whole swiss-prot
SWSSSEQ=FileUtility.read_fasta_sequences('swiss_prot.fasta')

# look at the changes for 1000 sequences with respect to the sampling sizes
randseq=random.sample(SWSSSEQ, 1000)
size_change=dict()
for vocab in tqdm.tqdm(np.arange(10000,1000000,10000)):
    size_change[vocab]=[]
    f=open('../data_config/swissprot_cpe','r')
    CPE_Applier=CPE(f,separator='', merge_size=int(vocab))
    for seq in randseq:
         size_change[vocab].append(len(CPE_Applier.segment(seq).split()))

all_samples=[]
for i in tqdm.tqdm(range(0,1000)):
    sample=[]
    for vocab in np.arange(10000,1000000,10000):
        sample.append(size_change[vocab][i])
    all_samples.append(-np.diff(sample))


sample_mat=np.mean(normalize_mat(all_samples),axis=0)
sample_mat_std=np.std(normalize_mat(all_samples),axis=0)

# check distributions
data = sample_mat
distributions = [st.laplace, st.norm, st.beta, st.alpha, st.expon, st.gamma]
mles = []

for distribution in distributions:
    pars = distribution.fit(data)
    mle = distribution.nnlf(pars, data)
    mles.append(mle)

results = [(distribution.name, mle) for distribution, mle in zip(distributions, mles)]
best_fit = sorted(zip(distributions, mles), key=lambda d: d[1])[0]
print ('Best fit reached using {}, MLE value: {}'.format(best_fit[0].name, best_fit[1]))



# plotting
fig, ax = plt.subplots(figsize=(9, 8))
ax.plot(np.arange(20000,1000000,10000)-10000,sample_mat)
ax.fill_between(np.arange(20000,1000000,10000)-10000, sample_mat-sample_mat_std,
                            sample_mat+sample_mat_std, alpha=0.2, linewidth=4,
                            linestyle='dashdot', antialiased=True)
plt.xticks(np.arange(20000,1000000,100000)-10000)
plt.xlim([10000,980000])
#ax.plot(np.arange(20000,1000000,10000),sample_mat+sample_mat_std,dashes=[1,2])
#ax.plot(np.arange(20000,1000000,10000),sample_mat-sample_mat_std,dashes=[1,2])
plt.xlabel('Merging steps')
plt.ylabel('Normalized average number of alternative segmentations for 1000 sequences')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'')
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
plt.rc('text', usetex=True)
plt.show()
