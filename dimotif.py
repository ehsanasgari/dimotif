import argparse
import os
import os.path
import sys
sys.path.append('../')
from make_representations.cpe_apply import CPE
from utility.file_utility import FileUtility
from multiprocessing import Pool
import tqdm
from sklearn.feature_extraction.text import TfidfVectorizer
import numpy as np
from nltk import FreqDist
from proteinseq_util.biophysical import ProtSeqProp
from utility.math_utility import normalize_mat
import scipy.stats as st
from chi2analysis.chi2analysis import Chi2Analysis
from utility.math_utility import get_sym_kl_rows
from clustering.hierarchical import HierarchicalClutering
#from proteinseq_util.motif_tree_visualization import VisualizeTreeOfMotifs

class DiMotif(object):
    def __init__(self, pos_fasta, neg_fasta, output_path, segmentation_schemes=10, topN=100):
        '''

        '''
        if not isinstance(pos_fasta, str):
            self.pos=pos_fasta
        elif pos_fasta.split('.')[-1]=='txt':
            self.pos=FileUtility.load_list(pos_fasta)[0:10]
        elif pos_fasta.split('.')[-1]=='fasta':
            self.pos=FileUtility.read_fasta_sequences(pos_fasta)
        if not isinstance(neg_fasta, str):
            self.neg=neg_fasta
        elif neg_fasta.split('.')[-1]=='txt':
            self.neg=FileUtility.load_list(neg_fasta)[0:10]
        elif neg_fasta.split('.')[-1]=='fasta':
            self.neg=FileUtility.read_fasta_sequences(neg_fasta)
        self.seqs=[seq.lower() for seq in self.pos+self.neg]
        self.labels=[1]*len(self.pos)+[0]*len(self.neg)
        self.segmentation_schemes=segmentation_schemes
        self.load_alpha_distribution()
        self.prepare_segmentations()
        print (output_path)
        FileUtility.ensure_dir(output_path)
        self.output_path=output_path
        self.motif_extraction(topN)

    def load_alpha_distribution(self):
        swiss_size_change=FileUtility.load_obj('data_config/swiss_1000_samples.pickle')
        all_samples=[]
        for i in tqdm.tqdm(range(0,1000)):
            sample=[]
            for vocab in np.arange(10000,1000000,10000):
                sample.append(swiss_size_change[vocab][i])
            all_samples.append(-np.diff(sample))

        sample_mat=np.mean(normalize_mat(all_samples),axis=0)
        sample_mat_std=np.std(normalize_mat(all_samples),axis=0)
        self.alpha_param = st.alpha.fit(sample_mat)

    def get_alpha_samples(self):
        r = st.alpha.rvs(self.alpha_param[0], size=self.segmentation_schemes)
        idx=np.array(np.round(10000+(r*10000)),dtype=np.int32).tolist()
        idx.sort()
        return idx

    def prepare_segmentations(self):
        segmented_seqs=[]
        vocab_sizes=self.get_alpha_samples()
        for i, vocab in tqdm.tqdm(enumerate(vocab_sizes)):
            f=open('data_config/swissprot_ppe','r')
            CPE_Applier=CPE(f,separator='', merge_size=vocab)
            for idx, seq in enumerate(self.seqs):
                if i ==0:
                    segmented_seqs.append([CPE_Applier.segment(seq)])
                else:
                    segmented_seqs[idx]+=[CPE_Applier.segment(seq)]
        self.extended_sequences=[' '.join(l) for l in segmented_seqs]
        self.possible_segmentations=['@@@'.join(l) for l in segmented_seqs]


    def motif_extraction(self, topn=100):
        cpe_vectorizer = TfidfVectorizer(use_idf=False, analyzer='word',
                                              norm=None, stop_words=[], lowercase=True, binary=False, tokenizer=str.split)
        tf_vec=cpe_vectorizer.fit_transform(self.extended_sequences)
        vocab=cpe_vectorizer.get_feature_names()
        CH=Chi2Analysis(tf_vec,self.labels,vocab)
        vocab_binary=[x[0] for x in CH.extract_features_fdr(self.output_path+'/motifs.txt', N=topn, alpha=5e-2, direction=True, allow_subseq=True, binarization=True, remove_redundant_markers=False) if x[1]>0]
        vocab_binary=vocab_binary[0:min(100,len(vocab_binary))]
        idxs=[vocab.index(v) for v in vocab_binary]
        pos_matrix=tf_vec.toarray()[0:len(self.pos),idxs]
        DIST=get_sym_kl_rows(pos_matrix.T)
        FileUtility.save_obj(self.output_path+'/sym_KL', DIST)
        #HC=HierarchicalClutering(DIST,vocab_binary)
        self.motifs=vocab_binary
        #self.tree=HC.nwk
        #FileUtility.save_list(self.output_path+'/motif_tree.txt', [HC.nwk])

def checkArgs(args):
    '''
        This function checks the input arguments and returns the errors (if exist) otherwise reads the parameters
    '''
    # keep all errors
    err = "";
    # Using the argument parser in case of -h or wrong usage the correct argument usage
    # will be prompted
    parser = argparse.ArgumentParser()

    def file_choices(choices,fname):
        ext = os.path.splitext(fname)[1][1:]
        if ext not in choices:
           parser.error("file doesn't end with one of {}".format(choices))
        return fname

    ## to do : chi2 print

    # positive file #################################################################################################
    parser.add_argument('--pos', action='store', dest='pos_file', type=lambda s:file_choices(("txt","fasta"),s),
                        help='positive fasta or txt sequence file')

    # negative file #######################################################################################################
    parser.add_argument('--neg', action='store', dest='neg_file', type=lambda s:file_choices(("txt","fasta"),s),
                        help='negative fasta or txt sequence file')

    # output directory #################################################################################################
    parser.add_argument('--outdir', action='store', dest='output_dir', default=False, type=str,
                        help="directory for storing the output files, if doesn't exist will be created.")

    # to override the previous files or to continue ####################################################################
    parser.add_argument('--topn', action='store', dest='topn',default=100, type=int,
                        help='How many motifs to extract if possible?')

    # to override the previous files or to continue ####################################################################
    parser.add_argument('--segs', action='store', dest='segs',default=10, type=int,
                        help='How many segmentation samples for each seq')

    parsedArgs = parser.parse_args()

    if (not os.access(parsedArgs.pos_file, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the positive file!"
        return err
    if (not os.access(parsedArgs.neg_file, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the negative file!"
        return err
    try:
        print('Extract motifs..')
        DMF=DiMotif(parsedArgs.pos_file,parsedArgs.neg_file,parsedArgs.output_dir, topN=parsedArgs.topn, segmentation_schemes=parsedArgs.segs)
        print('Visualize motifs..')
        #VisualizeTreeOfMotifs(DMF.tree, DMF.motifs)
    except:
        print ('error occured')

if __name__ == '__main__':
    err = checkArgs(sys.argv)
    if err:
        print(err)
        exit()

