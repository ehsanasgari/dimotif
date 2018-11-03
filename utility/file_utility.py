__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - DeepSeq2Sec"
__website__ = "https://llp.berkeley.edu/DeepSeq2Sec/"


import sys
sys.path.append('../')

import _pickle as pickle
import codecs
import fnmatch
import os
from multiprocessing import Pool
import numpy as np
import tqdm
from Bio import SeqIO, Entrez
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy import sparse
import h5py
import shutil
from urllib.request import urlopen

class FileUtility(object):
    def __init__(self):
        print('File utility object created..')

    @staticmethod
    def download_fasta_files_from_NCBI(accession_ids, db="protein"):
        Entrez.email='asgari@berkeley.edu'
        seq_dict=dict()
        missing_seqs=[]
        print ('start downloading the sequences..')

        for seq_id in tqdm.tqdm(accession_ids):
            try:
                handle = Entrez.efetch(db=db, id=seq_id, rettype="fasta", retmode="text")
                seq=''.join(handle.read().split('\n')[1::])
                seq_dict[seq_id]=seq
            except:
                missing_seqs.append(seq_id)
        SID=list(seq_dict.keys())
        SID.sort()
        seqs=[seq_dict[seqID] for seqID in SID]
        FileUtility.create_fasta_file(output_file,seqs,SID, num=False, type=db)
        return missing_seqs
    @staticmethod
    def getUni(seq_id):
        try:
            handle = urlopen("http://www.uniprot.org/uniprot/"+seq_id+".xml")
            return seq_id, str(SeqIO.read(handle, "uniprot-xml").seq)
        except:
            return seq_id, False
    @staticmethod
    def download_fasta_files_from_Uniprot(output_file, accession_ids, num_p=10):
        seq_dict=dict()
        missing_seqs=[]
        pool = Pool(processes=num_p)
        print ('start downloading the sequences..')
        for seq_id, seq in tqdm.tqdm(pool.imap_unordered(FileUtility.getUni, accession_ids, chunksize=num_p),total=len(accession_ids)):
            if seq:
                seq_dict[seq_id]=seq
            else:
                missing_seqs.append(seq_id)
        pool.close()
        SID=list(seq_dict.keys())
        SID.sort()
        seqs=[seq_dict[seqID] for seqID in SID]
        FileUtility.create_fasta_file(output_file,seqs,SID, num=False, type='protein')
        return missing_seqs
    @staticmethod
    def getELM(seq_id):
        try:
            handle = urlopen("http://elm.eu.org/instances.fasta?q="+seq_id)
            return seq_id, str(handle.readlines()[1].strip().decode("utf-8"))
        except:
            return seq_id, False
    @staticmethod
    def download_fasta_files_from_ELM(output_file, accession_ids, num_p=10):
        seq_dict=dict()
        missing_seqs=[]
        pool = Pool(processes=num_p)
        print ('start downloading the sequences..')
        for seq_id, seq in tqdm.tqdm(pool.imap_unordered(FileUtility.getELM, accession_ids, chunksize=num_p),total=len(accession_ids)):
            if seq:
                seq_dict[seq_id]=seq
            else:
                missing_seqs.append(seq_id)
        pool.close()
        SID=list(seq_dict.keys())
        SID.sort()
        seqs=[seq_dict[seqID] for seqID in SID]
        FileUtility.create_fasta_file(output_file,seqs,SID, num=False, type='protein')
        return missing_seqs

    @staticmethod
    def create_fasta_file(file_address, corpus, label=None, num=True, type='protein'):
        seq_id_pairs = [('.'.join([str(idx + 1), label[idx] if label else ''] ) if num else label[idx], x) for idx, x in enumerate(corpus)]
        seq_recs = [SeqRecord(Seq(seq, generic_protein if type=='protein' else generic_dna), id=id, description='') for id, seq in seq_id_pairs]
        SeqIO.write(seq_recs, file_address, "fasta")



    @staticmethod
    def read_sequence_file(file_name_sample):
        '''
        :param file_name_sample:
        :return:
        '''
        corpus = []
        if file_name_sample[-1] == 'q':
            for cur_record in SeqIO.parse(file_name_sample, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name_sample, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return file_name_sample.split('/')[-1], corpus

    @staticmethod
    def read_sequence_file_length(file_name_sample):
        '''
        :param file_name_sample:
        :return:
        '''
        corpus = []
        if file_name_sample[-1] == 'q':
            for cur_record in SeqIO.parse(file_name_sample, "fastq"):
                corpus.append(str(cur_record.seq).lower())
        else:
            for cur_record in SeqIO.parse(file_name_sample, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return file_name_sample.split('/')[-1], len(corpus)


    @staticmethod
    def read_fasta_directory(file_directory, file_extenstion, only_files=[]):
        '''
        :param file_directory:
        :param file_extenstion:
        :param only_files:
        :return: list of fasta files, and a dic to map file to index
        '''
        if len(only_files) > 0:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion) if
                           x.split('/')[-1] in only_files]
        else:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion)]

        fasta_files.sort()
        mapping = {v: k for k, v in enumerate(fasta_files)}
        return fasta_files, mapping


    @staticmethod
    def save_obj(filename, value):
        with open(filename + '.pickle', 'wb') as f:
            pickle.dump(value, f)

    @staticmethod
    def load_obj(filename):
        return pickle.load(open(filename, "rb"))

    @staticmethod
    def ensure_dir(file_path):
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)

    @staticmethod
    def exists(file_path):
        return os.path.exists(file_path)

    @staticmethod
    def remove(file_path):
        os.remove(file_path)


    @staticmethod
    def remove_dir(file_path):
        shutil.rmtree(file_path)

    @staticmethod
    def save_list(filename, list_names):
        #FileUtility.ensure_dir(filename)
        f = codecs.open(filename, 'w', 'utf-8')
        for x in list_names:
            f.write(x + '\n')
        f.close()

    @staticmethod
    def load_list(filename):
        return [line.strip() for line in codecs.open(filename, 'r', 'utf-8').readlines()]

    @staticmethod
    def save_sparse_csr(filename, array):
        np.savez(filename, data=array.data, indices=array.indices,
                 indptr=array.indptr, shape=array.shape)

    @staticmethod
    def load_sparse_csr(filename):
        loader = np.load(filename)
        return sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])

    @staticmethod
    def _float_or_zero(value):
        try:
            return float(value)
        except:
            return 0.0

    @staticmethod
    def recursive_glob(treeroot, pattern):
        '''
        :param treeroot: the path to the directory
        :param pattern:  the pattern of files
        :return:
        '''
        results = []
        for base, dirs, files in os.walk(treeroot):
            good_files = fnmatch.filter(files, pattern)
            results.extend(os.path.join(base, f) for f in good_files)
        return results

    @staticmethod
    def read_fasta_sequences(file_name):
        corpus=[]
        for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus.append(str(cur_record.seq).lower())
        return corpus

    @staticmethod
    def read_fasta_sequences_ids(file_name):
        corpus=dict()
        for cur_record in SeqIO.parse(file_name, "fasta"):
                corpus[str(cur_record.id)]=(str(cur_record.seq).lower(),str(cur_record.description))
        return corpus


    @staticmethod
    def loadH5file(filename):
        f = h5py.File(filename, 'r')
        a_group_key = list(f.keys())[0]
        return list(f[a_group_key])
