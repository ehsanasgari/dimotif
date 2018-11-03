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

def run_on_pair(inp):
    typ,k, seqs=inp
    FileUtility.ensure_dir('new_eval/'+typ+'/')
    FileUtility.ensure_dir('new_eval/'+typ+'/'+k+'/')
    train_cpe(seqs,'new_eval/'+typ+'/'+k+'/ppe.txt', 10000,'new_eval/'+typ+'/'+k+'/ppe_freq.txt', min_frequency=10)
    return True

# download the file from http://hh-motif.ibdm.univ-mrs.fr/downloads/elm_annotation_complete_26.03.2016.tsv
df_elm=pd.read_table('../elm_eval/elm_annotation_complete_26.03.2016.tsv')

seq_dict=FileUtility.read_fasta_sequences_ids('../elm_eval/all_ELM.fasta')
SID=list(seq_dict.keys())
SID.sort()
elm_seqs=[seq_dict[seqID][0] for seqID in SID]
elm_lengths=np.array([len(x) for x in elm_seqs])

sequence_dict=dict()
fasta_files=FileUtility.recursive_glob('../../000_datasets/prot/ELM/','*.fasta')

for fasta in tqdm.tqdm(fasta_files):
    name=fasta.split('/')[-1].split('.')[0]
    sequence_dict[name]=[seq.replace('-','') for seq in FileUtility.read_fasta_sequences(fasta)]

ELM_instance_dict=dict()
ELM_instance_loc=dict()
ELM_motif=dict()
ELM_instances=dict()
incomplete=[]
for x in df_elm.iterrows():
    if x[1].ELMType not in ELM_instance_dict:
        ELM_instance_dict[x[1].ELMType]=dict()
        ELM_instance_loc[x[1].ELMType]=dict()
        ELM_motif[x[1].ELMType]=dict()
        ELM_instances[x[1].ELMType]=dict()
    if x[1].ELMIdentifier not in ELM_instance_dict[x[1].ELMType]:
        ELM_instances[x[1].ELMType][x[1].ELMIdentifier]=[]
        ELM_instance_dict[x[1].ELMType][x[1].ELMIdentifier]=[]
        ELM_instance_loc[x[1].ELMType][x[1].ELMIdentifier]=[]
        ELM_motif[x[1].ELMType][x[1].ELMIdentifier]=[]
    try:
        ELM_instances[x[1].ELMType][x[1].ELMIdentifier]+=sequence_dict[x[1].Primary_Acc]
    except:
        incomplete.append((x[1].Primary_Acc,x[1].ELMType,x[1].ELMIdentifier ))
    ELM_instance_dict[x[1].ELMType][x[1].ELMIdentifier].append(x[1].Primary_Acc)
    ELM_instance_loc[x[1].ELMType][x[1].ELMIdentifier].append((x[1].Primary_Acc,x[1].Start-1,x[1].End))
    ELM_motif[x[1].ELMType][x[1].ELMIdentifier].append((x[1].Primary_Acc,seq_dict[x[1].Primary_Acc][0][(x[1].Start-1):x[1].End], (x[1].Start-1),x[1].End))

data_to_par=[]
for typ,dic in ELM_instance_dict.items():
    print (typ)
    for k,v in tqdm.tqdm(dic.items()):
        seqs=[x.lower() for x in ELM_instances[typ][k]]
        if len(seqs)>100:
            data_to_par.append((typ,k,seqs))

pool = Pool(processes=20)
res = []
for x in tqdm.tqdm(pool.imap_unordered(run_on_pair, data_to_par, chunksize=20), total=len(data_to_par)):
    res.append(x)
pool.close()
