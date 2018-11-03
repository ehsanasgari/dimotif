import sys
sys.path.append('../')
from utility.file_utility import FileUtility
import collections
import pandas as pd
import tqdm
import itertools
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from make_representations.cpe_apply import CPE
from multiprocessing import Pool
from chi2analysis.chi2analysis import Chi2Analysis
from utility.math_utility import get_sym_kl_rows

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


for typ,dic in ELM_instance_dict.items():
    print(typ)
    for k,v in dic.items():
        if not False:#FileUtility.exists('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/pos.txt'):
            try:
                # get positive sequences
                pos_seq=list(set([x.lower() for x in ELM_instances[typ][k]]))
                mean_l_pos=np.mean([len(x) for x in pos_seq])
                std_l_pos=np.std([len(x) for x in pos_seq])

                # get negative sequences
                selected_rand_seq=FileUtility.load_list('../elm_eval/rand_neg_samples.txt')
                length_rand_seq = np.array([len(seq) for seq in selected_rand_seq])
                ## length
                selected_rand_seq=[selected_rand_seq[idx] for idx,l in enumerate(length_rand_seq) if mean_l_pos-std_l_pos<l and l<mean_l_pos+std_l_pos]
                ## dissimilarity
                tp=TfidfVectorizer(use_idf=False, analyzer='char', ngram_range=(3,3), norm='l2', stop_words=[], lowercase=True, binary=False)
                corpus=pos_seq+selected_rand_seq
                X=tp.fit_transform(corpus).toarray()
                sim=X.dot(X.T)
                distances=1-sim
                random_idx_neg=np.random.choice(list(range(0,len(selected_rand_seq))),len(pos_seq),(np.sum(distances[0:len(pos_seq),len(pos_seq):],axis=0)/np.sum(distances[0:len(pos_seq),len(pos_seq):])).tolist())
                neg_seq=[selected_rand_seq[idx] for idx in random_idx_neg]

                # metadata
                lables=[1]*len(pos_seq)+[0]*len(neg_seq)
                FileUtility.create_fasta_file('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/pos.fasta',pos_seq)
                FileUtility.create_fasta_file('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/neg.fasta',neg_seq)
                # f=open('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/ppe.txt')
                # res = dict()
                # for size in [1000,2000,5000]:
                #     CPEO=CPE(f, separator=' ',merge_size=size)
                #     sequences=[(idx,x) for idx,x in enumerate(pos_seq+neg_seq)]
                #     pool = Pool(processes=20)
                #     for idx,x in tqdm.tqdm(pool.imap_unordered(CPEO.segment_with_keys, sequences, chunksize=20), total=len(sequences)):
                #         if idx not in res:
                #             res[idx]=[]
                #         res[idx].append(x)
                #     pool.close()
                # res=[' '.join(res[i]) for i in range(len(res))]
                #
                # tp=TfidfVectorizer(use_idf=False, analyzer='word', ngram_range=(1,1), stop_words=[], lowercase=True, binary=False)
                # X=tp.fit_transform(res)
                # CH=Chi2Analysis(X,lables,tp.get_feature_names())
                # a=CH.extract_features_fdr('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/motifs.txt',direction=True)
                # motifs=[l.split()[0] for l in FileUtility.load_list('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/motifs.txt')[1::]]
                # motifs=motifs[0:min(100,len(motifs))]
                # idxs=[tp.get_feature_names().index(v) for v in motifs]
                # pos_matrix=X.toarray()[0:len(pos_seq),idxs]
                # DIST=get_sym_kl_rows(pos_matrix.T)
                # FileUtility.save_obj('/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/'+typ+'/'+k+'/sym_KL', [DIST,motifs])
            except:
                print (k,typ)
