import sys

sys.path.append('../')
from utility.file_utility import FileUtility
import pandas as pd
import tqdm
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from make_representations.cpe_apply import CPE
from multiprocessing import Pool
from chi2analysis.chi2analysis import Chi2Analysis
from utility.math_utility import get_sym_kl_rows

# download the file from http://hh-motif.ibdm.univ-mrs.fr/downloads/elm_annotation_complete_26.03.2016.tsv
df_elm = pd.read_table('../elm_eval/elm_annotation_complete_26.03.2016.tsv')

seq_dict = FileUtility.read_fasta_sequences_ids('../elm_eval/all_ELM.fasta')
SID = list(seq_dict.keys())
SID.sort()
elm_seqs = [seq_dict[seqID][0] for seqID in SID]
elm_lengths = np.array([len(x) for x in elm_seqs])

sequence_dict = dict()
fasta_files = FileUtility.recursive_glob('../../000_datasets/prot/ELM/', '*.fasta')

for fasta in tqdm.tqdm(fasta_files):
    name = fasta.split('/')[-1].split('.')[0]
    sequence_dict[name] = [seq.replace('-', '') for seq in FileUtility.read_fasta_sequences(fasta)]

ELM_instance_dict = dict()
ELM_instance_loc = dict()
ELM_motif = dict()
ELM_instances = dict()
incomplete = []
for x in df_elm.iterrows():
    if x[1].ELMType not in ELM_instance_dict:
        ELM_instance_dict[x[1].ELMType] = dict()
        ELM_instance_loc[x[1].ELMType] = dict()
        ELM_motif[x[1].ELMType] = dict()
        ELM_instances[x[1].ELMType] = dict()
    if x[1].ELMIdentifier not in ELM_instance_dict[x[1].ELMType]:
        ELM_instances[x[1].ELMType][x[1].ELMIdentifier] = []
        ELM_instance_dict[x[1].ELMType][x[1].ELMIdentifier] = []
        ELM_instance_loc[x[1].ELMType][x[1].ELMIdentifier] = []
        ELM_motif[x[1].ELMType][x[1].ELMIdentifier] = []
    try:
        ELM_instances[x[1].ELMType][x[1].ELMIdentifier] += sequence_dict[x[1].Primary_Acc]
    except:
        incomplete.append((x[1].Primary_Acc, x[1].ELMType, x[1].ELMIdentifier))
    ELM_instance_dict[x[1].ELMType][x[1].ELMIdentifier].append(x[1].Primary_Acc)
    ELM_instance_loc[x[1].ELMType][x[1].ELMIdentifier].append((x[1].Primary_Acc, x[1].Start - 1, x[1].End))
    ELM_motif[x[1].ELMType][x[1].ELMIdentifier].append(
        (x[1].Primary_Acc, seq_dict[x[1].Primary_Acc][0][(x[1].Start - 1):x[1].End], (x[1].Start - 1), x[1].End))

for k in ['DOC_USP7_MATH_1', 'DOC_AGCK_PIF_3', 'DOC_PIKK_1', 'DOC_MAPK_2', 'LIG_SH2_STAT5', 'LIG_14-3-3_2', 'LIG_Mtr4_Air2_1', 'LIG_TYR_ITIM', 'DEG_CRL4_CDT2_1', 'DEG_Kelch_KLHL3_1', 'DEG_APCC_KENBOX_2', 'DEG_SPOP_SBC_1', 'MOD_SUMO_rev_2', 'MOD_CDK_1', 'MOD_NEK2_1', 'MOD_SPalmitoyl_4', 'TRG_NES_CRM1_1', 'TRG_ER_KDEL_1', 'TRG_AP2beta_CARGO_1', 'TRG_LysEnd_APsAcLL_3']:
    typ=k[0:3]
    # metadata
    pos_seq = [x for x in FileUtility.read_fasta_sequences(
        '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ + '/' + k + '/pos.fasta')]
    neg_seq = [x for x in FileUtility.read_fasta_sequences(
        '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ + '/' + k + '/neg.fasta')]
    lables = [1] * len(pos_seq) + [0] * len(neg_seq)

    f = open('/mounts/data/proj/asgari/final_proj/dimotif/data_config/swissprot_ppe')
    res = dict()
    for size in [10000, 20000, 50000]:
        CPEO = CPE(f, separator=' ', merge_size=size)
        sequences = [(idx, x) for idx, x in enumerate(pos_seq + neg_seq)]
        pool = Pool(processes=20)
        for idx, x in tqdm.tqdm(pool.imap_unordered(CPEO.segment_with_keys, sequences, chunksize=20),
                                total=len(sequences)):
            if idx not in res:
                res[idx] = []
            res[idx].append(x)
        pool.close()
    res = [' '.join(res[i]) for i in range(len(res))]

    tp = TfidfVectorizer(use_idf=False, analyzer='word', ngram_range=(1, 1), stop_words=[], lowercase=True,
                         binary=False)
    X = tp.fit_transform(res)
    CH = Chi2Analysis(X, lables, tp.get_feature_names())
    a = CH.extract_features_fdr(
        '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ + '/' + k + '/motifs_sws.txt',
        direction=True)
    motifs = [l.split()[0] for l in FileUtility.load_list(
        '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ + '/' + k + '/motifs_sws.txt')[
                                    1::]]
    motifs = motifs[0:min(100, len(motifs))]
    idxs = [tp.get_feature_names().index(v) for v in motifs]
    pos_matrix = X.toarray()[0:len(pos_seq), idxs]
    DIST = get_sym_kl_rows(pos_matrix.T)
    FileUtility.save_obj(
        '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ + '/' + k + '/sym_KL_sws',
        [DIST, motifs])
