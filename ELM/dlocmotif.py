import sys

sys.path.append('../')
import subprocess
import tqdm
from multiprocessing import Pool
from utility.file_utility import FileUtility



def runit(typ):
    pos = '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ[0:3] + '/' + typ + '/pos.fasta'
    neg = '/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ[0:3] + '/' + typ + '/neg.fasta'
    pc=[x.replace('.','') for x in FileUtility.load_list(pos)]
    nc=[x.replace('.','') for x in FileUtility.load_list(neg)]
    FileUtility.save_list(pos,pc)
    FileUtility.save_list(neg,nc)
    #cmd1 = "java -jar /mounts/data/proj/asgari/final_proj/dimotif_eval/DLocalMotif/DLocalMotif.jar -p " + pos + " -n " + neg + " -w1 2 -w2 25 -delta 2 -o " +'/mounts/data/proj/asgari/final_proj/dimotif/ELM/new_eval/' + typ[0:3] + '/' + typ + '/out'
    #print(cmd1)
    #out=subprocess.getoutput(cmd1)
    return False


types = ['DOC_USP7_MATH_1',
         'DOC_AGCK_PIF_3',
         'DOC_PIKK_1',
         'DOC_MAPK_2',
         'LIG_SH2_STAT5',
         'LIG_14-3-3_2',
         'LIG_Mtr4_Air2_1',
         'LIG_TYR_ITIM',
         'DEG_SCF_TIR1_1',
         'DEG_Kelch_KLHL3_1',
         'DEG_APCC_KENBOX_2',
         'LIG_WRPW_1',
         'MOD_SUMO_rev_2',
         'MOD_CDK_1',
         'MOD_NEK2_1',
         'MOD_NMyristoyl',
         'TRG_NES_CRM1_1',
         'TRG_ER_KDEL_1',
         'TRG_AP2beta_CARGO_1',
         'TRG_LysEnd_APsAcLL_3']

pool = Pool(processes=20)
for x in tqdm.tqdm(pool.imap_unordered(runit, types, chunksize=1), total=len(types)):
    print(x)
pool.close()
