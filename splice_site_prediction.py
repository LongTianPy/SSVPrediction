### IMPORT
import pandas as pd
from tqdm import tqdm
import sys
from pyfaidx import Fasta
# from sqlalchemy import create_engine
# import pymysql
# from maxentpy import maxent
# from maxentpy import maxent_fast
# from maxentpy.maxent import load_matrix5, load_matrix3
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set(color_codes=True)

### VARIABLES
# matrix5 = load_matrix5()
# matrix3 = load_matrix3()
genome_file = 'hg19.fa'

### FUNCTIONS
def get_mysql_sqlalchemy_engine_obj(conf_file='~/.my.cnf', **kwargs):
    engine = create_engine('mysql+pymysql://',
                        creator=lambda: pymysql.connect(
                                            read_default_file=conf_file,
                                            charset='utf8',
                                            local_infile=True, db='var_SPss'
                                            )
                    )
    return engine

def create_ss_pairs_by_window_size(line,fa,window_size,handler):
    [chrom, pos, ref, alt] = line[:4]
    chrom = 'chr' + str(chrom).upper()
    pos = int(pos)
    pairs = []
    for i in range(window_size):
        ref_seq = fa[chrom][pos - 1 - i:pos - 1 + window_size - i].seq
        if len(ref) == len(alt):
            alt_seq = ref_seq
            alt_seq = alt_seq[:i] + alt + alt_seq[i + len(ref):]
        else:
            if len(ref) > len(alt):
                alt_seq = fa[chrom][pos - 1 - i:pos - 1 + window_size - i + len(ref) - len(alt)].seq
                alt_seq = alt_seq[:i] + alt + alt_seq[i  + len(ref):]
            else:
                alt_seq = ref_seq
                alt_seq = alt_seq[:i] + alt + alt_seq[i + 1:]
                alt_seq = alt_seq[:window_size]
        # pairs.append([ref_seq, alt_seq])
        handler.write('{0}\n'.format('\t'.join([ref_seq,alt_seq])))


def check_var_SPss_5ss(pairs,engine):
    idx = []
    scores = []
    for pair in pairs:
        for i in pair:
            idx.append(i)
    df = pd.DataFrame(columns=['tendency-ratio'],index=idx)
    df.index.name='element'
    sql = "select element,`tendency-ratio` from score_9 where element in ({0})".format(','.join(["'"+i+"'" for i in idx]))
    from_db = pd.read_sql(sql,engine,index_col='element')
    merged = pd.merge(df,from_db,how='left',on='element')
    merged.fillna(-1,inplace=True)
    for i in pairs:
        if i[0][3:5]!='GT' and i[1][3:5]!='GT':
            scores.append(float(merged.loc[i[1],'tendency-ratio_y'])-float(merged.loc[i[0],'tendency-ratio_y']))
        else:
            scores.append(0)
    return scores

def check_var_SPss_3ss(pairs,engine):
    idx = []
    scores = []
    for pair in pairs:
        for i in pair:
            idx.append(i)
    df = pd.DataFrame(columns=['tendency-ratio'],index=idx)
    df.index.name='element'
    sql = "select element,`tendency-ratio` from score_16 where element in ({0})".format(','.join(["'"+i+"'" for i in idx]))
    from_db = pd.read_sql(sql,engine,index_col='element')
    merged = pd.merge(df,from_db,how='left',on='element')
    merged.fillna(-1,inplace=True)
    for i in pairs:
        if i[0][11:13]!='AG' and i[1][11:13]!='AG':
            scores.append(float(merged.loc[i[1],'tendency-ratio_y'])-float(merged.loc[i[0],'tendency-ratio_y']))
        else:
            scores.append(0)
    return scores

def check_MaxEntScan_5ss(pairs):
    scores = []
    for i in pairs:
        score = [0,0]
        [ref, alt] = i
        score[0] = maxent_fast.score5(ref,matrix=matrix5)
        score[1] = maxent_fast.score5(alt,matrix=matrix5)
        scores.append(score[1]-score[0])
    return scores

def check_MaxEntScan_3ss(pairs):
    scores = []
    for i in pairs:
        score = [0,0]
        [ref, alt] = i
        score[0] = maxent_fast.score3(ref,matrix=matrix3)
        score[1] = maxent_fast.score3(alt,matrix=matrix3)
        scores.append(score[1]-score[0])
    return scores

def SSP(file):
    fa = Fasta(genome_file)
    # engine = get_mysql_sqlalchemy_engine_obj()
    # conn = engine.connect()
    delta_5ss_var_SPss = []
    delta_3ss_var_SPss = []
    delta_5ss_maxentscan = []
    delta_3ss_maxentscan = []
    pairs_5ss_var_SPss_all = []
    pairs_3ss_var_SPss_all = []
    pairs_3ss_for_maxentscan_all = []
    out_5ss = open('5ss.txt','w')
    out_3ss_var_Spss = open('3ss_var_Spss.txt','w')
    out_3ss_maxentscan = open('3ss_MaxEntScan.txt','w')
    with open(file, 'r') as f:
        line = f.readline()
        while line:
            line = f.readline().strip().split('\t')
            if len(line)<=1:
                break
            if len(line[2]) == 1 or len(line[3]) == 1:
                create_ss_pairs_by_window_size(line,fa,9,out_5ss)
                create_ss_pairs_by_window_size(line,fa,16,out_3ss_var_Spss)
                # pairs_5ss_var_SPss_all = pairs_5ss_var_SPss_all + pairs_5ss
                # pairs_3ss_var_SPss_all = pairs_3ss_var_SPss_all + pairs_3ss_for_var_SPss
                create_ss_pairs_by_window_size(line,fa,23,out_3ss_maxentscan)
    # scores_5ss_var_SPss = check_var_SPss_5ss(pairs_5ss_var_SPss_all,conn)
    # scores_3ss_var_SPss = check_var_SPss_3ss(pairs_3ss_var_SPss_all,conn)
    out_5ss.close()
    out_3ss_var_Spss.close()
    out_3ss_maxentscan.close()
    # scores_5ss_maxentscan = check_MaxEntScan_5ss(pairs_5ss_var_SPss_all)
    # scores_3ss_maxentscan = check_MaxEntScan_3ss(pairs_3ss_for_maxentscan_all)
    # delta_5ss_maxentscan = scores_5ss_maxentscan
    # delta_3ss_maxentscan = scores_3ss_maxentscan
    # # delta_5ss_var_SPss = scores_5ss_var_SPss
    # # delta_3ss_var_SPss = scores_3ss_var_SPss
    # return  delta_5ss_maxentscan,delta_3ss_maxentscan



if __name__ == '__main__':
    SSP('HGMD_current.tsv')
    # delta_5ss_maxentscan,delta_3ss_maxentscan = SSP('/Users/longtian/Desktop/HGMD_current.tsv')
    # sns.distplot(delta_5ss_var_SPss)
    # plt.savefig('Var_SPss_5ss.png')
    # plt.clf()
    # sns.distplot(delta_3ss_var_SPss)
    # plt.savefig('Var_SPss_3ss.png')
    # plt.clf()
    # sns.distplot(delta_5ss_maxentscan)
    # plt.savefig('MaxEntScan_5ss.png')
    # plt.clf()
    # sns.distplot(delta_3ss_maxentscan)
    # plt.savefig('MaxEntScan_3ss.png')
    # plt.clf()



