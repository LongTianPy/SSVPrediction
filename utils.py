### IMPORT
import os
import pandas as pd
from tqdm import tqdm
import sys
from pyfaidx import Fasta
import shutil
from sqlalchemy import create_engine
import pymysql
from maxentpy import maxent
from maxentpy import maxent_fast
from maxentpy.maxent import load_matrix5, load_matrix3
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)


### VARIABLES
matrix5 = load_matrix5()
matrix3 = load_matrix3()
genome_file = 'hg19.fa'

### FUNCTIONS
def get_mysql_sqlalchemy_engine_obj(conf_file='~/.my.cnf', **kwargs):
    engine =  create_engine('mysql+pymysql://',
                        creator=lambda: pymysql.connect(
                                            read_default_file=conf_file,
                                            charset='utf8',
                                            local_infile=True, db='var_SPss'
                                            )
                    )
    return engine.connect()

def process_line(line,fa):
    [chrom, pos, ref, alt] = line[:4]
    chrom = 'chr' + str(chrom).upper()
    pos = int(pos)
    return "|".join([chrom, str(pos), ref, alt])

def create_ss_pairs_by_window_size(line,fa,window_size,handler):
    chromposrefalt = process_line(line,fa)
    [chrom, pos, ref, alt] = chromposrefalt.split('|')
    chrom = str(chrom)
    pos = int(pos)
    pairs = []
    if ref == '-':
        ref = ''
    if alt == '-':
        alt = ''
    for i in range(window_size):
        ref_seq = fa[chrom][pos - 1 - i:pos - 1 + window_size - i].seq
        if len(ref) == len(alt):
            alt_seq = ref_seq
            alt_seq = alt_seq[:i] + alt + alt_seq[i + len(ref):]
        else:
            if len(ref) > len(alt):
                alt_seq = fa[chrom][pos - 1 - i:pos - 1 + window_size - i + len(ref) - len(alt)].seq
                alt_seq = alt_seq[:i] + alt + alt_seq[i + len(ref):]
            else:
                alt_seq = ref_seq
                alt_seq = alt_seq[:i] + alt + alt_seq[i + len(ref):]
                alt_seq = alt_seq[:window_size]
        # pairs.append([ref_seq, alt_seq])
        handler.write('{0}\t{1}\n'.format(chromposrefalt,'\t'.join([ref_seq,alt_seq])))

def create_potential_ss(inputfile):
    fa = Fasta(genome_file)
    out_5ss = open('5ss.txt', 'w')
    out_3ss_var_Spss = open('3ss_var_Spss.txt', 'w')
    out_3ss_maxentscan = open('3ss_MaxEntScan.txt', 'w')
    with open(inputfile, 'r') as f:
        line = f.readline()
        while line:
            line = f.readline().strip().split('\t')
            if len(line)<=1:
                break
            if len(line[2]) == 1 or len(line[3]) == 1:
                create_ss_pairs_by_window_size(line,fa,9,out_5ss)
                create_ss_pairs_by_window_size(line,fa,16,out_3ss_var_Spss)
                create_ss_pairs_by_window_size(line,fa,23,out_3ss_maxentscan)
    out_5ss.close()
    out_3ss_var_Spss.close()
    out_3ss_maxentscan.close()
    with open('5ss.txt','r') as f:
        lines = [i.strip().split('\t') for i in f.readlines()]
    with open('5ss_seq.txt','w') as f:
        pool = {}
        for line in lines:
            for i in line[1:]:
                if i not in pool:
                    pool[i] = 1
                    f.write('{0}\n'.format(i))
    with open('3ss_var_Spss.txt', 'r') as f:
        lines = [i.strip().split('\t') for i in f.readlines()]
    with open('3ss_var_Spss_seq.txt', 'w') as f:
        pool = {}
        for line in lines:
            for i in line[1:]:
                if i not in pool:
                    pool[i] = 1
                    f.write('{0}\n'.format(i))
    with open('3ss_MaxEntScan.txt', 'r') as f:
        lines = [i.strip().split('\t') for i in f.readlines()]
    with open('3ss_MaxEntScan_seq.txt', 'w') as f:
        pool = {}
        for line in lines:
            for i in line[1:]:
                if i not in pool:
                    pool[i] = 1
                    f.write('{0}\n'.format(i))
    shutil.copy('5ss_seq.txt', '/Users/longtian/Desktop/MaxEntScan')
    shutil.copy('3ss_MaxEntScan_seq.txt', '/Users/longtian/Desktop/MaxEntScan')

def run_maxentscan():
    cwd = os.getcwd()
    os.chdir('/Users/longtian/Desktop/MaxEntScan')
    os.system('perl score5.pl 5ss_seq.txt > score_5ss_MaxEntScan.txt')
    os.system('perl score3.pl 3ss_MaxEntScan_seq.txt > score_3ss_MaxEntScan.txt')
    os.chdir(cwd)
    shutil.copy('/Users/longtian/Desktop/MaxEntSca/score_5ss_MaxEntScan.txt', './')
    shutil.copy('/Users/longtian/Desktop/MaxEntSca/score_3ss_MaxEntScan.txt', './')

def get_delta_maxentscan():
    with open('5ss.txt','r') as f:
        pairs_5ss = [i.strip().split('\t') for i in f.readlines()]
    # pairs_5ss = pd.read_csv('5ss.txt',sep='\t',header=None,index_col=0)
    score_5ss = pd.read_csv('score_5ss_MaxEntScan.txt',sep='\t',header=None,index_col=0)
    score_5ss.columns = ['score']
    with open('delta_5ss_maxentscan.txt', 'w') as f:
        for i in pairs_5ss:
            idx = i[0]
            ref = i[1]
            alt = i[2]
            f.write('{0}\t{1}\n'.format(i,float(score_5ss.loc[alt,'score']) - float(score_5ss.loc[ref,'score'])))
    with open('3ss_MaxEntScan.txt','r') as f:
        pairs_3ss = [i.strip().split('\t') for i in f.readlines()]
    # pairs_3ss = pd.read_csv('3ss_MaxEntScan.txt', sep='\t', header=None, index_col=0)
    score_3ss = pd.read_csv('score_3ss_MaxEntScan.txt',sep='\t',header=None,index_col=0)
    score_3ss.columns = ['score']
    with open('delta_3ss_maxentscan.txt', 'w') as f:
        for i in pairs_3ss:
            idx = i[0]
            ref = i[1]
            alt = i[2]
            f.write('{0}\t{1}\n'.format(i,float(score_3ss.loc[alt,'score']) - float(score_3ss.loc[ref,'score'])))

def query_var_SPss():
    engine = get_mysql_sqlalchemy_engine_obj()
    conn = engine.connect()
    with open('5ss.txt','r') as f:
        pairs_5ss = [i.strip().split('\t') for i in f.readlines()]
    with open('5ss_seq.txt','r') as f:
        lines_5ss = [i.strip() for i in f.readlines()]
    lines_5ss = list(set(lines_5ss))
    df_5ss = pd.DataFrame(columns=['tendency-ratio'],index=lines_5ss)
    df_5ss.index.name='element'
    sql = "select element,`tendency-ratio` from score_9 where element in ({0})".format(','.join(["'"+i+"'" for i in lines_5ss]))
    from_db_5ss = pd.read_sql(sql,conn,index_col='element')
    merged_5ss = pd.merge(df_5ss,from_db_5ss,how='left',on='element')
    merged_5ss.fillna(-1,inplace=True)
    scores_5ss = []
    with open('delta_5ss_var_Spss.txt', 'w') as f:
        for i in pairs_5ss:
            if i[1][3:5]=='GT' and i[2][3:5]=='GT':
                delta = float(merged_5ss.loc[i[2],'tendency-ratio_y'])-float(merged_5ss.loc[i[1],'tendency-ratio_y'])
                f.write(('{0}\t{1}\n'.format(i[0], delta)))
    with open('3ss_var_Spss.txt', 'r') as f:
        pairs_3ss = [i.strip().split('\t') for i in f.readlines()]
    with open('3ss_var_Spss_seq.txt', 'r') as f:
        lines_3ss = [i.strip() for i in f.readlines()]
    lines_3ss = list(set(lines_3ss))
    df_3ss = pd.DataFrame(columns=['tendency-ratio'], index=lines_3ss)
    df_3ss.index.name = 'element'
    sql = "select element,`tendency-ratio` from score_16 where element in ({0})".format(
        ','.join(["'" + i + "'" for i in lines_3ss]))
    from_db_3ss = pd.read_sql(sql, conn, index_col='element')
    merged_3ss = pd.merge(df_3ss, from_db_3ss, how='left', on='element')
    merged_3ss.fillna(-1, inplace=True)
    scores_3ss = []
    with open('delta_3ss_varspss.txt', 'w') as f:
        for i in pairs_3ss:
            if i[1][11:13]=='AG' and i[2][11:13]=='AG':
                delta = float(merged_3ss.loc[i[2],'tendency-ratio_y'])-float(merged_3ss.loc[i[1],'tendency-ratio_y'])
                f.write(('{0}\t{1}\n'.format(i[0], delta)))

def main():
    inputfile = sys.argv[1]
    create_potential_ss(inputfile)
    # run_maxentscan()
    # get_delta_maxentscan()
    query_var_SPss()

if __name__ == '__main__':
    main()


