import pandas as pd
import sys
from os.path import exists, join,expanduser,basename
from tqdm import tqdm
from api_tools.itol_func import *
from glob import glob
from ete3 import Tree
import plotly.express as px
from ete3 import NCBITaxa
from subprocess import check_call
import os
from collections import defaultdict
from Bio import SeqIO
from collections import Counter
from .tools import parse_blast
ncbi = NCBITaxa()

def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    return genome_ID.split('_')[-1].replace('.','v')

def convert_genome_ID_rev(locus_ID):
    # for 900078535v2
    # it will return
    if '_' in locus_ID:
        locus_ID = locus_ID.split('_')[0]
    return 'GCA_'+locus_ID.replace('v','.')

blastp_pth = '/home-user/software/blast/latest/bin/blastp'
def run(cmd):
    check_call(cmd,shell=True)

def reformat(s):
    a = s.split('_')[-1]
    if not '_' in s:
        return s
    try:
        float(a)
        return s
    except:
        if len(s.rpartition('_')[-1]) == 1:
            return s
        else:
            return s.rpartition('_')[0]


    
tmp_dir = join('./gene_annotation/')
os.makedirs(tmp_dir,exist_ok=True)

collected_gs = expanduser('~/project/nitrogen_cycle/curated_genes/')
genes = ['hzsA','hzsB','hzsC','hao','hdh','hzo',]
do_genes = glob(join(collected_gs,'*.faa'))
do_genes = [_ for _ in do_genes 
            if basename(_).replace('.faa','') in genes]
all_records = {r.id:basename(_).replace('.faa','') 
               for _ in do_genes 
               for r in SeqIO.parse(_,'fasta')}
with open(join(tmp_dir,'db_curated.list'),'w') as f1:
    for k,v in all_records.items():
        f1.write(f'{k}\t{v}\n')
with open(join(tmp_dir,'db.faa'),'w') as f1:
    records = [r
               for _ in do_genes 
               for r in SeqIO.parse(_,'fasta')]
    SeqIO.write(records,f1,'fasta-2line')

all_ids = open('./over20p_genomes_add_cyano.list').read().strip('\n').split('\n')
p_faa_dir = '/mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files'
redo = False

genome2collect_genes = defaultdict(list)
gene2locus = defaultdict(set)

#gene_name = basename(fa).replace('.faa','').strip()
for aid in tqdm(all_ids):
    db_faa = join(p_faa_dir,f"{aid}.faa")
    otab = join(tmp_dir,'total',basename(db_faa).replace('.faa','').strip()+'.tab')
    if not exists(dirname(otab)):
        os.makedirs(dirname(otab))
    if (not exists(otab)) or redo:
        cmd = f"{blastp_pth} -db {tmp_dir}/db.faa -query {db_faa} -out {otab} -num_threads 20 -outfmt 6 -max_target_seqs 100000 "
        run(cmd)
    # if getsize(otab)!=0:
    #     parse_blast()
    #     df = pd.DataFrame(otab)
    # for row in open(otab):
    #     if float(row.split('\t')[-2]) < 1e-20 and float(row.split('\t')[2]) > 60:
    #         genome2collect_genes[aid].append(gene_name)
    #         gene2locus[gene_name].add(row.split("\t")[0])
# genome2collect_genes = {k:set(v) for k,v in genome2collect_genes.items()}
parse_blast(indir=join(tmp_dir,'total'),
            outdir=tmp_dir,
            gene_lists=join(tmp_dir,'db_curated.list'),
            dbfaa=join(tmp_dir,'db.faa'))


