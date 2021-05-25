import pandas as pd
import sys
from os.path import exists, join,expanduser,basename
from tqdm import tqdm
from api_tools.itol_func import *
from glob import glob
from ete3 import Tree
import plotly.express as px
from ete3 import NCBITaxa
from global_search.classification_script import _classificated
from subprocess import check_call
import os
from collections import defaultdict
from Bio import SeqIO
from collections import Counter

def parse_blast(indir=None,outdir=None,infiles=[],
                gene_lists=None,dbfaa=None):
    gene_lists = open(gene_lists).read().split('\n')
    locus2gene = defaultdict(list)
    for row in gene_lists:
        if not row:
            continue
        rows = row.split('\t')
        locus2gene[rows[0]] = rows[1]
        
    def get_genome2genes(sub_df,locus2gene,locus_any):
        _sub_df = sub_df.loc[sub_df[1].isin(locus2gene)]
        l2g = {}
        for l,idx_list in _sub_df.groupby(0).groups.items():
            genes = [locus2gene[_] for _ in sub_df.loc[idx_list,1]]
            gene = sorted(genes,key=lambda x:Counter(genes)[x])[-1]
        # genes = [locus2gene[_] for _ in _sub_df[1]]
            l2g[l] = gene
        # l = _sub_df[0]
        # l2g = dict(zip(l,genes))
        # locus2gene_sub.update()
        # sub_df.loc[:,'gene'] = [locus2gene.get(_,'') for _ in sub_df[1]]
        # genes = set(genes)
        genes = list(set(l2g.values()))
        return {"GCA_"+locus_any.split('_')[0].replace('v','.'): genes},l2g

    # infile = "/mnt/home-backup/thliao/AOB/whole_tree_20200412/rebuild_gene_tree/nxrA/prepared.faa"
    curated_gene2length = {_.id: len(str(_.seq))
                        for _ in SeqIO.parse(dbfaa, format='fasta')}
    genome2genes = {}
    locus2gene_focal_phyla = {}
    if not infiles:
        infiles = glob(join(indir, '*.tab'))
        
    for tab in tqdm(infiles):
        if os.path.getsize(tab) == 0:
            continue
        df = pd.read_csv(tab,sep='\t',header=None)
        # df = df.loc[~df[1].isin(removed_false_positive_ref_faa),:]
        if df.shape[0]==0:
            continue
        df.loc[:,'total_length'] = [curated_gene2length[_] for _ in df[1]]
        sub_df = df.loc[(df[10]<=1e-20) & (df[2] >=60) & (abs(df[7]-df[6])/df['total_length'] >0.65) ]
        sub_df = sub_df.groupby(0).head(10)
        g = df.iloc[0,0]
        _dict,l2g_d = get_genome2genes(sub_df,locus2gene,g)
        genome2genes.update(_dict)
        locus2gene_focal_phyla.update(l2g_d)
        
    # write out the results incase repeat run above codes
    with open(join(outdir,'genome2genes.txt'),'w') as f1:
        for genome,genes in genome2genes.items():
            f1.write('\t'.join([genome] + list(genes)) + '\n')
    with open(join(outdir,'locus2genes.list'),'w') as f1:
        for locus,gene in locus2gene_focal_phyla.items():
            f1.write('\t'.join([locus,gene])+ '\n')
    
    info2style = {'hzsA':{'color':'#ff0000',
                      'info':'hzsA',
                      "shape":"1"},
              'hzsB':{'color':'#ff0000',
                      'info':'hzsB',
                      "shape":"1"},
              'hzsC':{'color':'#ff0000',
                      'info':'hzsC',
                      "shape":"1"},
              'hao':{'color':'#b68100',
                      'info':'hao',
                      "shape":"1"},
              'hdh':{'color':'#b68100',
                      'info':'hdh',
                      "shape":"1"},
            'hzo':{'color':'#b68100',
                      'info':'hzo',
                      "shape":"1"},
            "potential_hzs":{'color':'#616161',
                      'info':'potential_hzs',
                      "shape":"3"}
              }
    text = to_binary_shape({k:[v] for k,v in locus2gene_focal_phyla.items()},
                        info2style=info2style,
                        unfilled_other=True)
    with open(join(outdir,'locus2gene.txt'),'w') as f1:
        f1.write(text)
    

    # write out annotation
    _info2style = {'hzsA':{'color':'#ff0000',
                      'info':'hzsA',
                      "shape":"1"},
              'hzsB':{'color':'#ff0000',
                      'info':'hzsB',
                      "shape":"1"},
              'hzsC':{'color':'#ff0000',
                      'info':'hzsC',
                      "shape":"1"},
              'hao':{'color':'#b68100',
                      'info':'hao',
                      "shape":"1"},
              'hdh':{'color':'#b68100',
                      'info':'hdh',
                      "shape":"1"},
    }
    
    sub_genome2genes = {k:[_ for _ in v if _ in _info2style] for k,v in genome2genes.items() }
    #text = to_color_range(_g2g, info2color=info2style)
    text = to_binary_shape(sub_genome2genes,info2style=_info2style,unfilled_other=True,
                           manual_v=['hzsA','hzsB','hzsC','hao','hdh'])
    with open(join(outdir,'target_genes_binary.txt'),'w') as f1:
        f1.write(text)


    all_locus = {_:_ for _ in locus2gene_focal_phyla}
    text = to_label(all_locus)
    with open(join(outdir,'reset_names.txt'),'w') as f1:
        f1.write(text)
            

    return genome2genes,locus2gene_focal_phyla
