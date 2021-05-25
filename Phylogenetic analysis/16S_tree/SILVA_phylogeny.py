
# download all 20142 rrna sequence from silva database affiliated with Brocadiae taxonomic level

import plotly.express as px
# from api_tools.metadata_for.auto_classify import *
import pandas as pd
from for_software.for_cd_hit.parse import parse_clstr
from api_tools.itol_func import *
from collections import defaultdict
from Bio import SeqIO
from subprocess import check_call
from os.path import *
from tqdm import tqdm
import os

os.chdir("/mnt/home-backup/thliao/plancto/rrna/silva")

ref_seq = './Brocadiae_20142.fasta'

# 1 run cd-hit to cluster and output representative sequence to following analysis
threshold = 0.97
n = int(threshold*100)
cmd = f"cd-hit -c {threshold} -i ./Brocadiae_20142.fasta -o ./Brocadiae_{n}.fasta -d 0"
cluster_log = f"./Brocadiae_{n}.fasta.clstr"
clustered_fasta = f"./Brocadiae_{threshold*100}.fasta"
if not exists(cluster_log):
    check_call(cmd, shell=1)


# 2. add ref
need_to_add = "../ref_seq/16S.fasta"
new_seq = list(SeqIO.parse(clustered_fasta, format='fasta'))
need_to_add = list(SeqIO.parse(need_to_add, format='fasta'))

name = "prepared"
with open(f'./{name}.fasta', 'w') as f1:
    for r in new_seq:
        f1.write(f">{r.id}\n{r.seq}\n")
    for r in need_to_add:
        f1.write(
            f">GCA_{r.id.split('_')[0].replace('v','.')}\n{r.seq.transcribe()}\n")

cmd = f"mafft ./{name}.fasta > ./{name}.aln; "
check_call(cmd, shell=1)

cmd = "iqtree -nt AUTO -wbtl -bb 1000 -m MFP -redo -mset GTR,HKY85,K80 -s ./prepared.aln -pre ./prepared_iqtree"
check_call(cmd, shell=1)
cmd = "python3 ~/bin/format_newick.py rename -i ./prepared_iqtree.contree -o ./prepared_iqtree.newick "
check_call(cmd, shell=1)
cmd = "python3 ~/bin/format_newick.py reroot -i ./prepared_iqtree.newick -o ./prepared_rerooted.newick -r GCA_000019965.1,GCA_900097105.1 -f 3"
check_call(cmd, shell=1)
cmd = "python3 ~/bin/format_newick.py sort -i ./prepared_rerooted.newick -o ./prepared_sorted.newick -f 3"
check_call(cmd, shell=1)
text = to_node_symbol("./prepared_iqtree.newick")
with open('./prepared_bp.txt', 'w') as f1:
    f1.write(text)

# 3.1 annotate outgroup of anammox
outgroup = ['GCA_001303885.1', 'GCA_003551305.1',
            'GCA_003551565.1', 'GCA_003567495.1']
text = to_binary_shape({k: ['outgroup'] for k in outgroup},
                       {"outgroup": {'color': '#33691E'}})
with open('./anammox_outgroup_annotate.txt', 'w') as f1:
    f1.write(text)

# 3.2 annotate number of clustered sequences
cluster2seqs, cluster2repr = parse_clstr(cluster_log)
id2val = {}
for c, repr_seq in cluster2repr.items():
    num_seqs = cluster2seqs[c]
    id2val[repr_seq] = len(num_seqs)
text = to_simple_bar(id2val)
with open('./rrna_silva_cluster_size_bar.txt', 'w') as f1:
    f1.write(text)

text = to_simple_bar(id2val)


# 4. read metadata of SILVA
records = SeqIO.parse(ref_seq, format='fasta')
acc_list = set([_.id.split('.')[0] for _ in records])
# download from https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/full_metadata/SILVA_138_SSUParc.full_metadata.gz
all_m = open(
    '/mnt/home-backup/thliao/silva/SILVA_138/SILVA_138_SSUParc.full_metadata')
seq2info = defaultdict(dict)
header = all_m.readline().strip('\n').split('\t')
for row in tqdm(all_m):
    rows = row.strip('\n').split('\t')
    if rows[0] in acc_list:
        seq2info[rows[0]] = dict(zip(header[1:],
                                     rows[1:]))
df = pd.DataFrame.from_dict(seq2info, orient='index')
df.to_csv('./Brocadiae_20142.metadata', sep='\t', index=1, index_label="acc")


# manual assigns habitats

genome2habitat = pd.read_excel("/mnt/home-backup/thliao/plancto/rawdata/classified_metadata curated.xlsx", index_col=0)
genome2habitat = dict(zip(genome2habitat.index,
                          genome2habitat.loc[:, "habitat old?"]))
manual_assigns = pd.read_excel(
    "./silva habitats manual assigns.xlsx", index_col=None)
isolation2types = dict(zip(manual_assigns['isolation'],
                           manual_assigns["manual assign"]))
tab = pd.read_csv('./Brocadiae_20142.metadata', sep='\t', index_col=0)


refined_habitats = {'sediment': 'sediment',
                    'marine sediment': 'marine',
                    'others': 'unknown',
                    'unknown': 'unknown',
                    'terrestrial': 'terrestrial',
                    'marine': 'marine',
                    'human-made': 'human-made',
                    'groundwater': 'marine',
                    'freshwater': 'freshwater',
                    'extreme': 'marine'}
sediment = [k for k,v in isolation2types.items() if v == 'sediment']
sediment_refined = [k for k in sediment if 'estuary' in k or 'saline' in k or 'estuarine' in k or 'intertidal' in k or 'mangrove' in k or 'salt mash' in k]
for k in sediment_refined:
    isolation2types[k] = 'marine'

tab.loc[:, 'Ecosystem type'] = [refined_habitats[isolation2types.get(str(_).lower(), 'unknown')]
                                for _ in tab['isolation_source']]


seq2eco = dict(zip(tab.index,
                   tab['Ecosystem type']))

cluster2seqs, cluster2repr = parse_clstr(cluster_log)

seq_id2habitats = {}
for c, repr_seq in cluster2repr.items():
    habitats = [seq2eco[_.split('.')[0]] for _ in cluster2seqs[c]]
    seq_id2habitats[repr_seq] = habitats

seq_id2habitats.update({k: [v] 
                        for k, v in genome2habitat.items() 
                        if not pd.isna(v)})
colors = px.colors.qualitative.Dark24
e2c = {'sediment': {'color': '#E15F99'},
       #  'marine sediment': {'color': ''},
       #  'others': {'color': ''},
       'unknown': {'color': '#222A2A'},
       'terrestrial': {'color': '#1CA71C'},
       'marine': {'color': '#511CFB'},
       'human-made': {'color': '#B68100'},
       'groundwater': {'color': '#750D86'},
       'freshwater': {'color': '#2E91E5'},
       #  'extreme': {'color': ''}
       }
text = to_binary_shape({k: [_.capitalize() for _ in v] for k, v in seq_id2habitats.items()},
                       {k.capitalize(): v for k, v in e2c.items()},
                       info_name="Ecosystem type", 
                       unfilled_other=True,
                       manual_v=['Marine','Groundwater',
                                 "Human-made",'Sediment','Unknown',
                                 'Terrestrial','Freshwater']
                       )
with open('./Ecosystem_type.txt', 'w') as f1:
    f1.write(text)


l2c = {s.split('.')[0]: c for c, seqs in cluster2seqs.items() for s in seqs}
tab.loc[:, 'repr of the corresponding cluster'] = [
    cluster2repr[l2c[s]].split('.')[0] for s in tab.index]
repr_seqs = [_.split('.')[0] for _ in Tree(
    './prepared_sorted.newick', format=3).get_leaf_names()]
repr_seqs = [_ for _ in repr_seqs if _ in tab.index]

tab.loc[:, 'Size of clusters'] = [len(cluster2seqs[l2c[s]]) for s in tab.index]
remained_columns = ['acc', 'date', 'Size of clusters',
                    'Ecosystem type', 'description', 'isolation_source',
                    ]
tab = tab.reindex(index=repr_seqs,
                  columns=remained_columns)
tab.to_excel('./DataSet S2.xlsx')


#
tree = Tree('./prepared_sorted.newick', format=3)
n1 = "I425_S100"
n2 = "I229_S97"
n2node = {n.name: n for n in tree.traverse()}
left_nodes = set(n2node[n1].get_leaf_names()).difference(
    n2node[n2].get_leaf_names())
text = to_color_range({n: 'more basal' for n in left_nodes},
                      {"more basal": "#e53238"})
with open('./more_basal.txt', 'w') as f1:
    f1.write(text)

#
tab = pd.read_excel('./Brocadiae_20142.xlsx', index_col=0)
tab.loc[list(l2group), 'classified lineages'] = [v for k, v in l2group.items()]
tree = Tree(expanduser(
    '~/data/plancto/trees/iqtree/over20p_bac120.formatted.newick'), 3)

group_info = {  # 'I106_S100': 'Anammox',
    'I196_S100': 'Ca.Scalindua',
    "I248_S100": 'Ca.Kuenenia',
    "I304_S100": 'Ca.Jettenia',
    "I368_S100": "Ca.Brocadia",
    "I145_S100": "Basal lineage",
    "I369_S91": "hzsCBA-loss lineage"}
group_info = {[_
               for _ in tree.traverse()
               if _.name == I_id][0]: group
              for I_id, group in group_info.items()}

l2group = {}
for node, group in group_info.items():
    l2group.update({n: group for n in node.get_leaf_names()})

tree = Tree('./prepared_sorted.newick', format=3)
l2group = {k: v for k, v in l2group.items() if k in tree.get_leaf_names()}

info2style = {'Anammox': dict(color="#ff8000"),
              'Ca.Scalindua': dict(color="#7373ff"),
              'Ca.Kuenenia': dict(color="#008000"),
              'Ca.Jettenia': dict(color="#ffff00"),
              "Ca.Brocadia": dict(color="#ff0000"),
              "Basal lineage": dict(color="#B01455"),
              "hzsCBA-loss lineage": dict(color="#88b719"),
              }
text = to_color_range(l2group,
                      info2color={k: list(v.values())[0]
                                  for k, v in info2style.items()},
                      dataset_name='lineage annotations',
                      no_legend=True,)
with open(f"./lineage annotation.txt", 'w') as f1:
    f1.write(text)
