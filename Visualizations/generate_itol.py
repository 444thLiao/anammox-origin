from ete3 import Tree
import os
from tqdm import tqdm
from api_tools import to_label
import pandas as pd

from api_tools.itol_func import to_binary_shape, to_color_range


def read_summary(metadata):
    metadata_df = pd.read_csv(
        metadata, sep='\t', low_memory=False, comment='#', header=None)
    _ = open(metadata)
    header = _.readline()
    header = _.readline().strip('# \n').split('\t')
    metadata_df.columns = header
    metadata_df = metadata_df.set_index("assembly_accession")
    return metadata_df


genome_ids = open(
    '/mnt/home-backup/thliao/plancto/update_Apr2021/over20p_genomes_add_cyano.list').read().strip('\n').split('\n')
# metadata
home = os.getenv("HOME")
metadata = f"{home}/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"
gids = set(genome_ids)
gid2name = {}
for row in tqdm(open(metadata)):
    if not row.startswith("assembly_accession"):
        row = row.split('\t')
        if row[0] in gids:
            gid2name[row[0]] = row[7] + ' ' + row[9]
gid2name = {k: v.replace('.fa.gz', '') for k, v in gid2name.items()}
text = to_label(gid2name)
with open('./update_Apr2021/itol_txt/Genome2names.txt', 'w') as f1:
    f1.write(text)


# for completness
from api_tools.itol_func import color_gradient,get_text_anno
c_df = pd.read_csv('./update_Apr2021/checkM_result_phylum/merged_checkM.tab',sep='\t',index_col=0)
v_df = pd.read_csv("/mnt/home-backup/thliao/verruco/analysis/checkM_result_phylum/merged_checkM.tab",sep='\t',index_col=0)
cy_df = pd.read_csv("/home-user/thliao/data/cyano_basal/analysis/checkM_result_phylum/merged_checkM.tab",sep='\t',index_col=0)

remove_v = set(c_df.index).intersection(v_df.index)
remove_cy = set(c_df.index).intersection(cy_df.index)
c_df = c_df.drop(remove_v); c_df = c_df.drop(remove_cy)

c_df = pd.concat([c_df,v_df,cy_df],axis=0)
l2v = dict(zip(c_df.index,c_df['Completeness']))
t = Tree('./update_Apr2021/trees/iqtree/over20p_bac120_formatted.newick',format=3)
ids = t.get_leaf_names()
text = color_gradient({k:v for k,v in l2v.items() if k in ids})
with open('./update_Apr2021/itol_txt/completeness_phylum.txt','w') as f1:
    f1.write(text)

text = get_text_anno(l2v,extra_replace={"#HEIGHT_FACTOR,1":"HEIGHT_FACTOR\t1.5",
                                           "#HORIZONTAL_GRID,1":"HORIZONTAL_GRID\t0",
                                           "#VERTICAL_GRID,1":"VERTICAL_GRID\t0",
                                           })
with open('./update_Apr2021/itol_txt/completeness_phylum_text.txt','w') as f1:
    f1.write(text)

###
t = Tree('./update_Apr2021/trees/iqtree/over20p_bac120_formatted.newick', format=3)
name2node = {_.name: _ for _ in t.traverse()}
info2style = {  # 'Anammox': dict(color="#ff8000"),
    'Ca.Scalindua': dict(color="#7373ff"),
    'Ca.Kuenenia': dict(color="#008000"),
    'Ca.Jettenia': dict(color="#ffff00"),
    "Ca.Brocadia": dict(color="#ff0000"),
    "Basal lineage": dict(color="#B01455"),
    "hzsCBA-loss lineage": dict(color="#88b719"),
    'Unclassified Anammox bacteria': dict(color="#ff8000"),
}
lineage2n = {'Basal lineage': "I350_S100",
             'Ca.Brocadia': "I1266_S100",
             'Ca.Jettenia': "I1267_S100",
             'Ca.Kuenenia': "I779_S100",
             'Ca.Scalindua': "I565_S100",
             'hzsCBA-loss lineage': "I1135_S100",
             }
S1 = pd.read_excel('./update_Apr2021/DataSet S1.xlsx', index_col=0)
S1.loc[name2node['I263_S100'].get_leaf_names(
), 'classified lineages'] = 'Unclassified Anammox bacteria'
for l, n in lineage2n.items():
    S1.loc[name2node[n].get_leaf_names(), 'classified lineages'] = l
S1.loc[:, 'classified lineages'] = S1.loc[:, 'classified lineages'].fillna('-')
# S1.to_excel('./update_Apr2021/DataSet S1.xlsx',index=1,index_label=S1.index.name)
genome2lineages = dict(zip(S1.index,
                           S1['classified lineages']))
genome2lineages = {k: v for k, v in genome2lineages.items() if str(v) != '-'}

text = to_color_range(genome2lineages,
                      {k: v['color'] for k, v in info2style.items()},
                      dataset_name='Lineages')
with open('./update_Apr2021/itol_txt/Lineage_colorrange.txt', 'w') as f1:
    f1.write(text)


cmd = f"format_newick.py itol-bp -i ./trees/iqtree/over20p_bac120_formatted.newick -o ./itol_txt/over20p_bp"


##
new_df = pd.read_excel('./update_Apr2021/DataSet S1.xlsx',index_col=0)
genome2m = dict(zip(new_df.index,
                    new_df['habitat marine/non']))
genome2m = {k: v for k, v in genome2m.items() if str(v) != '-' and str(v) != 'nan'}

genome2m = {k: ['marine']
            for k,v in genome2m.items() if v in ['marine',
                                                 'groundwater'
                                                 ]}

text = to_binary_shape(genome2m,
                      {'marine':{'color':"#6495ed","shape":"2"}},
                      )
with open('./update_Apr2021/itol_txt/marine_symbols.txt', 'w') as f1:
    f1.write(text)
    
tax_df = pd.read_csv(
    "/home-user/thliao/.cache/ncbi-genome-download/bacteria2taxonomy.tab", sep='\t')
