import plotly.graph_objects as go
from collections import defaultdict
from for_software.for_bayestraits.api.summarize_hmm import filtration_part
import itertools
import os
import re
import string
from os.path import *
from glob import glob
from collections import Counter,defaultdict
import pandas as pd
from api_tools.itol_func import *
from bin.other_convertor.classify_kos import *
from ete3 import Tree
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

def get_ko_info(ko):
    df = pd.read_csv('/mnt/home-backup/thliao/protein_db/kegg/ko_info.tab',sep='\t',index_col=0,header=None)
    df.index = [_.split(':')[-1] for _ in df.index]
    koids = sep_ko(ko)
    
    sub_df = df.reindex(koids)
    text = ""
    for kid,row in sub_df.iterrows():
        des = row[1]
        name = des.split('; ')[0]
        text += '\t'.join([kid,name,'#ff77aa','2',des ]) + '\n'
    print(text)
        
def convert2dict(tab):
    df = pd.read_csv(tab, sep='\t', header=None)
    df = df.sort_values(10)
    df = df.groupby(0).head(1)
    final_dict = {}
    for g in df.loc[:, 1].unique():
        final_dict[g] = dict(zip(['GCA_' + _.split('_')[0].replace('v', '.')
                                  for _ in df.loc[:, 0]],
                                 [1] * df.shape[0]))
    return final_dict
def process_grouping(info2style,sub_df):
    group2ko = defaultdict(list)
    _i2style = {}
    for ko,_d in list(info2style.items()):
        if 'group' in _d and _d['group'].startswith('G:'):
            group_name = _d['group'].split(':')[-1].strip()
            group2ko[group_name].append(ko)
            _b = _d.copy()
            _b['info'] = group_name
            _i2style[group_name] = _b
        else:
            _i2style[ko] = _d.copy()
    for g,ko_list in group2ko.items():
        sub_df.loc[:,g] = 0
        sub_df.loc[sub_df.loc[:,ko_list].sum(1) >= len(ko_list)/2,g] = 1
    sub_df = sub_df.loc[:,_i2style]
    return _i2style,sub_df

tab = "/mnt/home-backup/thliao/plancto/protein_annotations/join_all_6dbs.tab"
all_df = pd.read_csv(tab, sep='\t', index_col=0)
def draw_pattern(intab="./gene_GandL/kegg_db/genera_specific/draw_genes.tab",
                 odir = "./gene_GandL/itol_txt/",
                 anammox_only=True,
                 ID2color=None):
    # intab = "./gene_GandL/kegg_db/genera_specific/draw_genes.tab"
    #tab = "/mnt/home-backup/thliao/plancto/protein_annotations/hmmsearch_merged/merged_hmm_binary.tab"
    # intab =
    if anammox_only:
        draw_gids = open(
            "/home-user/thliao/data/plancto/gene_GandL/anammox_ids.list").read().split('\n')
        draw_gids = [_ for _ in draw_gids if _]
    else:
        draw_gids = []
    f1 = open(intab).read()
    header = ['KO number', 'info', 'color', 'shape','group']
    for section in f1.split('\n#'):
        sections = section.split('\n')
        section_name = sections[0].strip()
        if not section_name:
            continue
        contents = [_ for _ in sections[1:] if _.strip()]
        info2style = {}
        for row in contents:
            rows = row.strip('\n').split('\t')
            _dict = dict(zip(header, rows))
            if 'info' not in _dict:
                _dict['info'] = _dict['KO number']
                _dict['shape'] = '2'
            if ID2color is not None and _dict['KO number'] in ID2color:
                _dict['color'] = ID2color[_dict['KO number']]
            _dict.pop('KO number')
            info2style[rows[0]] = _dict
            
        # if section_name=='oxidoreductase and ATPase':
        #     break
        # preprocess
        ko_list = list(info2style)
        if draw_gids:
            sub_df = all_df.reindex(index=draw_gids, columns=ko_list).fillna(0)
        else:
            sub_df = all_df.reindex(columns=ko_list).fillna(0)
        info2style,sub_df = process_grouping(info2style,sub_df)
        ko_list = list(info2style)
        id2info = {k: [_
                       for _, _v in g.items() if _v != 0]
                   for k, g in sub_df.to_dict(orient='index').items()}
        print(info2style)
        text = to_binary_shape(id2info,
                               info2style=info2style,
                               manual_v=list(ko_list),
                               info_name=section_name,
                               same_color='#0077aa',
                               no_legend=True,
                               extra_replace={'#MARGIN,0': 'MARGIN\t50',
                                              "#SHOW_LABELS,1": "SHOW_LABELS\t1\nSIZE_FACTOR\t2\n",
                                              "#SYMBOL_SPACING,10": "SYMBOL_SPACING\t15"})
        if not exists(odir):
            os.system(f"mkdir -p {odir}")
        with open(f"{odir}/{section_name.replace('/', ' or ')}.txt", 'w') as f1:
            f1.write(text)

def sep_ko(ko_list):
    kos = []
    if isinstance(ko_list,list):
        return ko_list
    for _ in ko_list.split('K'):
        if _.strip(string.printable[10:]).strip(' '):
            kos.append('K' + _.strip(string.printable[10:]).strip(' '))
    return kos


def get_count(df, tree, ko_list, group1, group2):
    LCA = tree.get_common_ancestor(group1.split(','))
    larger_LCA = tree.get_common_ancestor(group2.split(','))

    g = larger_LCA.get_leaf_names()
    g1 = LCA.get_leaf_names()
    # target group
    g2 = set(g).difference(set(g1))

    kos = sep_ko(ko_list)

    remained_kos = []
    for _ in kos:
        if _ in df.columns:
            remained_kos.append(_)
    s1 = df.loc[g1, remained_kos].sum()
    s2 = df.loc[g2, remained_kos].sum()
    print(f"total: {len(g1)}, present in ({s1 / len(g1) * 100}) genomes")
    print(f"total: {len(g2)}, present in ({s2 / len(g2) * 100}) genomes")


def fe_text(g1, g2, ko_df):
    ko2tab = {}
    ko2odd_ratio = {}
    ko2pvalue = {}
    for ko in tqdm(set(ko_df.columns)):
        ko_g1_p = sum(ko_df.loc[g1, ko] > 0)
        ko_g1_a = sum(ko_df.loc[g1, ko] == 0)
        ko_g2_p = sum(ko_df.loc[g2, ko] > 0)
        ko_g2_a = sum(ko_df.loc[g2, ko] == 0)
        tab_num = [[ko_g1_p, ko_g1_a],
                   [ko_g2_p, ko_g2_a]]
        ko2tab[ko] = tab_num
        oddratio, p = fisher_exact(tab_num, alternative="two-sided")
        ko2pvalue[ko] = p
        ko2odd_ratio[ko] = oddratio

    fe_corrected_p = multipletests([_
                                    for k, _ in ko2pvalue.items()],
                                   method='fdr_bh')
    fe_corrected_ko2p = dict(zip([k
                                  for k, _ in ko2pvalue.items()],
                                 fe_corrected_p[1]))
    sig_ko_list = {k: v
                   for k, v in fe_corrected_ko2p.items()
                   if v <= 0.05}
    return ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list

def main(annotation_df, tree, group1, group2, ofile=None,
         stop_at_annotation=False, just_g2=False,
         pass_id=False, return_df=False):
    """

    :param annotation_df:
    :param tree:
    :param group1: some specific group which is smaller
    :param group2: bigger group which normally comprise of group1
    :param ofile:
    :return:
    """
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))

    if pass_id:
        g1 = group1
        g2 = group2
    else:
        LCA = tree.get_common_ancestor(group1.split(','))
        larger_LCA = tree.get_common_ancestor(group2.split(','))

        g = larger_LCA.get_leaf_names()
        g1 = LCA.get_leaf_names()
        # target group
        if not just_g2:
            g2 = set(g).difference(set(g1))
        else:
            g2 = set(g)

    # remainning group
    ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = fe_text(
        g1, g2, annotation_df)

    if stop_at_annotation:
        return ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list

    # summarized df
    br_kos = ko_classified_br(sig_ko_list)
    md_kos = ko_classified_module(sig_ko_list)
    md2info = get_md_infos(md_kos.keys())
    info_kos = get_ko_infos(sig_ko_list)
    infos = get_br_info(br_kos)

    if len(infos) == 0:
        raise Exception
    
    df = pd.concat(infos, axis=0)
    df.loc[:, 'des'] = [info_kos.get(_, '') for _ in df.index]

    new_df = pd.DataFrame(index=sig_ko_list, columns=['num_p in target', 'num_a in target', "ratio in target (%)",
                                                      'num_p in remained', 'num_a in remained', "ratio in remained (%)",
                                                      ])
    for ko in new_df.index:
        row = ko2tab[ko][0][::]
        row.append(round(row[0] / (row[0] + row[1]) * 100, 4))
        row += ko2tab[ko][1]
        row.append(round(row[3] / (row[4] + row[3]) * 100, 4))
        new_df.loc[ko, :] = row
    new_df = new_df.reindex(index=df.index)

    final_df = pd.concat([new_df, df], axis=1)
    for md, kos in tqdm(md_kos.items()):
        kos = [_ for _ in kos if _ in final_df.index]
        try:
            final_df.loc[kos, "module ID"] = md
            final_df.loc[kos, "module name"] = md2info[md]['name']
            final_df.loc[kos, "module all ko"] = md2info[md]['all ko']
        except:
            pass
            print(kos)
    final_df = final_df.sort_values(
        ['num_p in target', 'top', 'A', 'B', 'C', 'des'])

    final_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in final_df.index]
    final_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_]
                                           for _ in final_df.index]

    if return_df:
        return final_df
    else:
        final_df.to_excel(ofile, index_label='K number')


# all gene table
tab = "/mnt/home-backup/thliao/plancto/protein_annotations/hmmsearch_merged/merged_hmm_binary.tab"
all_df = pd.read_csv(tab, sep='\t', index_col=0)


# for dating tree
intree = './trees/iqtree/over20p_bac120.formatted.newick'
t = Tree(intree, 3)
## anammox group
LCA = t.get_common_ancestor(['GCA_001824525.1','GCA_008933285.1'])
larger_LCA = t.get_common_ancestor(['GCA_004376375.1','GCA_008933285.1'])
g = larger_LCA.get_leaf_names()
g1 = LCA.get_leaf_names()
# target group
g2 = set(g).difference(set(g1))


all2p = {}
for row in tqdm(f1.split('\n')):
    if row and not row.startswith('#'):
        ID = row.split('\t')[0]
        if ID not in all_df.columns:
            all2p.update({ID:1})
            continue
        ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = fe_text(
        g1, g2, all_df.loc[:,[ID]])
        
        all2p.update(ko2pvalue)
ID2color = {k: '#1976D2' if k in all2cp else '#558B2F' for k,v in all2p.items()}
draw_pattern(intab='./gene_GandL/draw_genes.tab',
             odir='./gene_GandL/itol_txt/othercolors',ID2color=ID2color,anammox_only=False)

main(all_df,
     Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3),
     group1='GCA_001824525.1,GCA_008933285.1',
     group2='GCA_004376375.1,GCA_008933285.1',
     ofile='./gene_GandL/kegg_db/over20p_gain_loss.xlsx')

main(annotation_df=all_df,
     tree=Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3),
     group1='GCA_001824525.1,GCA_008933285.1',
     group2='GCA_005239935.1,GCA_003452155.1',
     ofile='./gene_GandL/kegg_db/over20p_biglineage_gain_loss.xlsx')

main(all_df,
     Tree('/mnt/home-backup/thliao/plancto/trees/final/83g_merged_sorted.newick', 3),
     group1='GCA_001828545.1,GCA_001828295.1',
     group2='GCA_003694635.1,GCA_003551565.1',
     ofile='./gene_GandL/dating_83g.xlsx')

## compared group

group_info_dict = {'anammox': 'GCA_001824525.1,GCA_008933285.1',
                   'Scalindua': 'GCA_004282745.1,GCA_007859995.1',
                   'Kuenenia': "GCA_002418245.1,GCA_900696675.1",
                   'Jettenia': "GCA_900696475.1,GCA_900696655.1",
                   "Brocadia": "GCA_001828415.1,GCA_001828605.1",
                   "Basal": "GCA_001828605.1,GCA_001828415.1"}
group_info = ['Scalindua',
              'Kuenenia',
              'Jettenia',
              "Brocadia",
              "Basal"]


# tree = Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3)
# for g1, g2 in itertools.combinations(group_info, 2):
#     main(all_df.reindex(columns=[_ for _ in all_df.columns if _.startswith('K')]),
#          Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3),
#          group1=group_info_dict[g1],
#          group2=group_info_dict[g2],
#          just_g2=True,
#          ofile=f'./gene_GandL/intergeneric/{g1}_{g2}.xlsx')


def add_count(ori_df, tree, group1, group2, id_list, extra_df: pd.DataFrame = None):
    if len(id_list) == 0:
        return pd.DataFrame()
    LCA = tree.get_common_ancestor(group1.split(','))
    larger_LCA = tree.get_common_ancestor(group2.split(','))

    g = larger_LCA.get_leaf_names()
    g1 = LCA.get_leaf_names()
    # target group
    g2 = set(g).difference(set(g1))
    others = set(tree.get_leaf_names()).difference(set(g1).union(set(g2)))
    # remainning group
    if extra_df is not None:
        count_df = extra_df.reindex(id_list)
    else:
        count_df = pd.DataFrame(index=id_list)
    count_df.loc[:, 'ratio in target (%)'] = 0
    count_df.loc[:, 'ratio in remained (%)'] = 0
    # count_df.loc[:, 'present in other'] = 0

    g1_df = ori_df.reindex(columns=id_list, index=g1)
    g2_df = ori_df.reindex(columns=id_list, index=g2)
    other_df = ori_df.reindex(columns=id_list, index=others)
    count_df.loc[:, 'num_p in target'] = g1_df.sum(0)
    count_df.loc[:, 'num_a in target'] = len(g1) - g1_df.sum(0)
    count_df.loc[:, 'ratio in target (%)'] = g1_df.sum(0) / len(g1)

    count_df.loc[:, 'num_p in remained'] = g2_df.sum(0)
    count_df.loc[:, 'num_a in remained'] = len(g2) - g1_df.sum(0)
    count_df.loc[:, 'ratio in remained (%)'] = g2_df.sum(0) / len(g2)
    # count_df.loc[:, 'present in other'] = other_df.sum(0) / len(others)
    # count_df.loc[:, 'sub1'] = count_df.loc[:, 'ratio in target (%)'] - count_df.loc[:, 'present in neighbour']
    # count_df.loc[:, 'sub2'] = count_df.loc[:, 'present in target'] - count_df.loc[:, 'present in other']
    # count_df = count_df.sort_values(['sub1', 'sub2'], ascending=False)
    return count_df


## go enrichment for kegg
def get_GO_enrichment(sig_kegg, pop_kegg, ofile, dtype="kegg", extra_parameters=''):
    tmp_dir = "./tmp/GO_enrich"
    if not exists(tmp_dir):
        os.system(f"mkdir -p {tmp_dir}")
    study_file = join(tmp_dir, "study")
    pop_file = join(tmp_dir, "populations")

    with open(study_file, "w") as f1:
        f1.write("\n".join(list(set(sig_kegg))))
    with open(pop_file, "w") as f1:
        f1.write("\n".join(list(set(pop_kegg))))
    if dtype == "kegg":
        asso_file = "/home-user/thliao/data/plancto/gene_GandL/asso_files/kegg_associations.txt"
    elif dtype == "interpro":
        asso_file = "/home-user/thliao/data/plancto/gene_GandL/asso_files/ipr2go_ass.txt"
    if not exists(dirname(ofile)):
        os.system(f"mkdir -p {dirname(ofile)}")
    cmd = f"python3 ~/software/goatools/scripts/find_enrichment.py {study_file} {pop_file} {asso_file} --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile={ofile} --obo /home-user/thliao/software/goatools/go-basic.obo " + extra_parameters
    os.system(cmd)


## GO enrichment between anammox group and non-anammox group
kegg_df = add_count(all_df,
                    tree=Tree(
                        './trees/iqtree/over20p_bac120.formatted.newick', 3),
                    group1='GCA_001824525.1,GCA_008933285.1',
                    group2='GCA_007131925.1,GCA_008933285.1',
                    id_list=all_df.columns, )
kegg_df_sub = pd.read_excel("./gene_GandL/kegg_db/over20p_gain_loss.xlsx")
kegg_df = kegg_df.reindex(kegg_df_sub.iloc[:, 0])
for cal_type in ["Gain", "Loss"]:
    if cal_type == "Gain":
        study_kegg = list(
            kegg_df.index[(kegg_df["sub1"] > 0) & (kegg_df["sub2"] > 0)])
    else:
        study_kegg = list(
            kegg_df.index[(kegg_df["sub1"] < 0) & (kegg_df["sub2"] < 0)])
    # significant Loss?
    draw_gids = open(
        "/home-user/thliao/data/plancto/gene_GandL/near_anammox_ids.list").read().split('\n')
    draw_gids = [_ for _ in draw_gids if _]
    sub_df = all_df.reindex(draw_gids)
    sub_df = sub_df.loc[:, sub_df.sum() != 0]
    pop_kegg = list(sub_df.columns)
    # get_GO_enrichment(study_kegg, pop_kegg, "./gene_GandL/kegg_db/Gain_anammox.xlsx")
    get_GO_enrichment(study_kegg, pop_kegg, f"./gene_GandL/kegg_db/{cal_type}_anammox.xlsx",
                      extra_parameters=" --sections=goatools.test_data.sections.dvk2018_genefamilyclusters --method=fdr_bh")
## compare group (design)
# Jettenia vs other anammox (example)
collect_dfs = {}
for genus in group_info:
    ofile = f'./gene_GandL/kegg_db/genera_specific/{genus}.xlsx'
    df = main(all_df.reindex(columns=[_ 
                                      for _ in all_df.columns 
                                      if _.startswith('K')]),
              Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3),
              group1=group_info_dict[genus],
              group2=group_info_dict['anammox'],
              ofile=ofile,
              return_df=True)
    collect_dfs[genus] = df
with pd.ExcelWriter(f"./gene_GandL/kegg_db/genera_specific.xlsx") as writer:
    for genus, df in collect_dfs.items():
        df.to_excel(writer, sheet_name=genus, index_label="K number", index=1)


# enrichment analysis
draw_gids = open(
    "/home-user/thliao/data/plancto/gene_GandL/anammox_ids.list").read().split('\n')
draw_gids = [_ for _ in draw_gids if _]
sub_df = all_df.reindex(draw_gids)
sub_df = sub_df.loc[:, sub_df.sum() != 0]
pop_kegg = list(sub_df.columns)

dfs = pd.read_excel(
    f'./gene_GandL/kegg_db/genera_specific.xlsx', sheet_name=None, index_col=0)
for genus in group_info:
    # tree = Tree('./trees/iqtree/over20p_bac120.formatted.newick', 3)
    # LCA = tree.get_common_ancestor(group_info_dict[genus].split(','))
    # gids = LCA.get_leaf_names()
    # sub_df = all_df.reindex(draw_gids)
    # sub_df = sub_df.loc[:, sub_df.sum() != 0]
    # pop_kegg = list(sub_df.columns)
    for cal_type in ["Gain", "Loss"]:
        kegg_df = dfs[genus].drop_duplicates('des')
        _t = kegg_df["ratio in target (%)"] 
        _t2 = kegg_df["ratio in remained (%)"]
        study_kegg = kegg_df.index[_t > _t2] if cal_type == "Gain" else kegg_df.index[_t <= _t2]
        print(genus,kegg_df.shape[0],cal_type,len(study_kegg))
        get_GO_enrichment(study_kegg, pop_kegg, f"./gene_GandL/kegg_db/genera_specific/GO_enrichment/{cal_type}_{genus}.xlsx",
                          extra_parameters=" --method=fdr_bh")

##### for others database
mixed_df = pd.read_csv(
    './protein_annotations/merged_interpro/mixed_binary.tab', sep='\t', index_col=0)
mixed_df2 = pd.read_csv(
    './protein_annotations/merged_interpro/mixed_interpro_binary.tab', sep='\t', index_col=0)
mixed_df = pd.concat([mixed_df, mixed_df2], axis=0)
mixed_df = mixed_df.T

# actually not kegg... it should be id from other database
ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = main(annotation_df=mixed_df,
                                                    tree=Tree(
                                                        './trees/iqtree/over20p_bac120.formatted.newick', 3),
                                                    group1='GCA_001824525.1,GCA_008933285.1',
                                                    group2='GCA_004376375.1,GCA_008933285.1',
                                                    ofile='./gene_GandL/over20p_gain_loss.xlsx',
                                                    stop_at_annotation=True)


def pack_up(sig_ko_list, ko2pvalue, ofile, ori_df=None, all_one=False,
            g1='GCA_001824525.1,GCA_008933285.1',  # anammox
            g2='GCA_004376375.1,GCA_008933285.1'):  # near anammox

    if all_one:
        db_df = add_count(ori_df,
                          tree=Tree(
                              './trees/iqtree/over20p_bac120.formatted.newick', 3),
                          group1=g1,
                          group2=g2, )
        return db_df

    pfam_list = [_ for _ in sig_ko_list if _.startswith('PF')]
    cdd_list = [_ for _ in sig_ko_list if _.startswith('cd')]
    interpro_list = [_ for _ in sig_ko_list if _.startswith('IPR')]
    TIGRFAM_list = [_ for _ in sig_ko_list if _.startswith('TIGR')]
    print(len(cdd_list), len(interpro_list), len(pfam_list), len(TIGRFAM_list))

    # pfam_df
    pfam_df = pd.read_csv(
        '/home-user/thliao/data/protein_db/mapping_files/PFAM_v32.tab', sep='\t', index_col=0)
    pfam_df.index = [_.split('.')[0]
                     for _ in pfam_df.index]  # remove version number
    pfam_df = pfam_df.reindex(pfam_list)
    pfam_df = add_count(mixed_df,
                        tree=Tree(
                            './trees/iqtree/over20p_bac120.formatted.newick', 3),
                        group1=g1,
                        group2=g2,
                        id_list=pfam_df.index,
                        extra_df=pfam_df)
    pfam_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in pfam_df.index]
    pfam_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_]
                                          for _ in pfam_df.index]

    cdd_df = pd.read_csv(
        '/home-user/thliao/software/interproscan-5.38-76.0/data/cdd/3.17/data/cddid.tbl', sep='\t', header=None)
    cdd_df = cdd_df.set_index(1)
    cdd_df = cdd_df.reindex(cdd_list)
    cdd_df = cdd_df.reindex(columns=[2, 3])
    cdd_df.columns = ["Name", "des"]
    cdd_df = add_count(mixed_df,
                       tree=Tree(
                           './trees/iqtree/over20p_bac120.formatted.newick', 3),
                       group1=g1,
                       group2=g2,
                       id_list=cdd_df.index,
                       extra_df=cdd_df)
    cdd_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in cdd_df.index]
    cdd_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_] for _ in cdd_df.index]

    tigrfam_df = pd.read_csv(
        '/home-user/thliao/data/protein_db/mapping_files/TIGRFAM_v15.tab', sep='\t', index_col=0)
    tigrfam_df = tigrfam_df.reindex(TIGRFAM_list)
    tigrfam_df = add_count(mixed_df,
                           tree=Tree(
                               './trees/iqtree/over20p_bac120.formatted.newick', 3),
                           group1=g1,
                           group2=g2,
                           id_list=tigrfam_df.index,
                           extra_df=tigrfam_df)
    tigrfam_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in tigrfam_df.index]
    tigrfam_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_]
                                             for _ in tigrfam_df.index]

    interpro_df = pd.read_csv(
        '/mnt/home-backup/thliao/protein_db/interpro/entry.list', sep='\t', index_col=0)
    interpro_df = interpro_df.reindex(interpro_list)
    interpro_df = add_count(mixed_df,
                            tree=Tree(
                                './trees/iqtree/over20p_bac120.formatted.newick', 3),
                            group1=g1,
                            group2=g2,
                            id_list=interpro_df.index,
                            extra_df=interpro_df)
    interpro_df.loc[:, 'corrected pvalue'] = [ko2pvalue[_]
                                              for _ in interpro_df.index]
    interpro_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_]
                                              for _ in interpro_df.index]

    # enrichment analysis
    # python3 ~/software/goatools/scripts/find_enrichment.py ~/data/plancto/gene_GandL/GO_enrichment/distinct_interpro_study.go.list ~/data/plancto/gene_GandL/GO_enrichment/anammox_populations.go.list ~/data/plancto/gene_GandL/GO_enrichment/ipr2go_ass.txt   --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_id2gos.xlsx --obo ~/software/goatools/scripts/go-basic.obo

    # odir = './gene_GandL/other_db'
    odir = dirname(ofile)
    os.makedirs(f'{odir}', exist_ok=True)

    with pd.ExcelWriter(ofile) as writer:
        pfam_df.to_excel(writer, sheet_name="Pfam", index_label='Pfam ID')
        cdd_df.to_excel(writer, sheet_name="CDD", index_label='CDD ID')
        tigrfam_df.to_excel(writer, sheet_name="TIGRFAM",
                            index_label='TIGRFAM ID')
        interpro_df.to_excel(writer, sheet_name="interpro",
                             index_label='Interpro ID')


# anammox itself
pfam_list = [_ for _ in mixed_df.columns if str(_).startswith('PF')]
cdd_list = [_ for _ in mixed_df.columns if str(_).startswith('cd')]
interpro_list = [_ for _ in mixed_df.columns if str(_).startswith('IPR')]
TIGRFAM_list = [_ for _ in mixed_df.columns if str(_).startswith('TIGR')]

all_sig_ko_list = {}
all_ko2pvalue = {}
for db in [pfam_list, cdd_list, interpro_list, TIGRFAM_list]:
    _1, _2, ko2pvalue, sig_ko_list = main(mixed_df.reindex(columns=db),
                                          Tree(
                                              './trees/iqtree/over20p_bac120.formatted.newick', 3),
                                          group1='GCA_001824525.1,GCA_008933285.1',
                                          group2='GCA_004376375.1,GCA_008933285.1',
                                          ofile='./',  # actually useless here
                                          stop_at_annotation=True)
    all_sig_ko_list.update(sig_ko_list)
    all_ko2pvalue.update(ko2pvalue)
pack_up(all_sig_ko_list, all_ko2pvalue,
        ofile=f'./gene_GandL/other_db/anammox.xlsx', ori_df=mixed_df)

## other db (genus specific)
for genus in group_info:
    ofile = f'./gene_GandL/other_db/genus_specific/{genus}'
    _1, _2, ko2pvalue, sig_ko_list = main(mixed_df,
                                          Tree(
                                              './trees/iqtree/over20p_bac120.formatted.newick', 3),
                                          group1=group_info_dict[genus],
                                          group2=group_info_dict['anammox'],
                                          ofile='./',  # actually useless here
                                          stop_at_annotation=True)
    pfam_list = [_ for _ in sig_ko_list if str(_).startswith('PF')]
    cdd_list = [_ for _ in sig_ko_list if str(_).startswith('cd')]
    interpro_list = [_ for _ in sig_ko_list if str(_).startswith('IPR')]
    TIGRFAM_list = [_ for _ in sig_ko_list if str(_).startswith('TIGR')]

    # print(genus,len(cdd_list),len(interpro_list),len(pfam_list),len(TIGRFAM_list))
    pack_up(sig_ko_list, ko2pvalue,
            ofile=ofile + '.xlsx',
            g1=group_info_dict[genus],
            g2=group_info_dict['anammox'],
            )

# construct relationships among dbs

dbs = ['Pfam', 'interpro', 'CDD', 'TIGRFAM', 'kegg', "TCDB"]
cross_db_df = pd.read_csv(
    './protein_annotations/all_merged_20200328.tab', sep='\t', index_col=0)
# cross_db_df = cross_db_df.loc[cross_db_df.Pfam.isin(pfam_list) | cross_db_df.CDD.isin(cdd_list) | cross_db_df.TIGRFAM.isin(tigrfam_df) | cross_db_df.interpro.isin(interpro_df)]

nodes = [_ for _db in dbs for _ in cross_db_df.loc[:, _db].unique()]
edges = []
_df = cross_db_df.loc[cross_db_df.loc[:, dbs].count(1) > 1, dbs]
_df = _df.drop_duplicates()
for _, row in tqdm(_df.iterrows(), total=_df.shape[0]):
    v = [_ for _ in row if not pd.isna(_)]
    edges += ['\t'.join(sorted(_)) for _ in list(itertools.combinations(v, 2))]

edges = set(edges)
with open('./protein_annotations/all_edges.tab', 'w') as f1:
    f1.write('\n'.join(edges))

for genus in tqdm(group_info):
    ofile = f'./gene_GandL/other_db/genus_specific/{genus}.xlsx'
    collect_dfs = pd.read_excel(ofile, sheet_name=None, index_col=0)

    for db, df in list(collect_dfs.items()):
        if df.shape[0] == 0:
            continue
        if all([type(_) == int for _ in df.index]):
            df = df.set_index(df.columns[0])
        df.loc[:, 'dup_kegg'] = ''
        for idx, row in df.iterrows():
            s = [_ for _ in edges if idx in _.split('\t')]
            s = [_v for _ in s for _v in _.split('\t') if _v.startswith('K')]
            if s:
                df.loc[idx, 'dup_kegg'] = ';'.join(sorted(list(set(s))))
        collect_dfs[db] = df
        
    with pd.ExcelWriter(ofile) as writer:
        for db, df in collect_dfs.items():
            df.to_excel(writer, sheet_name=db, index=True)


# fix the add count problem
ofile = f'./gene_GandL/kegg_db/genera_specific.xlsx'
collect_dfs = pd.read_excel(ofile, sheet_name=None, index_col=0)
for db, df in list(collect_dfs.items()):
    new_df = add_count(new_mixed_df,
                       tree=Tree(
                           './trees/iqtree/over20p_bac120.formatted.newick', 3),
                       group1=group_info_dict[db],
                       group2=group_info_dict["anammox"],
                       id_list=df.index,
                       extra_df=df)
with pd.ExcelWriter(ofile) as writer:
    for db, df in collect_dfs.items():
        df.to_excel(writer, sheet_name=db, index=True)

new_mixed_df = pd.concat([mixed_df, merged_TCDB_df], axis=1)
for genus in tqdm(group_info):
    ofile = f'./gene_GandL/other_db/genus_specific/{genus}.xlsx'
    collect_dfs = pd.read_excel(ofile, sheet_name=None, index_col=0)
    for db, df in list(collect_dfs.items()):
        new_df = add_count(new_mixed_df,
                           tree=Tree(
                               './trees/iqtree/over20p_bac120.formatted.newick', 3),
                           group1=group_info_dict[genus],
                           group2=group_info_dict["anammox"],
                           id_list=df.index,
                           extra_df=df)
        collect_dfs[db] = new_df
    with pd.ExcelWriter(ofile) as writer:
        for db, df in collect_dfs.items():
            df.to_excel(writer, sheet_name=db, index=True)

# enrichment analysis
draw_gids = open(
    "/home-user/thliao/data/plancto/gene_GandL/anammox_ids.list").read().split('\n')
draw_gids = [_ for _ in draw_gids if _]
sub_df = mixed_df.reindex(draw_gids)
sub_df = sub_df.loc[:, sub_df.sum() != 0]
pop_IPR = list([_ for _ in sub_df.columns if str(_).startswith("IPR")])
for genus in group_info:
    odir = f'./gene_GandL/other_db/genus_specific/{genus}'
    df_file = f'{odir}/interpro_sig.xlsx'
    df = pd.read_excel(df_file, index_col=0)
    for cal_type in ["Gain", "Loss"]:
        IPR_df = pd.read_excel(df_file, index_col=0)
        study_IPR = IPR_df.index[IPR_df["sub1"] >
                                 0] if cal_type == "Gain" else IPR_df.index[IPR_df["sub1"] < 0]
        get_GO_enrichment(study_IPR, pop_IPR, f"{odir}/{cal_type}_Interpro.xlsx",
                          dtype="interpro",
                          extra_parameters=" --sections=goatools.test_data.sections.dvk2018_genefamilyclusters --method=fdr_bh")

## use TCDB to annotated (add new database)
# cmd:
# python3 ~/bin/batch_run/batch_any.py -i ./tmp_faa -o ./TCDB_out_20200327/ -s faa -ns tab -np 10 -cmd 'blastp -query {infile} -db /mnt/home-backup/thliao/protein_db/TCDB/tcdb -outfmt 6 -evalue 1e-10 -num_threads 3 > {ofile} '


# cmd (diamond)  abandon
# python3 ~/bin/batch_run/batch_any.py -i ./tmp_faa -o ./TCDB_out_20200327_diamond/ -s faa -ns tab -np 10 -cmd 'diamond blastp -q {infile} --db /mnt/home-backup/thliao/protein_db/TCDB/tcdb -f 6 --evalue 1e-10 --threads 3 -o {ofile} '

all_tabs = [_ for _ in glob("./*.tab") if _.startswith('./GCA')]
dfs = [pd.read_csv(_, sep='\t', header=None)
       for _ in all_tabs if os.path.getsize(_) != 0]
dfs = [df.loc[df[10] <= 1e-20, :] for df in dfs]
dfs = [df.sort_values(10, ascending=True).groupby(0).head(1)
       for df in dfs]
genome2info = defaultdict(list)
for df in tqdm(dfs):
    genome = "GCA_" + df.iloc[0, 0].split('_')[0].replace('v', '.')
    for _, row in df.iterrows():
        genome2info[genome].append(
            (row[0], row[1].split('|')[-1], float(row[10])))
post_filtered = filtration_part(genome2info, evalue=1e-20)
odir = './merged_result/'
if not exists(odir):
    os.makedirs(odir)
ofile_info = join(odir, "merged_TCDB_info.tab")
ofile_binary = join(odir, "merged_TCDB_binary.tab")
ofile_num = join(odir, "merged_TCDB_num.tab")
final_df = pd.DataFrame.from_dict(post_filtered, orient='index')
bin_df = final_df.applymap(lambda x: 0 if pd.isna(x) else 1)
num_df = final_df.applymap(lambda x: 0 if pd.isna(x)
                           else len(str(x).split(',')))
final_df.to_csv(ofile_info, sep='\t', index=1)
bin_df.to_csv(ofile_binary, sep='\t', index=1)
num_df.to_csv(ofile_num, sep='\t', index=1)
TCDB_df = pd.read_csv(
    '/home-user/thliao/data/plancto/protein_annotations/TCDB_out_20200327/merged_result/merged_TCDB_binary.tab', sep='\t', index_col=0)

# collapse df
mapping = {c: c.rpartition('.')[0].rpartition('.')[0]
           for c in TCDB_df.columns}
merged_df = TCDB_df.groupby(mapping, axis=1).sum()
merged_TCDB_df = merged_df.where(merged_df == 0, 1)

extra_info = open(
    '/home-user/thliao/data/protein_db/TCDB/tcdb.dr').read().split('\n')
id2info = {k.split(';')[1].strip(): k.split(';')[-1].strip()
           for k in extra_info if k}
_df = pd.DataFrame.from_dict(id2info, orient='index')
_df.columns = ['des']

id2info = {k.split(';')[1].strip().rpartition('.')[0].rpartition('.')[0]: k.split(';')[-1].strip()
           for k in extra_info if k}
_df = pd.DataFrame.from_dict(id2info, orient='index')
_df.columns = ['des']

# genus
odir = "/home-user/thliao/data/plancto/gene_GandL/other_db/genus_specific"
for genus in group_info:
    collect_dfs = {}
    ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = main(merged_TCDB_df,
                                                        tree=Tree(
                                                            './trees/iqtree/over20p_bac120.formatted.newick', 3),
                                                        group1=group_info_dict[genus],
                                                        group2=group_info_dict['anammox'],
                                                        ofile='./',
                                                        stop_at_annotation=True)
    print(genus, len(sig_ko_list))
    if len(sig_ko_list) == 0:
        print(genus)
        continue
    db_df = add_count(merged_TCDB_df,
                      tree=Tree(
                          './trees/iqtree/over20p_bac120.formatted.newick', 3),
                      group1=group_info_dict[genus],
                      group2=group_info_dict['anammox'],
                      id_list=list(sig_ko_list),
                      extra_df=_df)
    db_df.index.name = "TC number"

    db_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in db_df.index]
    db_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_] for _ in db_df.index]
    ofile = f'./gene_GandL/other_db/genus_specific/{genus}.xlsx'
    if exists(ofile):
        collect_dfs.update(pd.read_excel(ofile, sheet_name=None, index_col=0))
    collect_dfs["TCDB"] = db_df
    with pd.ExcelWriter(ofile) as writer:
        for db, df in collect_dfs.items():
            df.to_excel(writer, sheet_name=db, index=True)
# anammox itself

ko2tab, ko2odd_ratio, ko2pvalue, sig_ko_list = main(merged_TCDB_df,
                                                    tree=Tree(
                                                        './trees/iqtree/over20p_bac120.formatted.newick', 3),
                                                    group1='GCA_001824525.1,GCA_008933285.1',
                                                    group2='GCA_004376375.1,GCA_008933285.1',
                                                    ofile='./',
                                                    stop_at_annotation=True)
db_df = add_count(merged_TCDB_df,
                  tree=Tree(
                      './trees/iqtree/over20p_bac120.formatted.newick', 3),
                  group1='GCA_001824525.1,GCA_008933285.1',
                  group2='GCA_004376375.1,GCA_008933285.1',
                  id_list=list(sig_ko_list),
                  extra_df=_df)
db_df.index.name = "TC number"
db_df.loc[:, 'raw pvalue'] = [ko2pvalue[_] for _ in db_df.index]
db_df.loc[:, 'corrected pvalue'] = [sig_ko_list[_] for _ in db_df.index]
ofile = f'./gene_GandL/other_db/anammox.xlsx'
if exists(ofile):
    collect_dfs.update(pd.read_excel(ofile, sheet_name=None, index_col=0))
collect_dfs["TCDB"] = db_df
with pd.ExcelWriter(ofile) as writer:
    for db, df in collect_dfs.items():
        df.to_excel(writer, sheet_name=db, index=True)
## draw kegg pathway

df = pd.read_excel("./gene_GandL/kegg_db/over20p_gain_loss.xlsx")
df = df.set_index("KO number")
t1 = df["ratio in target (%)"]
t2 = df["ratio in remained (%)"]
kos = set(df.index[(t1 >= 0.7) & (t2 <= 0.3)])

# gained genes
ko2p = ko2pathway(kos)
tmp = {k: v for _, _k in ko2p.items()
       for k, v in _k.items()}
counter_pathway = [k for _1, _ in ko2p.items() for k in list(_)]

used_list = sorted(list(set(tmp)), key=lambda x: Counter(counter_pathway)[x])
for _ in used_list:
    print(_, Counter(counter_pathway)[
          _], tmp[_], '+'.join([_] + [k for k, ko in ko2p.items() if _ in list(ko)]))

fig = go.Figure()
fig.add_bar(x=[Counter(counter_pathway)[_] for _ in used_list],
            y=[tmp[_] for _ in used_list],
            orientation='h')
fig.write_html('./test.html')

kos = set(df.index[(t1 <= 0.3)])

# lost genes
ko2p = ko2pathway(kos)
tmp = {k: v for _, _k in ko2p.items()
       for k, v in _k.items()}
counter_pathway = [k for _1, _ in ko2p.items() for k in list(_)]

used_list = sorted(list(set(tmp)), key=lambda x: Counter(counter_pathway)[x])
for _ in used_list:
    print(_, Counter(counter_pathway)[
          _], tmp[_], '+'.join([_] + [k for k, ko in ko2p.items() if _ in list(ko)]))

fig = go.Figure()
fig.add_bar(x=[Counter(counter_pathway)[_] for _ in used_list],
            y=[tmp[_] for _ in used_list],
            orientation='h')
fig.write_html('./test.html')




## ladderance
from subprocess import check_call
cmd = "mkdir -p /mnt/home-backup/thliao/plancto/gene_GandL/itol_txt/usages"
check_call(cmd,shell=1)
for section_name in ['Anammox',
 'acquire nitrogen related',
 'ladderane synthesis',
 'key gene in Anaerobic biosynthesis of VB12',
 'oxidative stress',
 'Sulfur metabolism',
 'iron utilization amd hame',
 'bioenergetic',
 'oxidoreductase and ATPase',
 'mutualism']:
    f = f"./gene_GandL/itol_txt/{section_name.replace('/', ' or ')}.txt"
    cmd = f"cp '{f}' /mnt/home-backup/thliao/plancto/gene_GandL/itol_txt/usages/"
    check_call(cmd,shell=1)



# cmd (diamond)  abandon
# python3 ~/bin/batch_run/batch_any.py -i /mnt/home-backup/thliao/plancto/protein_annotations/tmp_faa -o /mnt/home-backup/thliao/plancto/special_genes/trihame_nitrite_reductase/annotated/ -s faa -ns tab -np 10 -cmd 'diamond blastp -q {infile} --db /mnt/home-backup/thliao/plancto/special_genes/trihame_nitrite_reductase/ref -f 6 --evalue 1e-50 --threads 3 -o {ofile} '
indir = "/mnt/home-backup/thliao/plancto/special_genes/trihame_nitrite_reductase/annotated/"
all_tabs = [_ for _ in glob(f"{indir}/*.tab") if basename(_).startswith('GCA')]
dfs = [pd.read_csv(_, sep='\t', header=None)
       for _ in all_tabs if os.path.getsize(_) != 0]
dfs = [df.loc[df[10] <= 1e-50, :] for df in dfs]
dfs = [df.sort_values(10, ascending=True).groupby(0).head(1)
       for df in dfs]
genome2info = defaultdict(list)
for df in tqdm(dfs):
    genome = "GCA_" + df.iloc[0, 0].split('_')[0].replace('v', '.')
    for _, row in df.iterrows():
        genome2info[genome].append(
            (row[0], row[1].split('|')[-1], float(row[10])))
post_filtered = filtration_part(genome2info, evalue=1e-50)

odir = f'{indir}/..'
if not exists(odir):
    os.makedirs(odir)
ofile_info = join(odir, "merged_info.tab")
ofile_binary = join(odir, "merged_binary.tab")
ofile_num = join(odir, "merged_num.tab")
final_df = pd.DataFrame.from_dict(post_filtered, orient='index')
bin_df = final_df.applymap(lambda x: 0 if pd.isna(x) else 1)
num_df = final_df.applymap(lambda x: 0 if pd.isna(x)
                           else len(str(x).split(',')))
final_df.to_csv(ofile_info, sep='\t', index=1)
bin_df.to_csv(ofile_binary, sep='\t', index=1)
num_df.to_csv(ofile_num, sep='\t', index=1)