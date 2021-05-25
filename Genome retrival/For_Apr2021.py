from collections import defaultdict
from glob import glob
import os
from os.path import *
from subprocess import check_call
from Bio import SeqIO
from tqdm import tqdm 
import pandas as pd
from ete3 import Tree
from typing_extensions import final
from api_tools.itol_func import to_binary_shape
from dating_workflow.step_script.quick_sampling import *

os.chdir("/mnt/home-backup/thliao/plancto/update_Apr2021")

intree = "trees/iqtree/over20p_bac120_formatted.newick"
st = Tree(intree,format=3)
name2n = {_.name:_ for _ in st.traverse()}
name2n['I27_S100'].up.remove_child(name2n['I27_S100'])
name2n['I21_S100'].up.remove_child(name2n['I21_S100'])

cmd = f"/home-user/thliao/script/TreeCluster/TreeCluster.py -i {intree} -tf argmax_clusters -t 1 > {intree}.cluster"
check_call(cmd,shell=1)
cluster_file = f"{intree}.cluster"
cluster2genomes = get_cluster(cluster_file)
g2cluster = {v:c for c,d in cluster2genomes.items() for v in d}
            
from dating_workflow.step_script.extract_cog25 import parse_annotation 
genome2cog25 = parse_annotation(expanduser('~/data/NCBI/modified_data/cog25_annotate/'), top_hit=True, evalue=1e-20)

target_nodes = ["I263_S100",   # anammox
                ]

must_in_gids = [_ for _ in open(expanduser("~/data/cyano/ref_tree_new/used_genomes.list")).read().split('\n') if _]
must_in_gids.extend(name2n['I4_S100'].get_leaf_names())  # verruco
# manually add more anammox 
must_in_gids.extend(["GCA_000315115.1","GCA_005524015.1","GCA_008933285.1",
                     "GCA_007860005.1",])

final_leaves = sampling(st, target_nodes, 
                        must_in=must_in_gids, 
                        node2cluster=g2cluster,
                        up_level=5, 
                        max_num_up_each_level=4,
                        max_num_down=90, 
                        genome2cog25=genome2cog25)
final_leaves.remove('GCA_905339135.1')
print(len(final_leaves))
text = to_binary_shape({_:['remained'] for _ in final_leaves},
                       {'remained':{'shape':'2','color':'#536DFE'}})
with open('./tmp.txt','w') as f1:
    f1.write(text)

odir = './dating/id_list/'
if not exists(odir):
    os.system(f"mkdir -p {odir}")
with open(f"{odir}/{len(final_leaves)}g_ids.list",'w') as f1:
    f1.write('\n'.join(final_leaves))

name = f"{len(final_leaves)}g"
gl = f"{odir}/{name}_ids.list"

cmd = f"""
python3 ~/script/evol_tk/dating_workflow/step_script/extract_bac120.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/bac120_annotate -o ./bac120_extract/{name} -evalue 1e-50 -gl {gl}
python3 ~/bin/batch_run/batch_mafft.py -i  ./bac120_extract/{name} -s faa -o ./bac120_extract/{name} -gl {gl} -m mafft -f ;
python3 ~/bin/batch_run/batch_trimal.py -i  ./bac120_extract/{name} -s aln -o ./bac120_extract/{name}
python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i ./bac120_extract/{name} -o ./trees/concat/{name}.trimal -s trimal -gl {gl} -ct partition -no_graph

mkdir -p ./trees/fasttree; iqtree -nt 30 -m LG -redo -pre ./trees/fasttree/{name}_guide -s ./trees/concat/{name}.trimal
mkdir -p ./trees/iqtree; iqtree -nt 30 -wbtl -bb 1000 -m LG+C20+F+G -redo -ft ./trees/fasttree/{name}_guide.treefile -s ./trees/concat/{name}.trimal -pre ./trees/iqtree/{name}_bac120_complex
"""
cmd = f"""
format_newick.py set1 -i ./trees/iqtree/{name}_bac120_complex.contree -r "GCA_000019965.1,GCA_016192525.1" -o ./trees/final/{name}_template.newick -f 0
nw_prune ./trees/final/{name}_template.newick 'GCA_905339135.1'  > ./trees/final/{len(final_leaves)-1}g_template.newick

"""
check_call(cmd,shell=1)

r = open('/mnt/home-backup/thliao/plancto/update_calibrations/dating_reorder/calibration_files/cal_C7.txt').read()
nr = r.replace('GCA_003576915.1','GCA_000019665.1')
with open('./dating/cal_C7.txt','w') as f1:
    f1.write(nr)
cmd = f"""format_newick.py add-cal -i ./trees/final/{len(final_leaves)-1}g_template.newick -o ./dating/cal_tree/{len(final_leaves)-1}g_C7.newick -c ./dating/cal_C7.txt -f 3"""
check_call(cmd,shell=1)



st = f"./trees/final/{len(final_leaves)-1}g_template.newick"
t = Tree(st,3)
gids = t.get_leaf_names()
odir = './dating/id_list/'
gl = f"{odir}/{len(gids)}_ids.list"
with open(gl,'w') as f1:
    f1.write('\n'.join(gids))
##
cmd = f"""
python3 /home-user/thliao/script/evol_tk/dating_workflow/step_script/extract_cog25.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate -pd /mnt/home-backup/thliao/NCBI/modified_data/prokka_o -o ./cog25_single/nucl_e20 -evalue 1e-20 -gl {gl} -ot nucl
"""
check_call(cmd,shell=1)
    
cmd = f"""
python3 ~/bin/batch_run/batch_mafft.py -i ./cog25_single/nucl_e20 -s ffn -o ./cog25_single/{len(gids)}g_nucl -f -m mafft -gl {gl}
python3 ~/bin/batch_run/batch_trimal.py -i ./cog25_single/{len(gids)}g_nucl -o ./cog25_single/{len(gids)}g_nucl;
python3 ~/script/evol_tk/dating_workflow/bin/concat_aln.py -i ./cog25_single/{len(gids)}g_nucl -ct phy -gl {gl} -o ./trees/phy_files/{len(gids)}g_cog25_nucl.trimal -s trimal -no_graph"""
check_call(cmd,shell=1)


f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i ./trees/phy_files/{len(gids)}g_cog25_nucl.phy -it ./dating/cal_tree/{len(gids)}g_C7.newick -o ./dating/{len(gids)}g/nucl/C7 -p 1 -rg '1 100 1' -sf 30 -c 2  &"

cmd = f"python3 ~/script/evol_tk/dating_workflow/figtree2itol.py -i ./trees/final/{len(gids)}g_template.newick -i2 ./dating/{len(gids)}g/nucl/C7/mcmc_for/FigTree.tre -o ./dating/{len(gids)}g/nucl/C7/posterior.newick"
check_call(cmd,shell=1)