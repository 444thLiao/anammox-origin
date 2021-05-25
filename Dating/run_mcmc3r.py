import multiprocessing as mp
import os
from os.path import exists

import numpy as np
import pandas as pd
from tqdm import tqdm

from dating_workflow.bin.dating_pro import run, modify

target_sets = ['C7', 'C13', 'C16', 'C21']
target_dir = '/mnt/home-backup/thliao/plancto/update_calibrations/dating_reorder/85g/nucl/clock2_diff_cal/repeat_nucl_C1_run1'

odir = "/mnt/home-backup/thliao/plancto/update_calibrations/dating_reorder/85g/nucl/mcmc3r"
for s in target_sets:
    for model in ['IR', 'AR']:
        if not exists(f"{odir}/{model}_{s}"):
            os.makedirs(f"{odir}/{model}_{s}")
        os.system(f'cp {target_dir}/04_mcmctree.ctl {odir}/{model}_{s}/')
        param = {'treefile': "/mnt/home-backup/thliao/plancto/update_calibrations/dating_reorder/cal_tree/85g_C1.newick".replace('C1', s),
                 'usedata': "2 ../../in.BV 1",
                 'clock': 2 if model == 'IR' else 3}
        text = modify(f'{odir}/{model}_{s}/04_mcmctree.ctl',
                      **param)
        with open(f'{odir}/{model}_{s}/04_mcmctree.ctl', 'w') as f1:
            f1.write(text)
        cmd = f"""R -e "setwd('{odir}/{model}_{s}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='04_mcmctree.ctl', betaf='beta.txt')" """
        os.system(cmd)

cmd = f"ln -s `realpath {target_dir}/in.BV` {odir}"
os.system(cmd)
    
from os.path import *
from glob import glob
cmds = []
for ctl in glob(f"{odir}/**/04_mcmctree.ctl",recursive=True):
    dname = dirname(ctl)
    if exists(join(dname,'beta.txt')):
        continue
    if not exists(join(dname,'mcmc.txt')):
        cmd = f"cd {dname}; /home-user/thliao/software/paml4.10/bin/mcmctree ./04_mcmctree.ctl > ./run.log"
        cmds.append(cmd)
        

from bin.multiple_sbatch import sbatch_all
sbatch_all(cmds,True,1)



def get_v(rout):
    outs = rout.split('\n')
    idx1, idx2 = outs.index('$logml'), outs.index('$se')
    logml, se = map(lambda x: float(x.split(' ')[-1].strip()),
                    (outs[idx1 + 1], outs[idx2 + 1]))
    return logml, se

collect_df = pd.DataFrame(columns=['calibration set', 'model', 'Log marginal (s. d)', 'BF'])
count = 0
for t in target_sets:
    cmd = f"""R -e "setwd('AR_{t}'); AR<- mcmc3r::stepping.stones(); AR " """
    AR = os.popen(cmd).read()
    AR_logml, AR_se = get_v(AR)
    cmd = f"""R -e "setwd('IR_{t}'); IR<- mcmc3r::stepping.stones(); IR " """
    IR = os.popen(cmd).read()
    IR_logml, IR_se = get_v(IR)
    c = np.array([AR_logml, IR_logml])
    BF = np.exp(c - np.max(c))
    collect_df.loc[count, :] = [t, 'AR', f'{AR_logml} ({AR_se})', BF[0]]
    collect_df.loc[count + 1, :] = [t, 'IR', f'{IR_logml} ({IR_se})', BF[1]]
    count += 2
        
# cmd = """ `which R` -e "setwd('{mcmc3r_dir}/mcmctree/{mn}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='mcmctree.ctl', betaf='beta.txt')" """



# cmd = """ `which R` -e "setwd('AR_C7'); AR <- mcmc3r::stepping.stones();setwd('../IR_C7'); IR <- mcmc3r::stepping.stones();mlnl <- c(AR$logml, IR$logml);BF <- exp(mlnl - max(mlnl)); BF" """
# mlnl

# # Bayes factors
# ( BF <- exp(mlnl - max(mlnl)) )
# BF

# # Posterior model probabilities
# ( Pr <- BF / sum(BF) )

# # or alternatively:
# mcmc3r::bayes.factors(IR, AR)
