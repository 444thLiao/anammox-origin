# anammox-origin
Dated the origin time of anammox bacteria with latest genomical sequences from NCBI

Molecuar dating analysis of the origin time of anammox bacteria includes multiple parts and a raw workflow for dating a focal lineage with all published genomes.
It includes but not only
1. download all genomes available (not limited at NCBI)
2. build species tree with markers genes and gene trees and following reconciliation
3. annotation with multiple databases for following comparative analysis
4. molecular dating analysis based on MCMCTree impletement in PAML and *Bayes Factor calculation*


# Contents of each directory

* Comparative analysis 
  * Scripts and full results of genome comparative analysis and enrichment analsyis
* Dating analysis
  * template for mcmctree
  * script for mcmc3r 
  * Log and required files for dating
  * Cog25 genes of corresponding genomes
* Genome retrival
  * script for genome retrival and list of accession id of used genomes
* Phylogenetic analysis
  * 16S tree
  * gene trees (hzsCBA)
  * species tree
* Reconcilications
  * configuration file and input data for GenRax
* Visualizations
  * script generating itol annotation files 

Above steps mostly are depend on the another toolkit ([evol-tk](https://github.com/444thLiao/evol_tk)) 
If you want to utilize scripts involved in above steps, maybe you should import functions or run some scripts from this repo.


# Potential question

The calibration used within the scripts might have little bits different to the paper said. It is because the ordering has changed once time. Such as, C7 in the scipt represent the C1 in the paper.
More you could see, the file `./Dating/calibrations_sets_new.xlsx` (sheet 'rename')


# Publication
Not yet

# Contact Us
If you have any questions or suggestions, you are welcome to contact us via email: l0404th@gmail.com