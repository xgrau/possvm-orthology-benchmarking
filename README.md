# *Possvm* accuracy benchmarks

This is a companion repository to [*Possvm*](https://github.com/xgrau/possvm-orthology), which contains scripts and data to reproduce the benchmarking analyses done for the *Possvm* manuscript ([Grau-Bové and Sebé-Pedrós, bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.05.03.442399v1)).

## TALE homeobox classification

Test the accuracy of *Possvm* using manually curated homeobox families from HomeoDB. 

From the `**homeobox-test**` folder:

1. Get TALE, PRD and ANTP sequences from [HomeoDB](http://homeodb.zoo.ox.ac.uk/) (as of 7th Feb 2021).

2. Blast (diamond).

```bash
# to orthobench dataset (17 metazoan species)
# bash s01_get_trees-diamond.sh seed_tale.fasta tale ../orthobench-test/proteomes/
# bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTP ../orthobench-test/proteomes/
# bash s01_get_trees-diamond.sh seed_PRD.fasta PRD ../orthobench-test/proteomes/
# to bilateria+cnidaria dataset
bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTP proteomes/
bash s01_get_trees-diamond.sh seed_PRD.fasta  PRD  proteomes/
bash s01_get_trees-diamond.sh seed_TALE.fasta TALE proteomes/
```

3. Run *Possvm*.

```bash
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom -itermidroot 10
possvm -i results_trees/PRD.genes.iqtree.treefile -p PRD.possom -itermidroot 10
possvm -i results_trees/TALE.genes.iqtree.treefile -p TALE.possom -itermidroot 10
```

4. Evaluate using classification from blast to HomeoDB.

```bash
Rscript s02_evaluate_homeodb.R
```

5. Plot and annotate trees?

```bash
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom -itermidroot 10 -r reference_ANTP.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
possvm -i results_trees/PRD.genes.iqtree.treefile -p PRD.possom -itermidroot 10 -r reference_PRD.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
possvm -i results_trees/TALE.genes.iqtree.treefile -p TALE.possom -itermidroot 10 -r reference_TALE.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
```

6. Downsample ANTP trees, by downsampling dataset (always include Homo so that we have at least one reference species to evaluate, though):

```bash
# create 20 versions of the original dataset downsampled to 10%, 20%, ... 70% of the species
Rscript s10_downsample_fasta.R
# launch trees
for i in results_trees/downsampling_alignment/downsample_0.*.genes.lt.fasta ; do qsub -N tree-$(basename ${i}) -pe smp 4 qsub_iqtree.sh ${i} 4 ; done
# call OGs
for i in results_trees/downsampling_alignment/downsample*treefile ; do possvm -i $i -p $(basename ${i}) -itermidroot 10 -skipprint ; done
# evaluate precision and recall
# modify s05_evaluate_permutations_ref.R
```

7. Evaluate effect of tree accuracy in ANTP classification, by permutating tip labels in the original tree

```bash
# permutations
Rscript s20_randomise_tips_all.R
# call OGs
for i in results_trees/downsampling/downsample*newick ; do possvm -i $i -p $(basename ${i}) -itermidroot 10 -skipprint ; done
# evaluate precision and recall
Rscript s21_evaluate_permutations_all.R
```

## Orthobench testing

Test the accuracy of *Possvm* using manually curated orthologs from the [Orthobench 2.0 repository](https://github.com/davidemms/Open_Orthobench) (see [Emms et al. GBE 2020](https://academic.oup.com/gbe/article/12/12/2258/5918455)).

From the `**orthobench-test**` folder:

1. Get Orthobench reference groups (`refOG`) from the [Open_Orthobench Github](https://github.com/davidemms/Open_Orthobench/tree/master/BENCHMARKS); total: 70 RefOGs), as well as raw and tight trees (ignore intermediate trees).

```bash
# get Orthobench repository
git clone git@github.com:davidemms/Open_Orthobench.git

# refOGs assignments:
ls <path>/Open_Orthobench/BENCHMARKS/RefOGs/

# tight and raw trees:
ls <path>/Open_Orthobench/Supporting_Data/Data_for_RefOGs/trees
ls <path>/Open_Orthobench/Supporting_Data/Data_for_RefOGs/trees_tight

# # HMM profiles in this folder
# ls <path>/Open_Orthobench/Supporting_Data/Additional_Files/hmm_profiles/
# # copy the HMM files into the `hmm_profiles_strict` and `hmm_profiles_weak` folders in the present directory
```

2. NOTE: Fix duplicated assignment of `FBpp0309618` to the `RefOG021` and `RefOG068` clusters (remove from latter), format into `refOGs.csv`.

3. Get primary transcripts from 17 metazoan species from the [Orthobench repository](https://github.com/davidemms/Open_Orthobench/tree/master/Supporting_Data/Additional_Files/proteomes/primary_transcripts), and HMM profiles for each refOG:

```bash
# primary transcripts for each species, in this folder:
ls <path>/Open_Orthobench/Supporting_Data/Additional_Files/proteomes/primary_transcripts
# copy the fasta files into the `proteomes` folder in in the present directory and concatenate them...
```

4. Run homology searches, MSAs and phylogenies:

```bash
# run from the present directory
bash s01_get_trees-diamond.sh
```

5. Run *Possvm* as follows:

```bash
# find OGs in each tree with possom:
for i in orthobench_trees/tight/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done
for i in orthobench_trees/raw/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done
for i in results_orthology/*.treefile ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done
```

6. Calculate accuracy relative to `refOGs.csv`:

```bash
Rscript s02_evaluate_orthobench.R
```
