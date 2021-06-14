# *Possvm* accuracy benchmarks

This is a companion repository to [*Possvm*](https://github.com/xgrau/possvm-orthology), and contains scripts and data to reproduce the benchmarking analyses done for the *Possvm* manuscript ([Grau-Bové and Sebé-Pedrós, bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.05.03.442399v1)).

You can find the *Possvm* source code, manual and installation instructions in the [main repository for this project](https://github.com/xgrau/possvm-orthology).

## Data sources

If you use the datasets provided in this repository, cite the original sources:

* [HomeoDB](http://homeodb.zoo.ox.ac.uk/) database (v2): **[Zhong et al. Evolution & Development, 2011](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1525-142X.2011.00513.x)**.
* [Orthobench](https://github.com/davidemms/Open_Orthobench) repository (v2.0): **[Emms et al. GBE 2020](https://academic.oup.com/gbe/article/12/12/2258/5918455)**.

If you use *Possvm*, please cite the following papers:

* *Possvm*: **[Grau-Bové and Sebé-Pedrós, bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.05.03.442399v1)**.
* *ETE* toolkit: **[Huerta-Cepas *et al.* Molecular Biology and Evolution 2016](https://academic.oup.com/mbe/article/33/6/1635/2579822)**.
* Species overlap algorithm: **[Huerta-Cepas *et al.* Genome Biolgy 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)**.
* *MCL* clustering: **[Enright *et al.* Nucleic Acids Research 2002](https://pubmed.ncbi.nlm.nih.gov/11917018/)**.

## Homeobox classification

Test the accuracy of *Possvm* using manually curated homeobox families from HomeoDB.

From the **`homeobox-test`** folder:

1. Get ANTP sequences from [HomeoDB](http://homeodb.zoo.ox.ac.uk/) (as of 7th Feb 2021).

2. Blast (diamond).

```bash
# to bilateria+cnidaria dataset
bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTP proteomes/
# UNUSED:
# bash s01_get_trees-diamond.sh seed_PRD.fasta  PRD  proteomes/
# bash s01_get_trees-diamond.sh seed_TALE.fasta TALE proteomes/

# UNUSED:
# to orthobench dataset (17 metazoan species)
# bash s01_get_trees-diamond.sh seed_tale.fasta tale ../orthobench-test/proteomes/
# bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTP ../orthobench-test/proteomes/
# bash s01_get_trees-diamond.sh seed_PRD.fasta PRD ../orthobench-test/proteomes/
```

3. Run *Possvm*.

```bash
# base run with iterative rooting
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom -itermidroot 10
# run with iterative rooting, excluding cnidarians from the orthology graph
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_nocni -itermidroot 10 --outgroup outgroups.txt
# run with iterative rooting and LPA clustering
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_lpa -itermidroot 10 -method lpa
# run with iterative rooting and Louvain clustering
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_lou -itermidroot 10 -method louvain
# run with iterative rooting and Clauset-Newman-Moore greedy modularity maximization
python ../scripts/possvm-greedy.py -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_gre -itermidroot 10
# run with midpoint rooting
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_mid

# possvm -i results_trees/PRD.genes.iqtree.treefile -p PRD.possom -itermidroot 10
# possvm -i results_trees/TALE.genes.iqtree.treefile -p TALE.possom -itermidroot 10

```

4. Evaluate using classification from blast to HomeoDB.

```bash
# using one-to-one orthogroup assignments (i.e. only best possvm OG is considered)
Rscript s02_evaluate_homeodb_one2one.R
Rscript s02_evaluate_homeodb_one2one_mid.R
Rscript s02_evaluate_homeodb_one2one_lpa.R
Rscript s02_evaluate_homeodb_one2one_louvain.R
Rscript s02_evaluate_homeodb_one2one_greedymod.R
# using one-to-many orthogroup assignments (i.e. all possvm OGs are considered)
Rscript s02_evaluate_homeodb_one2many.R
Rscript s02_evaluate_homeodb_one2many_mid.R
```

5. Annotate trees with human genes:

```bash
# using iterative rooting
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom -itermidroot 10 -r reference_ANTP.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
# possvm -i results_trees/PRD.genes.iqtree.treefile -p PRD.possom -itermidroot 10 -r reference_PRD.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
# possvm -i results_trees/TALE.genes.iqtree.treefile -p TALE.possom -itermidroot 10 -r reference_TALE.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50

# using midpoint rooting
possvm -i results_trees/ANTP.genes.iqtree.treefile -p ANTP.possom_mid -r reference_ANTP.type.csv -refsps Hsap -o results_annotation -printallpairs -min_support_transfer 50
```

6. Downsample ANTP trees, by downsampling dataset (always include Homo so that we have at least one reference species to evaluate, though):

```bash
# UNUSED
# create 20 versions of the original dataset downsampled to 10%, 20%, ... 70% of the species
# Rscript s10_downsample_original_alignment.R
# # launch trees
# for i in results_trees/downsampling_alignment/downsample_0.*.genes.lt.fasta ; do qsub -N tree-$(basename ${i}) -pe smp 4 qsub_iqtree.sh ${i} 4 ; done
# # call OGs
# for i in results_trees/downsampling_alignment/downsample*treefile ; do possvm -i $i -p $(basename ${i}) -itermidroot 10 -skipprint ; done
# # evaluate precision and recall
# Rscript s11_evaluate_downsampling.R
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

## Orthobench

Test the accuracy of *Possvm* using manually curated orthologs from the [Orthobench 2.0 repository](https://github.com/davidemms/Open_Orthobench) (see [Emms et al. GBE 2020](https://academic.oup.com/gbe/article/12/12/2258/5918455)).

From the **`orthobench-test`** folder:

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

4. Run *Possvm* as follows:

```bash
# find OGs in each tree with possom, using the original trees from orthobench:
for i in orthobench_trees/raw/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done
for i in orthobench_trees/raw/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_iter -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done


# UNUSED
# for i in orthobench_trees/tight/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done # UNUSED
```

5. Calculate accuracy relative to `refOGs.csv`:

```bash
# using one-to-one orthogroup assignments (i.e. only best possvm OG is considered)
Rscript s02_evaluate_orthobench_one2one.R
# using one-to-many orthogroup assignments (i.e. all possvm OGs are considered)
Rscript s02_evaluate_orthobench_one2many.R
```

Run homology searches, MSAs and new phylogenies:

```bash
# # run from the present directory
# bash s01_get_trees-diamond.sh

# # find OGs in each tree with possom, using new expanded trees (necessary to assess the effect of iterative rooting?)
# cp results_trees/*.treefile results_orthology/
# for i in results_orthology/*.treefile ; do possvm -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." ; done
# for i in results_orthology/*.treefile ; do possvm -i $i -p $(basename ${i%%.*}).possom_iter -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done

# # evaluate iterative rooting accuracy:
# # DEPRECATED
# Rscript s03_evaluate_orthobench_one2one_iter.R
# Rscript s04_compare_midroot_iterroot.R
```

### Effect of iterative rooting

Based on raw Orthobench trees, let's increase the length of random internal branches to evaluate the effect of the iterative tree rooting procedure:

1. Create collection of inflated trees (20 per original tree):

```bash
s05_create_inflated_trees.R
```

2. Run *Possvm* with and without iterative rooting:

```bash
# inflation in one branch
for i in results_rooting_inflation_one/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_mid -ogprefix "$(basename ${i%%.*})." ; done
for i in results_rooting_inflation_one/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_ite -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done

# low inflation (bound between 5x and 20x, for 5% of edges)
for i in results_rooting_inflation/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_mid -ogprefix "$(basename ${i%%.*})." ; done
for i in results_rooting_inflation/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_ite -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done

# high inflation (bound between 5x and 50x, for 10% of edges)
for i in results_rooting_inflation_high/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_mid -ogprefix "$(basename ${i%%.*})." ; done
for i in results_rooting_inflation_high/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_ite -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done

# higher inflation (bound between 5x and 50x, for 20% of edges)
for i in results_rooting_inflation_higher/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_mid -ogprefix "$(basename ${i%%.*})." ; done
for i in results_rooting_inflation_higher/*.tre ; do possvm -i $i -p  $(basename ${i%%.*}).possom_ite -ogprefix "$(basename ${i%%.*})." -itermidroot 10 ; done

# TODO: decide whether to use high or low inflation!
# Probably high is better
```

3. Evaluate difference:

```bash
Rscript s06_evaluate_iterroot_one2one.R
# TODO: decide whether to use high or low inflation!
```

## Alternative methods

### Phylome-style

Let's define orthogroups from a tree using naive approaches based on species overlap, and compare these orthogroups with those defined by *Possvm*.

#### Phylome-style in Orthobench

From the **`orthobench-test`** folder:

* OGs based on the list of genes descending from a given duplication node (absurd):

```bash
# UNUSED
# for i in orthobench_trees/raw/*.tre ; do python ../scripts/phylome-dups.py -i $i -p  $(basename ${i%%.*}).dups -ogprefix "$(basename ${i%%.*})." ; done
```

* OGs based on connected components starting from a focus species (human), using Orthobench trees:

```bash
for i in orthobench_trees/raw/*.tre ; do python ../scripts/phylome-ccs.py -i $i -p  $(basename ${i%%.*}).ccs -ogprefix "$(basename ${i%%.*})." ; done
Rscript s10_evaluate_orthobench_one2one_CCfocus.R
```

* OGs based on connected components, using PhylomeDB trees from a focus species (human) — *this seems fair, but it's not a great idea: many reasons for sequence dropout may artifically reduce recall*.

```bash
# DEPRECATED
# get data from http://www.phylomedb.org/phylome_514 (meta human phylome with 66 species, including 15 species also present in Orthobench)
cat ../phylome-data/phylome-0514-human/proteomes/*.fa > ../phylome-data/phylome-0514-human/proteomes_phylome.fasta
diamond makedb --in ../phylome-data/phylome-0514-human/proteomes_phylome.fasta -d ../phylome-data/phylome-0514-human/proteomes_phylome.fasta

# find matching sequences in phylomeDB
mkdir -p results_phylomedb/
for i in $(cut -f1 refOGs.csv | sort -u ) ; do
diamond blastp --more-sensitive -d ../phylome-data/phylome-0514-human/proteomes_phylome.fasta -q results_searches/${i}.seed.fasta -o results_phylomedb/${i}.diamond_phylome.csv -p 4
done

# do it with all seqs
diamond blastp -q ../phylome-data/phylome-0514-human/proteomes_phylome.fasta -d proteomes/all_proteomes.fa -o results_phylomedb/dictionary_diamond_phylome.csv -p 10
awk '$3 == 100 { print $1,$2 }' results_phylomedb/dictionary_diamond_phylome.csv | tr ' ' '\t' > results_phylomedb/dictionary_diamond_phylome_filt.csv

# get ortholog pairs from phylomeDB
for i in $(cut -f1 refOGs.csv | sort -u ) ; do
fgrep -f <(grep -w ${i} refOGs.csv | cut -f2) results_phylomedb/dictionary_diamond_phylome_filt.csv | grep "_HUMAN" | cut -f1 | sort -u > results_phylomedb/${i}.diamond_phylome.human.txt
fgrep -f results_phylomedb/${i}.diamond_phylome.human.txt ../phylome-data/phylome-0514-human/orthologs.txt > results_phylomedb/${i}.diamond_phylome.human_orthologs.csv
done

# get orthogroups and evaluate estimates
Rscript s20_orthogroups_from_phylomeDB.R
```

#### Phylome-style in HomeoDB

From the **`homeobox-test`** folder:

* OGs based on connected components starting from a focus species (human), using Orthobench trees:

```bash
i=results_trees/ANTP.genes.iqtree.treefile
python ../scripts/phylome-ccs-focus.py -i $i -p  $(basename ${i%%.*}).ccs -ogprefix ccs
Rscript s10_evaluate_orthobench_one2one_CCfocus.R
```

### BranchClust

1. Download:

```bash
#wget http://www.bioinformatics.org/branchclust/BranchClust_1-01.txt && mv BranchClust_1-01.txt ../scripts/BranchClust_1.01.pl
wget http://www.bioinformatics.org/branchclust/BranchClust_1-00.txt && mv BranchClust_1-00.txt ../scripts/BranchClust_1.00.pl
wget http://www.bioinformatics.org/branchclust/BranchClust_Tutorial.tgz && mv BranchClust_Tutorial.tgz ../scripts/
tar xvfz ../scripts/BranchClust_Tutorial.tgz && mv BranchClust_Tutorial/ ../scripts
wget http://bioinformatics.org/branchclust/BranchClust_all.tgz && tar xvfz BranchClust_all.tgz && mv BranchClust_all.tgz BranchClust_all ../scripts
```

#### BranchClust in Orthobench

From the **`orthobench-test`** folder:

1. Prepare taxa dictionary files:

```bash
mkdir -p results_branchclust/
for i in orthobench_trees/raw/RefOG0*.ortholog_groups.csv ; do awk 'NR>1' $i | cut -f1 -d '_' ; done | sort -u | awk '{ print $1" | " $1"_.*" }' > results_branchclust/gi_numbers.out
```

2. Run BranchClust at various possible `<MANY>` values: 50% of the datset, 60%, 70%, and 80% (to assess tradeoff between precision and recall):

```bash
cd results_branchclust/
for i in $(cut -f1 ../refOGs.csv | sort -u ) ; do
# let's define ntax as 60% of the species in tree (recommended values: 50-80%)
# ntax=$(cat ../orthobench_trees/raw/${i}.tre| tr ',' '\n' | tr -d '()' |grep "^[A-Z]" | cut -f1 -d '_'| sort -u| wc -l | awk 'n=($1 * 0.6) { print int(n) }')
# 50%
perl ../../scripts/BranchClust_1.00.pl ../orthobench_trees/raw/${i}.tre 9
mv clusters.out ${i}.bc_clusters_50.out
mv families.list ${i}.bc_families_50.out
grep -A 1 ' CLUSTER ' ${i}.bc_clusters_50.out | grep -v " CLUSTER " | grep -v "\-\-" | awk '{ for(i=1; i<=NF; i++) { print $i"\t'${i}'."NR }}' > ${i}.bc_clusters_50.csv
# 60%
perl ../../scripts/BranchClust_1.00.pl ../orthobench_trees/raw/${i}.tre 10
mv clusters.out ${i}.bc_clusters_60.out
mv families.list ${i}.bc_families_60.out
grep -A 1 ' CLUSTER ' ${i}.bc_clusters_60.out | grep -v " CLUSTER " | grep -v "\-\-" | awk '{ for(i=1; i<=NF; i++) { print $i"\t'${i}'."NR }}' > ${i}.bc_clusters_60.csv
# 70%
perl ../../scripts/BranchClust_1.00.pl ../orthobench_trees/raw/${i}.tre 12
mv clusters.out ${i}.bc_clusters_70.out
mv families.list ${i}.bc_families_70.out
grep -A 1 ' CLUSTER ' ${i}.bc_clusters_70.out | grep -v " CLUSTER " | grep -v "\-\-" | awk '{ for(i=1; i<=NF; i++) { print $i"\t'${i}'."NR }}' > ${i}.bc_clusters_70.csv
# 80%
perl ../../scripts/BranchClust_1.00.pl ../orthobench_trees/raw/${i}.tre 14
mv clusters.out ${i}.bc_clusters_80.out
mv families.list ${i}.bc_families_80.out
grep -A 1 ' CLUSTER ' ${i}.bc_clusters_80.out | grep -v " CLUSTER " | grep -v "\-\-" | awk '{ for(i=1; i<=NF; i++) { print $i"\t'${i}'."NR }}' > ${i}.bc_clusters_80.csv
done
```

3. Test BranchClust:

```bash
Rscript s21_evaluate_branchclust.R
```

#### BranchClust in HomeoDB

From the **`homeobox-test`** folder:

1. Prepare taxa dictionary files:

```bash
mkdir -p results_branchclust/
for i in orthobench_trees/raw/RefOG0*.ortholog_groups.csv ; do awk 'NR>1' $i | cut -f1 -d '_' ; done | sort -u | awk '{ print $1" | " $1"_.*" }' > results_branchclust/gi_numbers.out
```

2. Run BranchClust at various possible `<MANY>` values: 50% of the datset, 60%, 70%, and 80% (to assess tradeoff between precision and recall):

```bash
cd results_branchclust/
i=../results_trees/ANTP.genes.iqtree.treefile
# let's define ntax as 60% of the species in tree (recommended values: 50-80%)
for m in 50 60 70 80 ; do
ntax=$(cat ${i} | tr ',' '\n' | tr -d '()' |grep "^[A-Z]" | cut -f1 -d '_'| sort -u| wc -l | awk 'n=($1 * "'$m'" / 100) { print int(n) }')
perl ../../scripts/BranchClust_1.00.pl ${i} ${ntax}
mv clusters.out ANTP.bc_clusters_${m}.out
mv families.list ANTP.bc_families_${m}.out
grep -A 1 ' CLUSTER ' ANTP.bc_clusters_${m}.out | grep -v " CLUSTER " | grep -v "\-\-" | awk '{ for(i=1; i<=NF; i++) { print $i"\tOG."NR }}' > ANTP.bc_clusters_${m}.csv
done
```

3. Test BranchClust:

```bash
Rscript s21_evaluate_branchclust.R
```

