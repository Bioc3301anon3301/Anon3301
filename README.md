# Anon3301
# First step required is to split the library and count sequences
# The batch script below splits the Read 1 which is the forward read from Illumina,the Index file is the barcode file which contains unique barcodes per sample
# -i denotes input -m is map -o is output -b is barcode reference
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# splitting libraries
echo "splitting libraries"
time split_libraries_fastq.py --rev_comp_barcode --rev_comp_mapping_barcodes 
-i Read1.fastq.gz -b Index.fastq.gz -o slout -m map.tsv
# counting sequences
echo "Counting sequences"
time count_seqs.py -i slout/seqs.fna
# deactivating environment
source deactivate


# Next we choose the OTUs using SILVA closed reference
#!/bin/bash          
#SBATCH -t 1:00:00                                                          
#SBATCH -n 16                                                               
#SBATCH -p short                                          
#Load modules to run python                                                            

module load eb

module load Miniconda2

source deactivate

# loading virtualenv                                                                   

echo "loading virtualenv"

source activate qiime1

# setting temporary directory                                                          

export TMPDIR=~/qiime_tmp

# picking OTUs                                                                         

echo "Picking OTUs with closed reference"

pick_closed_reference_otus.py
-i seqs.fna -o otus -r SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t SILVA_128_QIIME_release/taxonomy/16S_only/97/majority_taxonomy_all_levels.txt -a -O 16 -o pickedOTUs

# deactivating environment
source deactivate


# Then we perform a core diversity analyses which provides alpha and beta diversity

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#SBATCH --reservation=cbmucl
#Load modules to run python
module load eb
module load Miniconda2
source deactivate
# loading virtualenv
echo "loading virtualenv"
source activate qiime1
# setting temporary directory
export TMPDIR=~/qiime_tmp

# Running core diversity analyses, 30k sampling depth, 97% silva reference tree"

core_diversity_analyses.py --recover_from_failure -o ./cdout2018 -i pickedOTUs/OTU_table.biom -m map.tsv -t ~/silva/trees/97/97_otus.tre -e 30000 -a -O 24

# deactivating environment
source deactivate


# After generating a core diversity analyses of 2017 samples, a similar workflow is performed for 2016.

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# splitting libraries
echo "splitting libraries"
time split_libraries_fastq.py --barcode_type 12 -i Read1_YearAgo.fastq.gz -b IndexYearAgo.fastq.gz -o 2017slout -m /2017_03_smb/map_YearAgo.tsv
# counting sequences
echo "Counting sequences"
time count_seqs.py -i 2017slout/seqs.fna
# deactivating environment
source deactivate


# Similarly, SILVA closed referece OTUs were picked for the 2016 samples

#!/bin/bash
#SBATCH -t 1:00:00                                                          
#SBATCH -n 24                                               
#SBATCH -p short                                               
#Load modules to run python                                                              

module load eb

module load Miniconda2

source deactivate

# loading virtualenv                                                                   

echo "loading virtualenv"

source activate qiime1

# setting temporary directory                                                          

export TMPDIR=~/qiime_tmp

# picking OTUs                                                                         

echo "Picking OTUs with closed reference"

pick_closed_reference_otus.py -i ./slout2017/seqs.fna -o otus2017 -r silva/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t silva/taxonomy/16S_only/97/majority_taxonomy_all_levels.txt -a -O 16 -o pickedOTUs2017


# CDA is performed on 2016 samples, 2016 samples are tagged as 2017 as the experiments were run in 2017

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules to run python
module load eb
module load Miniconda2
source deactivate
# loading virtualenv
echo "loading virtualenv"
source activate qiime1
# setting temporary directory
export TMPDIR=~/qiime_tmp

# Running core diversity analyses, 30k sampling depth, 97% silva reference tree

time core_diversity_analyses.py --recover_from_failure -o ./cdout2017 -i ./pickedOTUs2017/otu_table.biom -m ~/2017_03_smb/map_YearAgo.tsv -t ~/silva/trees/97/97_otus.tre -e 30000 -a -O 24

#  2016 and 2017 maps were meregd for downstream analyses

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# mergemaps
merge_mapping_files.py -m 2018_02_smb/map.tsv,2017_03_smb/map_YearAgo.tsv -o mergedmap38
# deactivating
source deactivate

# OTU tables are merged for downstream analyses

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# merge OTU tables
echo merge otu tables
merge_otu_tables.py
-i pickedOTUs/otu_table.biom, pickedOTUs2017/otu_table.biom -o mergedOTU
#deactivate
source deactivate

# CDA on merged samples

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules to run python
module load eb
module load Miniconda2
source deactivate
# loading virtualenv
echo "loading virtualenv"
source activate qiime1
# setting temporary directory
export TMPDIR=~/qiime_tmp
# Running core diversity analyses, 30k sampling depth, 97% silva reference tree 

time core_diversity_analyses.py --recover_from_failure -o ./cdoutmergedfixed38 -i mergedOTU_table.biom -m mergedmap38 -t ~/silva/trees/97/97_otus.tre -e 30000 -a -O 24
# deactivating environment
source deactivate

# Obtain a distance matrix based on category Year

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# distance matrix year

distance_matrix_from_mapping.py -i mergedmap38 -c Year -o Yeardm38

#deactivating
source deactivate

# Run Mantel Statistics on unweighted Beta diversity

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# comparing distance mantel
compare_distance_matrices.py --method=mantel -i cdoutmergedfixed38/bdiv_even30000/unweighted_unifrac_dm.txt,Yeardm38 -o unweightedyearmantel_out38 -n 999
# deactivating
source deactivate

# Run Mantel Statistics on weighted Beta diversity

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# comparing distance mantel
compare_distance_matrices.py --method=mantel -i cdoutmergedfixed38/bdiv_even30000/weighted_unifrac_dm.txt,Yeardm38 -o weightedyearmantel_out38 -n 999
# deactivating
source deactivate


# Run ANOSIM statistics on unweighted unifrac

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# compare distance matrices
echo compare distance matrices
compare_categories.py --method anosim -i cdoutmergedfixed38/bdiv_even30000/unweighted_unifrac_dm.txt -m mergedmap38 -c Year -o unweightedanosim_out38 -n 99
# deactivating
source deactivate

# Run ANOSIM statistics on weighted unifrac

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
# compare distance matrices
echo compare distance matrices
compare_categories.py --method anosim -i cdoutmergedfixed38/bdiv_even30000/unweighted_unifrac_dm.txt -m mergedmap38 -c Year -o unweightedanosim_out38 -n 99
# deactivating
source deactivate

# Compute core microbiome of each individual year

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#SBATCH --reservation=cbmucl
#Load modules to run python
module load eb
module load Miniconda2
source deactivate
# loading virtualenv
echo "loading virtualenv"
source activate qiime1
# setting temporary directory
export TMPDIR=~/qiime_tmp
# compute core microbiome
echo "compute core microbiome"
compute_core_microbiome.py
-i ~pickedOTUs/otu_table.biom -o otu_table_core2018 --min_fraction_for_core 0.1
# deactivate
echo "deactivate"
source deactivate

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
#SBATCH --reservation=cbmucl
#Load modules to run python
module load eb
module load Miniconda2
source deactivate
# loading virtualenv
echo "loading virtualenv"
source activate qiime1
# setting temporary directory
export TMPDIR=~/qiime_tmp
# compute core microbiome
echo "compute core microbiome"
compute_core_microbiome.py
-i ~pickedOTUs2017/otu_table.biom -o otu_table_core2017 --min_fraction_for_core 0.1
# deactivate
echo "deactivate"
source deactivate

# Summarize taxa information of 2017 samples

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
echo summarize taxonomy plots
#summarise taxonomy plots based on taxonomic data generated
summarize_taxa.py -i pickedOTUs/otu_table.biom -m 2018_02_smb/map.tsv -o taxasummary2018
#deactivate
source deactivate

# summarize taxa information of 2016 samples

#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 24
#SBATCH -p short
##SBATCH --reservation=cbmucl
#Load modules
module load eb
module load Miniconda2
# setting temporary directory
export TMPDIR=~/qiime_tmp
# loading virtualenv
source activate qiime1
echo summarize taxonomy plots
#summarise taxonomy plots based on taxonomic data generated
summarize_taxa.py -i pickedOTUs2017/otu_table.biom -m 2017_03_smb/map_YearAgo.tsv -o taxasummary2017
#deactivate
source deactivate
