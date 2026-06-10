#!/bin/bash -l

# Name: assembly.sh
# Purpuse: perform bacterial assembly to multiple isolates
# Author: Savvas Paragkamian, Mara Plakogiannaki
# Date: 01/04/2025

######################### User Input ######################
time_start=`date +%s`

#cd /media/sarlab/DATA/
#md5sum -c BGI_result/md5.txt 
# if all is ok proceed to uncompress

# sparagkamian@sarrislab:/media/sarlab/DATA/biosolutions_project/genomes$ ls -d -1 */ > ~/Documents/biosolutions/assemblies/sequences_directories.txt 
# other way to create a list of directories
#dirs=$(find . -maxdepth 1 -type d ! -name "." | sed 's|^\./||')

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <filename> <path>"
    exit 1
fi

# Absolute paths
dirs=$(realpath "$1")
path=$(realpath "$2")

# Check if file exists
if [ ! -f "$dirs" ]; then
    echo "Error: File '$dirs' not found!"
    exit 1
fi

# Check if directory exists
if [ ! -d "$path" ]; then
    echo "Error: Directory '$path' not found!"
    exit 1
fi

# Prompt...

echo "User input is:"
echo "File of Directories of SRL microbes: $dirs"
echo "Path: $path"
microbes=$(wc -l < "$dirs")
echo "Number of SRL microbes: $microbes"

# Read the file line by line
#while IFS= read -r line; do
#    echo "$line"
#    #echo $dirs
#done < "$file"

######################### Automated Unicycler Assembly ######################
cd $path

# initiate conda in the script
source /opt/miniconda3/etc/profile.d/conda.sh
	
# activate the environment of assembly
conda activate autocycler
#exit 0

#----------------------------------------------------------------------#
#----------------------- Quality Filtering ----------------------------#
#----------------------------------------------------------------------#
while IFS= read -r dir_file; do

	dir=$(printf "%s" "$dir_file" | sed 's:[/[:space:]]*$::')
	echo "Processing directory: $dir"

	mkdir -p $dir/reads_qc
	cd $dir

	# fastp quality filtering of short reads
	echo "fastp of $dir"
	fastp --in1 1.Cleandata/$dir*1.fq.gz \
		--in2 1.Cleandata/$dir*2.fq.gz \
		--out1 reads_qc/$dir.QC_1.fq.gz \
		--out2 reads_qc/$dir.QC_2.fq.gz \
		--unpaired1 reads_qc/$dir.QC_1_u.fq.gz \
		--unpaired2 reads_qc/$dir.QC_2_u.fq.gz \
		--thread 16

	# fastplong for filtering the long reads
	echo "fastplong of $dir"
	fastplong \
		-i 1.Cleandata/$dir.filtered_reads.fq.gz \
		-o reads_qc/$dir.QC_long.fq.gz \
		--length_required 1000 \
		--qualified_quality_phred 10 \
		--thread 16

	echo "Finish Processing directory: $dir"
	cd ../

done < $dirs

#----------------------------------------------------------------------#
#----------------------- Unicycler Assembly ---------------------------#
#----------------------------------------------------------------------#
#cd $path

mkdir -p ../logs

rm -f ../logs/assemblies_unicycler_jobs.txt

# build the txt file with all the jobs
while IFS= read -r dir_file; do

	dir=$(printf "%s" "$dir_file" | sed 's:[/[:space:]]*$::')
	# Unicycler for short-reads only:
	echo "unicycler -1 $dir/reads_qc/$dir.QC_1.fq.gz -2 $dir/reads_qc/$dir.QC_2.fq.gz -l $dir/reads_qc/$dir.QC_long.fq.gz -o $dir/unicycler_assembly -t 8" >> ../logs/assemblies_unicycler_jobs.txt

done < $dirs

# run 3 jobs with GNU parallel
jobs=3
max_time=100000 # 100 thousand seconds is 27,7 hours
echo "max_time='$max_time'"

set +e
nice -n 19 parallel --jobs "$jobs" \
	--joblog ../logs/joblog.tsv \
	--results ../logs/logs \
	--delay 2 \
	--timeout "$max_time" < ../logs/assemblies_unicycler_jobs.txt
set -e

# --------------------- Assembly statistics --------------------#
# quast takes about a second per assembly

conda activate quast

while IFS= read -r dir_file; do

	echo $dir_file

	dir=$(echo $dir_file | awk '{print $1}')

	echo "Processing assembly: $dir"
	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
	echo "Directory of assembly: $dir_path"

	cd $dir_path

	assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")

	echo "assembly: $assembly"

	mkdir -p quast

	# quast
	echo "quast analysis of $dir"
	python /opt/miniconda3/envs/quast/bin/quast \
		$assembly \
		-t 12 \
		-o quast 

	cd ../../

	echo "Finish Processing directory: $dir_path"

done < $assemblies


	cd ../

# ------------------------ Annotation ----------------------------#
########################## With Bakta #############################

conda activate bakta

while IFS= read -r dir_file; do

	echo $dir_file

	dir=$(echo $dir_file | awk '{print $1}')

	echo "Processing assembly: $dir"
	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
	echo "Directory of assembly: $dir_path"

	cd $dir_path

	assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")

	echo "assembly: $assembly"

	# bakta
	echo "bakta analysis of $dir"
	bakta $assembly  \
		--db /media/sarlab/DATA/databases/bakta_v6.0/db \
		--threads 12 \
		--output bakta \
		--force

	cd ../../

	echo "Finish Processing directory: $dir_path"

done < $assemblies


# ------------------ de novo wf ----------------#
# this run takes 10 hours

conda activate gtdbtk-2.7.2

gtdbtk de_novo_wf \
	--extension fasta \
	--batchfile batchfile.txt \
	--out_dir gtdbtk_denovo \
	--bacteria \
	--outgroup_taxon p__Chloroflexota \
	--cpus 16 

# outgroup is hard to define https://github.com/Ecogenomics/GTDBTk/issues/390

# ------------------ classify ----------------#
# this run takes 1 hour

gtdbtk classify_wf \
	--batchfile batchfile.txt \
	--out_dir gtdbtk_classify \
	--cpus 16



################### end ##################
	
echo "Finished All Unicycler Assemblies"


time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "Execution time was $time_exec seconds"
