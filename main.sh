mkdir -p circularAmp/{data,scripts,ref,bwa,trim,merge}
cd circularAmp 
work_dir=$(pwd)
data=$work_dir/data
trim=$work_dir/trim
merge=$work_dir/merge
script_path=$work_dir/scripts
genome_dir=$work_dir/ref

cd $data
wget https://sequencing.mbl.edu/454/JJR_MDA_Amplicons.tar.gz
tar -zxvf JJR_MDA_Amplicons.tar.gz

### genome size assesment
## count k-mers and produce k-mer frequency histograms
#FastqFiles=$(ls ORI_R[1-2].fastq.gz)
#NumThreads=1
#module load jellyfish2/2.1.1
#for k in 17 21 25 29; do
#  jellyfish count -m $k -s 100M -t $NumThreads -C -o pdom-${k}mers.jf <(zcat $FastqFiles)
#  jellyfish histo pdom-${k}mers.jf > pdom-${k}mers.hist
#done
## estimate k-mer coverage, genome coverage, and genome size.
#cp $script_path/size-coverage-estimate.R .
#./size-coverage-estimate.R &> size.log
## Clean up huge data files.
##rm *.jf


## First trial assembly
module load velvet/1.2.10
# WARNING: This version of Velvet (1.2.10) is compiled for 1 thread, up to 5 categories, and up to Kmer 99.
# To use larger Kmers, please contact HPCC
# To use threading for Velvet, use 1.2.08
lib="Exo" ##"ORI"
mkdir $work_dir/asm_${lib} && cd $work_dir/asm_${lib}
velveth ${lib}.k71.covAuto 71 -shortPaired -fastq.gz -separate $data/${lib}_R1.fastq.gz $data/${lib}_R2.fastq.gz
velvetg ${lib}.k71.covAuto -cov_cutoff auto &> ${lib}.k71.covAuto.log
grep "^>" ${lib}.k71.covAuto/contigs.fa | wc -l ## 47785

cp -R ${lib}.k71.covAuto ${lib}.k71.cov3.len500
velvetg ${lib}.k71.cov3.len500 -cov_cutoff 3 -min_contig_lgth 500 &> ${lib}.k71.cov3.len500.log
grep "^>" ${lib}.k71.cov3.len500/contigs.fa | wc -l ## 137

cp -R ${lib}.k71.covAuto ${lib}.k71.cov3.len1000
velvetg ${lib}.k71.cov3.len1000 -cov_cutoff 3 -min_contig_lgth 1000 &> ${lib}.k71.cov3.len1000.log
grep "^>" ${lib}.k71.cov3.len1000/contigs.fa | wc -l ## 37

## Blast
module load BLAST+/2.2.29
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
DB="nt" ##"human_genomic_transcript" ##"human_genomic"
input=${lib}.k71.cov3.len500/contigs.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn
input=${lib}.k71.cov3.len1000/contigs.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn

## back mapping
## prepare BWA index (for Mapping reads)
mkdir -p $genome_dir/BwaIndex && cd $genome_dir/BwaIndex
cp ../sequence.fasta minicircles.fa
module load bwa/0.7.7.r441
bwa index -a bwtsw minicircles.fa
Bwa_ref="$genome_dir/BwaIndex/minicircles.fa"
cd $work_dir/bwa
bwa mem -t 4 -M $Bwa_ref $data/${lib}_R1.fastq.gz $data/${lib}_R2.fastq.gz > ${lib}_pe_aligned_reads.sam

cat ${lib}_pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -w '0' | wc -l       ## 842025
cat ${lib}_pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' | wc -l      ## 556792
cat ${lib}_pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' > frag_len

## error & adaptor trimming
cd $work_dir/trim
module load Trimmomatic/0.33
java -jar $TRIM/trimmomatic PE -threads 1 -phred33 $data/${lib}_R1.fastq.gz $data/${lib}_R2.fastq.gz ${lib}_R1.pe.fastq ${lib}_R1.se.fastq ${lib}_R2.pe.fastq ${lib}_R2.se.fastq ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:35

## merge overlapping PE reads # https://www.biostars.org/p/225683/
cd $work_dir/merge
module load BBMap/35.34
bbmerge.sh in1=$trim/${lib}_R1.pe.fastq in2=$trim/${lib}_R2.pe.fastq out=${lib}_R.merge.fastq outu1=${lib}_R1.pe.fastq outu2=${lib}_R2.pe.fastq
cat ${lib}_R.merge.fastq $trim/${lib}_R1.se.fastq $trim/${lib}_R2.se.fastq > ${lib}_R.fastq

### genome size assesment
## count k-mers and produce k-mer frequency histograms
#NumThreads=1
#module load jellyfish2/2.1.1
#for k in 17 21 25 29; do
#  jellyfish count -m $k -s 100M -t $NumThreads -C -o pdom-${k}mers.jf ORI_R.fastq
#  jellyfish histo pdom-${k}mers.jf > pdom-${k}mers.hist
#done
## estimate k-mer coverage, genome coverage, and genome size.
#cp $script_path/size-coverage-estimate.R .
#./size-coverage-estimate.R &> size_merge.log
## Clean up huge data files.
#rm *.jf *.hist

###################################################################################################
## second trial assembly
module load velvet/1.2.10
cd $work_dir/asm_${lib}
velveth ${lib}2.k71.covAuto 71 -short -fastq $merge/${lib}_R.fastq
velvetg ${lib}2.k71.covAuto -cov_cutoff auto &> ${lib}2.k71.covAuto.log
grep "^>" ${lib}2.k71.covAuto/contigs.fa | wc -l ## 503576  ## note that the auto changed the cutoff to 0.5 ## 310401

cp -R ${lib}2.k71.covAuto ${lib}2.k71.cov3
velvetg ${lib}2.k71.cov3 -cov_cutoff 3 &> ${lib}2.k71.cov3.log
grep "^>" ${lib}2.k71.cov3/contigs.fa | wc -l ## 545

cp -R ${lib}2.k71.covAuto ${lib}2.k71.cov3.len300
velvetg ${lib}2.k71.cov3.len300 -cov_cutoff 3 -min_contig_lgth 300 &> ${lib}2.k71.cov3.len300.log
grep "^>" ${lib}2.k71.cov3.len300/contigs.fa | wc -l ## 44 ## 30
grep "^>" ${lib}2.k71.cov3.len300/contigs.fa | awk -F"_" '{a+=$4}END{print a}'  ## 42418  ## 26331

cp -R ${lib}2.k71.covAuto ${lib}2.k71.cov3.len500
velvetg ${lib}2.k71.cov3.len500 -cov_cutoff 3 -min_contig_lgth 500 &> ${lib}2.k71.cov3.len500.log
grep "^>" ${lib}2.k71.cov3.len500/contigs.fa | wc -l ## 25

## change the kmer size
module load velvet/1.2.10
cd $work_dir/asm_${lib}
for k in 51 61 81 91 99;do
  velveth ${lib}2.k$k.covAuto $k -short -fastq $merge/${lib}_R.fastq
  velvetg ${lib}2.k$k.covAuto -cov_cutoff auto &> ${lib}2.k$k.covAuto.log
  grep "^>" ${lib}2.k$k.covAuto/contigs.fa | wc -l ## 523414, 512578, 493625, 480668, 463460 ## 320209, 314847, 304673, 296465, 285494

  cp -R ${lib}2.k$k.covAuto ${lib}2.k$k.cov3.len300
  velvetg ${lib}2.k$k.cov3.len300 -cov_cutoff 3 -min_contig_lgth 300 &> ${lib}2.k$k.cov3.len300.log
  grep "^>" ${lib}2.k$k.cov3.len300/contigs.fa | wc -l ## 44, 41, 52, 70, 82  ## 24, 30, 31, 39, 43 
  grep "^>" ${lib}2.k$k.cov3.len300/contigs.fa | awk -F"_" '{a+=$4}END{print a}'  ## 40606, 40405, 45127, 47018, 51744 ## 24295, 25290, 27418, 29547, 28771
done

## try to install your velvet to use longer kmer
conda create --name velvet
source activate velvet
conda install -c bioconda velvet  ## v1.2.10
for k in 105 111 125;do
  velveth ${lib}2.k$k.covAuto $k -short -fastq $merge/${lib}_R.fastq
  velvetg ${lib}2.k$k.covAuto -cov_cutoff auto &> ${lib}2.k$k.covAuto.log
  grep "^>" ${lib}2.k$k.covAuto/contigs.fa | wc -l ## 441945, 405894, 259791 ## 271375, 248317, 156989

  cp -R ${lib}2.k$k.covAuto ${lib}2.k$k.cov3.len300
  velvetg ${lib}2.k$k.cov3.len300 -cov_cutoff 3 -min_contig_lgth 300 &> ${lib}2.k$k.cov3.len300.log
  grep "^>" ${lib}2.k$k.cov3.len300/contigs.fa | wc -l ## 78, 71, 64 ## 32, 33, 26
  grep "^>" ${lib}2.k$k.cov3.len300/contigs.fa | awk -F"_" '{a+=$4}END{print a}'  ## 49606, 45895, 43560 ## 26189, 26996, 25863
done

## try to see the choices of velvet optimiser
conda install -c bioconda perl-velvetoptimiser
mkdir $HOME/miniconda3/envs/velvet/bin/VelvetOpt
wget https://raw.githubusercontent.com/tseemann/VelvetOptimiser/master/VelvetOpt/Assembly.pm
mv Assembly.pm $HOME/miniconda3/envs/velvet/bin/VelvetOpt/.
conda install -c bioconda perl-bioperl 
VelvetOptimiser.pl -s 91 -e 103 -f "-fastq -short ../merge/${lib}_R.fastq" -t 7   ## the best is 103

cp -R auto_data_103 auto_data_k103.cov3.len300  
velvetg auto_data_k103.cov3.len300 -cov_cutoff 3 -min_contig_lgth 300 &> auto_data_k103.cov3.len300.log
grep "^>" auto_data_k103.cov3.len300/contigs.fa | wc -l ## 77  ## 34
grep "^>" auto_data_k103.cov3.len300/contigs.fa | awk -F"_" '{a+=$4}END{print a}'  ## 50457  ## 26519

#VelvetOptimiser.pl -s 103 -e 125 -f "-fastq -short ../merge/${lib}_R.fastq" -t 6   ## the best is 125

## blast
module load BLAST+/2.2.29
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
DB="nt" ##"human_genomic_transcript" ##"human_genomic"
input=ori2.k71.cov3.len300/contigs.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn
input=ori2.k71.cov3.len500/contigs.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn

for input in *.cov3.len300/contigs.fa;do
  blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn
done

## compare assemblies
for f in *.cov3.len300/contigs.fa;do echo $f; grep "^>" $f | awk -F "_" '{print $4}' | sort -nr ;done > asm_len

echo "assembly total_contigs total_asm_nodeLen blast_contigs blast_contig_nodeLen sig_hits match_len" > compare_asm.txt
for asm in *.cov3.len300/contigs.fa;do
 total_contigs=$(grep "^>" $asm | wc -l)
 total_asm_nodeLen=$(grep "^>" $asm | awk -F"_" '{a+=$4}END{print a}')
 blastout=$asm.blastn 
 cat $blastout | awk 'BEGIN {FS=OFS="\t";}{if($4>=100)print;}' | grep -i 'Symbiodinium\|minicircle\|chloroplast' > $blastout.sig
 contigs=$(cat $blastout.sig | awk -F"\t" '{print $1}' | sort | uniq | wc -l)
 contig_nodeLen=$(cat $blastout.sig | awk -F"\t" '{print $1}' | sort | uniq | awk -F"_" '{a+=$4}END{print a}')
 hits=$(cat $blastout.sig | awk -F"\t" '{print $2}' | sort | uniq | wc -l)
 match_len=$(cat $blastout.sig | awk -F"\t" '{a+=$4}END{print a}')
 echo $(dirname $blastout) $total_contigs $total_asm_nodeLen $contigs $contig_nodeLen $hits $match_len
 #cat $blastout | awk 'BEGIN {FS=OFS="\t";}{if($4>=100)print;}' | grep -vi "Escherichia" | wc -l
 #cat $blastout | sort -k2,2 -k4,4nr > $HOME/temp/$(dirname $blastout).blastn
done >> compare_asm.txt

for blastout in *.cov3.len300/contigs.fa.blastn;do
  cat $blastout | sort -k2,2 -k4,4nr > $HOME/temp/$(dirname $blastout).blastn
done 


## back mapping
module load bwa/0.7.7.r441
cd $work_dir/bwa
bwa mem -t 4 -M $Bwa_ref $merge/ORI_R1.pe.fastq $merge/ORI_R2.pe.fastq > ori_unmerged.pe_aligned_reads.sam

cat ori_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -w '0' | wc -l                ## 113629
cat ori_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' | wc -l               ## 10919
cat ori_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' | sed 's/-//' > frag.unmerged_len
cat frag.unmerged_len | awk '{x+=$1}END{print x/NR}'   ## 308.794

## prepare BWA index (for Mapping reads)
mkdir -p $genome_dir/BwaIndex2 && cd $genome_dir/BwaIndex2
cp $work_dir/ori2.k71.cov3/contigs.fa ori2.k71.cov3.fa
module load bwa/0.7.7.r441
bwa index -a bwtsw ori2.k71.cov3.fa
Bwa_ref2="$genome_dir/BwaIndex2/ori2.k71.cov3.fa"
cd $work_dir/bwa
bwa mem -t 4 -M $Bwa_ref2 $merge/ORI_R1.pe.fastq $merge/ORI_R2.pe.fastq > ori2_unmerged.pe_aligned_reads.sam

cat ori2_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -w '0' | wc -l                ## 197030
cat ori2_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' | wc -l               ## 53763
cat ori2_unmerged.pe_aligned_reads.sam | grep -v "^@" | awk '{print $9}' | grep -vw '0' | sed 's/-//' > frag2.unmerged_len
cat frag2.unmerged_len | awk '{x+=$1}END{print x/NR}'  ## 149.521


## third trial assembly
module load velvet/1.2.10
cd $work_dir/asm
velveth ori3.k71.covAuto 71 -short -fastq merge/ORI_R.fastq -shortPaired -fastq -separate merge/ORI_R1.pe.fastq merge/ORI_R2.pe.fastq
velvetg ori3.k71.covAuto -cov_cutoff auto -ins_length 310 &> ori3.k71.covAuto.ins310.log
grep "^>" ori3.k71.covAuto/contigs.fa | wc -l ## 507218  ## note that the auto changed the cutoff to 0.5 

cp -R ori3.k71.covAuto ori3.k71.cov3
velvetg ori3.k71.cov3 -cov_cutoff 3 -ins_length 310 &> ori3.k71.cov3.ins310.log
grep "^>" ori3.k71.cov3/contigs.fa | wc -l ## 759

cp -R ori3.k71.covAuto ori3.k71.cov3.len500
velvetg ori3.k71.cov3.len500 -cov_cutoff 3 -ins_length 310 -min_contig_lgth 500 &> ori3.k71.cov3.ins310.len500.log
grep "^>" ori3.k71.cov3.len500/contigs.fa | wc -l ## 24


## fourth trial assembly
module load velvet/1.2.10
cd $work_dir/asm
velveth ori4.k71.covAuto 71 -short -fastq merge/ORI_R.fastq merge/ORI_R1.pe.fastq merge/ORI_R2.pe.fastq
velvetg ori4.k71.covAuto -cov_cutoff auto &> ori4.k71.covAuto.log
grep "^>" ori4.k71.covAuto/contigs.fa | wc -l ## 507212  ## note that the auto changed the cutoff to 0.5 

cp -R ori4.k71.covAuto ori4.k71.cov3
velvetg ori4.k71.cov3 -cov_cutoff 3 &> ori4.k71.cov3.log
grep "^>" ori4.k71.cov3/contigs.fa | wc -l ## 760

cp -R ori4.k71.covAuto ori4.k71.cov3.len300
velvetg ori4.k71.cov3.len300 -cov_cutoff 3 -min_contig_lgth 300 &> ori4.k71.cov3.len300.log
grep "^>" ori4.k71.cov3.len300/contigs.fa | wc -l ## 45

cp -R ori4.k71.covAuto ori4.k71.cov3.len500
velvetg ori4.k71.cov3.len500 -cov_cutoff 3 -min_contig_lgth 500 &> ori4.k71.cov3.len500.log
grep "^>" ori4.k71.cov3.len500/contigs.fa | wc -l ## 24

## Blast
module load BLAST+/2.2.29
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB
DB="nt" ##"human_genomic_transcript" ##"human_genomic"
input=ori4.k71.cov3.len300/contigs.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn

#####################

