#!/usr/bin/env bash
set -e
BENCHMARK_DIR=$PWD/BENCHMARK

##################################################################################################################################################################
##################################################################################################################################################################
######################### Analysis of proteins from Complete genome assemblies and representative genomes from proGenomes database.
##################################################################################################################################################################
##################################################################################################################################################################

mkdir -p $BENCHMARK_DIR/analysis
mkdir -p $BENCHMARK_DIR/analysis/genomes
mkdir -p $BENCHMARK_DIR/analysis/metagenomes
mkdir -p $BENCHMARK_DIR/analysis/comparison

##################################################################################################################################################################
# To download genomes

mkdir refseq
mkdir refseq/archaea
mkdir refseq/bacteria
mkdir refseq/viral
mkdir refseq/progenomes

wget --header 'Host: progenomes.embl.de' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:69.0) Gecko/20100101 Firefox/69.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'http://progenomes.embl.de/representatives.cgi' --header 'Cookie: _ga=GA1.2.6012777.1568177574; _pk_id.19.a302=e6b2e563065c012c.1568179528.4.1568630807.1568630773.; _pk_ref.19.a302=%5B%22%22%2C%22%22%2C1568630773%2C%22https%3A%2F%2Fwww.google.com%2F%22%5D; _pk_ses.19.a302=1' --header 'Upgrade-Insecure-Requests: 1' 'http://progenomes.embl.de/data/repGenomes/representatives.contigs.fasta.gz' --output-document 'representatives.contigs.fasta.gz'

##################################################################################################################################################################
## To download genome assemblies we used the ncbi-genome-download device (https://github.com/kblin/ncbi-genome-download)

ncbi-genome-download -F fasta --assembly-level complete bacteria -o refseq/bacteria --parallel 3
ncbi-genome-download -F fasta --assembly-level complete archaea -o refseq/archaea --parallel 3
ncbi-genome-download -F fasta --assembly-level complete viral -o refseq/viral --parallel 3

##################################################################################################################################################################

echo "Searching for AMP in the Representative Contigs"

time FACSv14.sh -m c \
    --fasta $BENCHMARK_DIR/genomes/representatives.contigs.fasta.gz \
    --outfolder $BENCHMARK_DIR/analysis/genomes/ \
    --outtag repcontigs \
    -t 3 \
    --block 100M \
    --log repcontigs.log > $BENCHMARK_DIR/analysis/genomes/bashlog.repcontigs.txt

echo "Doing metagenomes lists"

ls $BENCHMARK_DIR/metagenomes/fastq-files/*_1.fastq.gz > list1
ls $BENCHMARK_DIR/metagenomes/fastq-files/*_2.fastq.gz > list2
paste -d'\t' list1 list2 > readslist.txt; rm -rf list*

while read a b;
do
	echo "Searching for AMP in the metagenomes: $a -- $b"

    base_a=$(basename $a)
	time FACSv14.sh \
        -m r \
        --fwd $a \
        --rev $b \
        --outfolder $BENCHMARK_DIR/analysis/metagenomes/ \
        --outtag ${base_a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${base_a/_1.fastq.gz/}.log \
        --clust 0 > $BENCHMARK_DIR/analysis/metagenomes/bashlog.${base_a/_1.fastq.gz/}.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenomes: $a -- $b"
	time FACSv14.sh -m mr \
        --fwd $a \
        --rev $b \
        --ref $BENCHMARK_DIR/analysis/genomes/representative_contigs.tsv.gz \
        --outfolder $BENCHMARK_DIR/analysis/abundances/ \
        --outtag ${base_a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${base_a/_1.fastq.gz/}.log > $BENCHMARK_DIR/analysis/abundances/bashlog.${base_a/_1.fastq.gz/}.txt

done < readslist.txt

echo "Doing comparisons"

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.tsv.gz | awk '{print ">"$1"\n"$2}' > analysis/comparison/repcontigs.AMP.faa

for i in $(ls analysis/metagenomes/ | grep ".gz"); do zcat analysis/metagenomes/$i | awk '{print ">"$1"\n"$2}' > analysis/comparison/${i/.tsv.gz/.AMP.faa}; done

cd analysis/comparison/
diamond makedb --in repcontigs.AMP.faa --db REPCONTIGS
diamond makedb --in ADAM.faa --db ADAM
diamond makedb --in representativeproteins.faa --db REPprot

for i in from_*.AMP.faa; do diamond blastp -d REPCONTIGS -o ${i/.AMP.faa/.m8} -q $i --sensitive; done
for i in from_*.AMP.faa; do diamond blastp -d ADAM -o ${i/.AMP.faa/.m8.adam} -q $i --sensitive; done
for i in from_*.AMP.faa; do diamond blastp -d REPprot -o ${i/.AMP.faa/REPprot.m8} -q $i; done

for i in $(ls from_*AMP.faa); do grep -c ">" $i >> file1; done
for i in $(ls *.m8); do  cut -f1 $i | sort | uniq | wc -l >> file2; done
for i in $(ls *.m8.adam); do  cut -f1 $i | sort | uniq | wc -l >> file3; done

paste -d'\t' file1 file2 > tmp; rm -rf file1 file2
paste -d'\t' tmp file3 > tmp2; rm -rf tmp file3
echo -e "AMPs\tDetect_on_REP\tDetect_on_ADAM" > header
cat header tmp2 > comparison; rm -rf header tmp2
