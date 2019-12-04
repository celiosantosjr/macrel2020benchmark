#!/usr/bin/env bash

set -e

BENCHMARK_DIR=$PWD/BENCHMARK
FACS="/path/to/FACS20.sh"
addTEMP="/tmp/facstempdir"
ART="/path/to/art_bin_MountRainier/art_illumina"

##################################################################################################################################################################
##################################################################################################################################################################
# creating directories
mkdir -p $BENCHMARK_DIR/refseq
mkdir -p $BENCHMARK_DIR/ab/
mkdir -p $BENCHMARK_DIR/sepgen/
mkdir -p $BENCHMARK_DIR/simulated_data/{40M 60M 80M}
mkdir -p $BENCHMARK_DIR/FACS/{genomes metagenomes abundances}
mkdir -p $BENCHMARK_DIR/FACS/metagenomes/{40M 60M 80M}
add20="$BENCHMARK_DIR/simulated_data/40M"
add30="$BENCHMARK_DIR/simulated_data/60M"
add40="$BENCHMARK_DIR/simulated_data/80M"
out20M="$BENCHMARK_DIR/FACS/metagenomes/40M"
out30M="$BENCHMARK_DIR/FACS/metagenomes/60M"
out40M="$BENCHMARK_DIR/FACS/metagenomes/80M"

# To download genomes
cd $BENCHMARK_DIR/refseq/
wget --header 'Host: progenomes.embl.de' --referer 'http://progenomes.embl.de/representatives.cgi' 'http://progenomes.embl.de/data/repGenomes/representatives.contigs.fasta.gz' --output-document 'representatives.contigs.fasta.gz'

cd $BENCHMARK_DIR/ab/
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466916.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466953.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466965.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621107.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621229.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621247.abund

cat *.abund | sed 's/\.fna//g' | cut -f1 | sort | uniq > list

zgrep -A 1 -w -f list $BENCHMARK_DIR/refseq/representatives.contigs.fasta.gz > $BENCHMARK_DIR/sepgen/repcontigs.fasta
rm -rf $BENCHMARK_DIR/refseq/representatives.contigs.fasta.gz

# Separating reference genomes
cd $BENCHMARK_DIR/sepgen/
for i in $(cat $BENCHMARK_DIR/list)
do
	grep "${i/.fna*}" $BENCHMARK_DIR/sepgen/repcontigs.fasta > tmp
	xargs samtools $BENCHMARK_DIR/sepgen/repcontigs.fasta < tmp > $i
	rm -rf tmp
done
rm -rf $BENCHMARK_DIR/list

# Running FACS in reference genomes:
cat $BENCHMARK_DIR/ab/*.abund | cut -f1 | sort | uniq > .list

for i in $(cat .list)
do 
	echo "Processing genome $i"
	time $FACS -m c --fasta $BENCHMARK_DIR/sepgen/$i --outfolder $BENCHMARK_DIR/FACS/genomes/ --outtag ${i/.fna/} -t 3 --block 100M --log ${i/.fna/.log} --mem 0.75 --tmp $addTEMP --cls 1 --ep 0 > $addFACS/bashlog.${i/.fna/.txt}
done

# Simulating metagenomes
cd $BENCHMARK_DIR/simulated_data/
mkdir ab/
for i in $(ls $BENCHMARK_DIR/ab/); # number of reads
do

	awk '{print $1}' $BENCHMARK_DIR/ab/$i > tmp2
	sort -k1,1 $BENCHMARK_DIR/ab/$i > tmp; mv tmp $BENCHMARK_DIR/ab/$i
	grep -f tmp2 genomes_genomeInfo.txt | sort -k1,1 > sample_genomeInfo.txt # select genomes sizes
	rm -rf tmp2
	join $BENCHMARK_DIR/ab/$i sample_genomeInfo.txt | awk '{print $1"\t"($2*$4)/150"\t"$4}' > cov.list
	rm -rf tmp sample_genomeInfo.txt
	awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s\n",$1,$3,$2/s}' cov.list cov.list > tmp
	awk '{print $1"\t"$3*6000000000/$2}' tmp | sort -k1,1 > in20 # k=40*e6*150 - 20M
	awk '{print $1"\t"$3*9000000000/$2}' tmp | sort -k1,1 > in30 # k=60*e6*150 - 30M
	awk '{print $1"\t"$3*12000000000/$2}' tmp | sort -k1,1 > in40 # k=80*e6*150 - 40M
	mv cov.list ab/$i.1 # coverage per genome
	mv tmp ab/$i.2 # % of reads per genome

	echo "Doing $i"
	while read a b;
	do
		echo -e "\t\t\t$a"
		$ART -na -ss HS25 -i $BENCHMARK_DIR/sepgen/$a -p -l 150 -f $b -m 200 -s 10 -o ${a/.fna} > /dev/null
	done < in20
	cat *1.fq | pigz --best > $add20M/${i/.abund/_1.fastq.gz}
	rm -rf inc20 *1.fq
	cat *2.fq | pigz --best > $add20M/${i/.abund/_2.fastq.gz}
	rm -rf *2.fq

	echo "Doing $i"
	while read a b;
	do
		echo -e "\t\t\t$a"
		$ART -na -ss HS25 -i $BENCHMARK_DIR/sepgen/$a -p -l 150 -f $b -m 200 -s 10 -o ${a/.fna}
	done < in30
	cat *1.fq | pigz --best > $add30M/${i/.abund/_1.fastq.gz}
	rm -rf inc30 *1.fq
	cat *2.fq | pigz --best > $add30M/${i/.abund/_2.fastq.gz}
	rm -rf *2.fq

	echo "Doing $i"
	while read a b;
	do
		echo -e "\t\t\t$a"
		$ART -na -ss HS25 -i $BENCHMARK_DIR/sepgen/$a -p -l 150 -f $b -m 200 -s 10 -o ${a/.fna}
	done < in40
	cat *1.fq | pigz --best > $add40M/${i/.abund/_1.fastq.gz}
	rm -rf *1.fq
	cat *2.fq | pigz --best > $add40M/${i/.abund/_2.fastq.gz}
	rm -rf *2.fq
done

# Renaming simulated reads
zcat $add40/SAMEA2466916_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466916_1.fastq.gz
zcat $add40/SAMEA2466916_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466916_2.fastq.gz
zcat $add40/SAMEA2466953_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466953_1.fastq.gz
zcat $add40/SAMEA2466953_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466953_2.fastq.gz
zcat $add40/SAMEA2466965_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466965_1.fastq.gz
zcat $add40/SAMEA2466965_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2466965_2.fastq.gz
zcat $add40/SAMEA2621107_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621107_1.fastq.gz
zcat $add40/SAMEA2621107_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621107_2.fastq.gz
zcat $add40/SAMEA2621229_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621229_1.fastq.gz
zcat $add40/SAMEA2621229_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621229_2.fastq.gz
zcat $add40/SAMEA2621247_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:40:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621247_1.fastq.gz
zcat $add40/SAMEA2621247_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:40:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add40/SAMEA2621247_2.fastq.gz
zcat $add20/SAMEA2466916_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466916_1.fastq.gz
zcat $add20/SAMEA2466916_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466916_2.fastq.gz
zcat $add20/SAMEA2466953_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466953_1.fastq.gz
zcat $add20/SAMEA2466953_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466953_2.fastq.gz
zcat $add20/SAMEA2466965_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466965_1.fastq.gz
zcat $add20/SAMEA2466965_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2466965_2.fastq.gz
zcat $add20/SAMEA2621107_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621107_1.fastq.gz
zcat $add20/SAMEA2621107_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621107_2.fastq.gz
zcat $add20/SAMEA2621229_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621229_1.fastq.gz
zcat $add20/SAMEA2621229_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621229_2.fastq.gz
zcat $add20/SAMEA2621247_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:20:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621247_1.fastq.gz
zcat $add20/SAMEA2621247_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:20:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add20/SAMEA2621247_2.fastq.gz
zcat $add30/SAMEA2466916_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466916_1.fastq.gz
zcat $add30/SAMEA2466916_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466916:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466916_2.fastq.gz
zcat $add30/SAMEA2466953_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466953_1.fastq.gz
zcat $add30/SAMEA2466953_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466953:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466953_2.fastq.gz
zcat $add30/SAMEA2466965_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466965_1.fastq.gz
zcat $add30/SAMEA2466965_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2466965:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2466965_2.fastq.gz
zcat $add30/SAMEA2621107_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621107_1.fastq.gz
zcat $add30/SAMEA2621107_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621107:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621107_2.fastq.gz
zcat $add30/SAMEA2621229_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621229_1.fastq.gz
zcat $add30/SAMEA2621229_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621229:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621229_2.fastq.gz
zcat $add30/SAMEA2621247_1.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:30:1:1:2:" ++i ":" ++i "#0/1": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621247_1.fastq.gz
zcat $add30/SAMEA2621247_2.fastq.gz | awk '{print (NR%4 == 1) ? "@HS25-SAMEA2621247:30:1:1:2:" ++i ":" ++i "#0/2": $0}' | pigz --best > tmp2; mv tmp2 $add30/SAMEA2621247_2.fastq.gz

# Analyzing reads with FACS
for i in $(ls ab/ | grep ".abund.1")
do 
	time $FACS -m r --fwd $add20/${i/.abund.1/_1.fastq.gz} --rev $add20/${i/.abund.1/_2.fastq.gz} --outfolder $out20M --outtag ${i/.abund.1/.20M} -t 3 --block 100M --log ${i/.abund.1/.20M.log} --mem 0.70 --tmp $addTEMP >  $out20M/bashlog.screening.20M.${i/abund/}

	time $FACS -m r --fwd $add30/${i/.abund.1/_1.fastq.gz} --rev $add30/${i/.abund.1/_2.fastq.gz} --outfolder $out30M --outtag ${i/.abund.1/.30M} -t 3 --block 100M --log ${i/.abund.1/.30M.log} --mem 0.70 --tmp $addTEMP >  $out30M/bashlog.screening.30M.${i/abund/}

	time $FACS -m r --fwd $add40/${i/.abund.1/_1.fastq.gz} --rev $add40/${i/.abund.1/_2.fastq.gz} --outfolder $out40M --outtag ${i/.abund.1/.40M} -t 3 --block 100M --log ${i/.abund.1/.40M.log} --mem 0.70 --tmp $addTEMP >  $out40M/bashlog.screening.40M.${i/abund/}
done

# Converting fq to cram files
for i in $(ls ab/ | grep ".abund.1")
do 
	cut -f1 ab/$i > .tmp

	echo "[W ::: Generating ref file ]"
	for f in $(cat .tmp); do cat $BENCHMARK_DIR/sepgen/$f genomes.temp.fna > temp; mv temp genomes.temp.fna; done

	echo "[W ::: Indexing references ]"
	bwa index genomes.temp.fna
	samtools faidx genomes.temp.fna
 
	echo "[W ::: Mapping metagenome 40M ]"
	bwa mem genomes.temp.fna $add20/40M${i/.abund.1/_1.fastq.gz} $add20/40M${i/.abund.1/_2.fastq.gz} | samtools view -Sb | samtools sort > 40M${i/.abund.1/}.bam
	samtools view -T genomes.temp.fna -C -o 40M${i/.abund.1/}.cram 40M${i/.abund.1/}.bam
	rm -rf 40M${i/.abund.1/}.bam

	echo "[W ::: Mapping metagenome 60M ]"
	bwa mem genomes.temp.fna $add30/60M${i/.abund.1/_1.fastq.gz} $add30/60M${i/.abund.1/_2.fastq.gz} | samtools view -Sb | samtools sort > 60M${i/.abund.1/}.bam
	samtools view -T genomes.temp.fna -C -o 60M${i/.abund.1/}.cram 60M${i/.abund.1/}.bam
	rm -rf 60M${i/.abund.1/}.bam

	echo "[W ::: Mapping metagenome 80M ]"
	bwa mem genomes.temp.fna $add40/80M${i/.abund.1/_1.fastq.gz} $add40/80M${i/.abund.1/_2.fastq.gz} | samtools view -Sb | samtools sort > 80M${i/.abund.1/}.bam
	samtools view -T genomes.temp.fna -C -o 80M${i/.abund.1/}.cram 80M${i/.abund.1/}.bam
	rm -rf 80M${i/.abund.1/}.bam genomes.temp.* .tmp
done
##################################################################################################################################################################
##################################################################################################################################################################
## AMPs from reference genomes and simulated metagenomes were parsed manually
## All AMPs were analyzed manually to curate each cluster accordingly to the same sequence,
## a preliminary list of unique sequences obtained as clusters representatives was generated as follows:
##
## zcat *.tsv.gz | cut -f2 | sed '1,1d' | sort | uniq -c | awk '{print $2"\t"$1}' > t
## echo -e "Sequence\tCluster_size" > he; cat he t > preliminary_clusters; rm -rf he t
##################################################################################################################################################################
## Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>
## Basically command was:
## zcat FACS_output.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
## python3 spurio.py -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q tmp.fa -qt spurio_res
## rm -rf tmp.fa
## awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > spurious_list
##################################################################################################################################################################
##################################################################################################################################################################
