#!/usr/bin/env bash

set -e

BENCHMARK_DIR=$PWD/BENCHMARK
FACS="/path/to/FACS20.sh"
addTEMP="/tmp/facstempdir"
spurio="/path/to/spurio.py"
db="/path/to/spurio/db/fullsource_filter.fa"
express="/path/to/eXpress"

##################################################################################################################################################################
# Creating directories

mkdir -p $BENCHMARK_DIR/metagenomes/{data FACS}
mkdir -p $BENCHMARK_DIR/metagenomes/data/{metaT bam profiles}
mkdir -p $BENCHMARK_DIR/metagenomes/FACS/{BASHLOG IDS AMP LOG CONTIGS}
out="$BENCHMARK_DIR/metagenomes/FACS/"

##################################################################################################################################################################
# To download metagenomes and metatranscriptomes

for i in $(cat SRR_Acc_List.txt); 
do 
	fastq-dump -I --split-files --outdir $BENCHMARK_DIR/metagenomes/data/ --gzip $i
done

for i in $(cat SRR_Acc_List2.txt); 
do 
	fastq-dump -I --split-files --outdir $BENCHMARK_DIR/metagenomes/data/metaT/ --gzip $i
done

##################################################################################################################################################################
# Processing files with FACS
while read a
do

	echo -e "Doing $a metagenome"
	$FACS -m r --fwd $BENCHMARK_DIR/metagenomes/data/"$a"_1.fastq.gz --rev $BENCHMARK_DIR/metagenomes/data/"$a"_2.fastq.gz 
	\--outfolder $out --outtag $a
	\-t 3 --block 100M --log $a.log
	\--tmp $addTEMP > $out/bashlog.$a.txt
	
	mv $a.log $out/LOG/
	mv $a.ids.tsv.gz $out/IDS/
	mv $a.tsv.gz $out/AMP/
	mv $a.fna.gz $out/CONTIGS/
	
done < SRR_Acc_List.txt

##################################################################################################################################################################
# Then all AMPs were pooled:

cd $out/AMP/
touch wow
zcat * | sed '/Access/d' | cut -f2- | sort -k1,1 | uniq > tmp

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header

cat header tmp | pigz --best > heinz2016.AMP.tsv.gz

rm -rf header tmp

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

zcat heinz2016.AMP.tsv.gz | sed '1,1d' | awk '{print ">"NR"\n"$1}' > tmp.fa
python3 $spurio -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r $db -q tmp.fa -qt spurio_res
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > list
zgrep -v -w -A1 -f list tmp.fa | grep -v ">" > AMP
rm -rf list tmp.fa
zgrep -w -f AMP heinz2016.AMP.tsv.gz > heinz2016.AMP.tsv
rm -rf AMP heinz2016.AMP.tsv.gz
pigz --best heinz2016.AMP.tsv

##################################################################################################################################################################
##################################################################################################################################################################
# Extracting genes

zcat heinz2016.AMP.tsv.gz | cut -f1 > seqs

for i in $(ls $BENCHMARK_DIR/metagenomes/FACS/IDS/)
do
	new=${i/.tsv.gz/}
	zgrep -w -f seqs $BENCHMARK_DIR/metagenomes/FACS/IDS/$i > tmp
	cut -f3 tmp | sed 's/_#_/_/g' | awk -F"_" '{print $1"_"$2":"$4"-"$5"}' > tmp2
	pigz -dc $BENCHMARK_DIR/metagenomes/FACS/CONTIGS/${i/.tsv.gz/.fna.gz} > tmp.fa
	xargs samtools tmp.fa < tmp2 > tmp3
	rm -rf tmp2 tmp
	sed -i "s/>/>$new/g" tmp3
	mv tmp3 $new.tmp tmp.fa*
done
cat *.tmp > AMP.genes.fa
rm -rf seqs *.tmp

# Abundance
bwa index AMP.genes.fa

while read a
do
	bwa mem -p AMP.genes.fa $BENCHMARK_DIR/metagenomes/data/metaT/$a_1.fastq.gz
	\$BENCHMARK_DIR/metagenomes/data/metaT/$a_2.fastq.gz | samtools view -Sb > $BENCHMARK_DIR/metagenomes/data/bam/$a.bam
	samtools sort -n -o $BENCHMARK_DIR/metagenomes/data/bam/$a.bam.sorted $BENCHMARK_DIR/metagenomes/data/bam/$a.bam
	rm -rf $BENCHMARK_DIR/metagenomes/data/bam/$a.bam
done < SRR_Acc_List2.txt

mv AMP.genes.fa $BENCHMARK_DIR/metagenomes/data/FACS/AMP/

ls $BENCHMARK_DIR/metagenomes/data/bam/ | sed -n -e 'H;${x;s/\n/,/g;s/^,//;p;}' > list

## Verify is the list is complete then:
cd $BENCHMARK_DIR/metagenomes/data/bam/
$express --no-update-check --no-bias-correct --calc-covar
	\-o $BENCHMARK_DIR/metagenomes/data/profiles
	\$BENCHMARK_DIR/metagenomes/data/FACS/AMP/AMP.genes.fa
	\$(cat list)
	
rm -rf list



rm -rf AMP.genes.fa*
