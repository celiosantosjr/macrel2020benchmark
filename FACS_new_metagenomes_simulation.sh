#!/usr/bin/env bash
set -e
BENCHMARK_DIR=$PWD/BENCHMARK

##################################################################################################################################################################
##################################################################################################################################################################
######################### Analysis of new simulated datasets
##################################################################################################################################################################
##################################################################################################################################################################
# Simulation of new metagenomes

grep ">" $BENCHMARK_DIR/data/genomes/repcontigs.fasta | sed 's/>//g' > headers_rep_genomes.txt
cat headers_rep_genomes.txt | awk -F "." '{print $1"."$2}' | sort | uniq > uniq_genomes
mkdir genomes_separated
zgrep -A 1 -f uniq_genomes repcontigs.fasta > contigs.test.fna

for i in $(cat uniq_genomes);
do
        grep "$i" headers_rep_genomes.txt > tmp
        echo "Extracting sequence genome $i"
        grep -A 1 -w -f tmp contigs.test.fna > genomes_separated/$i.fna
        rm -rf tmp
done
rm -rf contigs.test.fna uniq_genomes headers_rep_genomes.txt

for i in $(ls genomes_separated/);
do
	val=`grep -v ">" genomes_separated/$i | wc | awk '{print $3-$1}'`
	echo -e "$i\t$val" >> genomes_genomeInfo.txt
	unset val
done

for i in $(ls $BENCHMARK_DIR/data/metagenomes/abund/);
do
	awk '{print $1"\t"int($2*100000000)}' $BENCHMARK_DIR/data/metagenomes/abund/$i | sort -k1,1 > ${i/.abund/.AbundanceFile.txt}
done # number of reads

mkdir $BENCHMARK_DIR/data/newmetagenomes

# To simulate reads it was used ArtBin_MountRainier software (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

for i in $(ls *.AbundanceFile.txt);
do
	awk '{print $1}' $i > tmp # select list
	grep -f tmp genomes_genomeInfo.txt | sort -k1,1 > sample_genomeInfo.txt # select genomes sizes
	join $i sample_genomeInfo.txt | awk '{print $1"\t"$2"\t"(150*$2)/$3}' > tmp
	while read a b c;
	do
		~/art_bin_MountRainier/art_illumina -ss HS25 -i separated_genomes/$a -p -l 150 -f $c -m 200 -s 10 -o ${a/.fna}
	done < tmp
	cat *1.fq | pigz --best > $BENCHMARK_DIR/data/newmetagenomes/${i/.AbundanceFile.txt/_1.fastq.gz}
	cat *2.fq | pigz --best > $BENCHMARK_DIR/data/newmetagenomes/${i/.AbundanceFile.txt/_2.fastq.gz}
	rm -rf tmp *.fq *.aln sample_genomeInfo.txt 
done

ls $BENCHMARK_DIR/data/newmetagenomes | grep "_1.fastq.gz" > list1
ls $BENCHMARK_DIR/data/newmetagenomes | grep "_2.fastq.gz" > list2
paste -d'\t' list1 list2 >  reads.list

##################################################################################################################################################################
# Getting AMPs

mkdir $BENCHMARK_DIR/analysis/metagenomes/newsim/
mkdir $BENCHMARK_DIR/analysis/abundances/newsim/

while read a b;
do
	echo "Searching for AMP in the metagenome ${a/_1.fastq.gz/}"

	./FACS.sh -m r \
        --fwd $BENCHMARK_DIR/data/newmetagenomes/$a \
        --rev $BENCHMARK_DIR/data/newmetagenomes/$b \
        --outfolder $BENCHMARK_DIR/analysis/metagenomes/newsim/ \
        --outtag ${a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${a/_1.fastq.gz/}.log \
        --clust 0 > $BENCHMARK_DIR/analysis/metagenomes/newsim/bashlog.${a/_1.fastq.gz/}.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenome ${a/_1.fastq.gz/}"
	./FACS.sh -m mr \
        --fwd $a \
        --rev $b \
        --ref $BENCHMARK_DIR/analysis/genomes/repcontigs.tsv.gz \
        --outfolder $BENCHMARK_DIR/analysis/abundances/newsim/ \
        --outtag ${a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${a/_1.fastq.gz/}.log > $BENCHMARK_DIR/analysis/abundances/newsim/bashlog.${a/_1.fastq.gz/}.txt

done < readslist.txt

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

cd $BENCHMARK_DIR/analysis/metagenomes/newsim/
touch wow

for i in $(ls *gz);
do
	zcat $i | sed '/Access/d' > tmp
	sed -i "s/^/${i/.tsv.gz/}|/g" tmp
	cat tmp wow > tw; rm -rf tmp
	mv tw wow
done

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header
cat header wow | pigz --best > $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.tsv.gz
rm -rf header wow

zcat $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
python3 spurio.py -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q tmp.fa -qt spurio_res
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > lista
zgrep -v -w -f lista $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.tsv.gz > AMP.tsv
mv AMP.tsv  $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.nonspu.tsv

##################################################################################################################################################################
# Comparison with prokaryotic genomes' AMP 

cd $BENCHMARK_DIR/analysis/comparison

zcat $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.nonspu.tsv | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp2.fa

blastp -db PAC.db -query tmp2.fa -out NEWSIM_clstrs_vsPAC.tsv\
	-evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0\
	-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

rm -rf tmp2.fa

##################################################################################################################################################################
# Comparison with representative contigs AMPs

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.clstrs.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
zcat $BENCHMARK_DIR/analysis/metagenomes/newsim/AMP.nonspu.tsv | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp2.fa

makeblastdb -in tmp.fa -dbtype prot -out repcontigs.db
rm -rf tmp

blastp -db repcontigs.db -query tmp2.fa -out NEWSIM_clstrs_vsRepContigs.tsv\
	-evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0\
	-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

rm -rf tmp2.fa

##################################################################################################################################################################
##################################################################################################################################################################
