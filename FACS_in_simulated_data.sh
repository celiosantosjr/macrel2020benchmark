#!/usr/bin/env bash
set -e
BENCHMARK_DIR=$PWD/BENCHMARK

##################################################################################################################################################################
##################################################################################################################################################################
######################### Analysis of simulated datasets
##################################################################################################################################################################
##################################################################################################################################################################

mkdir -p $BENCHMARK_DIR/analysis
mkdir -p $BENCHMARK_DIR/analysis/genomes
mkdir -p $BENCHMARK_DIR/analysis/metagenomes
mkdir -p $BENCHMARK_DIR/analysis/comparison
mkdir -p $BENCHMARK_DIR/data
mkdir -p $BENCHMARK_DIR/data/metagenomes
mkdir -p $BENCHMARK_DIR/data/metagenomes/abund/
mkdir -p $BENCHMARK_DIR/data/metagenomes/fastq/
mkdir -p $BENCHMARK_DIR/data/genomes

##################################################################################################################################################################
# Download representative genomes

wget --header 'Host: progenomes.embl.de' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:69.0) Gecko/20100101 Firefox/69.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'http://progenomes.embl.de/representatives.cgi' --header 'Cookie: _ga=GA1.2.6012777.1568177574; _pk_id.19.a302=e6b2e563065c012c.1568179528.4.1568630807.1568630773.; _pk_ref.19.a302=%5B%22%22%2C%22%22%2C1568630773%2C%22https%3A%2F%2Fwww.google.com%2F%22%5D; _pk_ses.19.a302=1' --header 'Upgrade-Insecure-Requests: 1' 'http://progenomes.embl.de/data/repGenomes/representatives.contigs.fasta.gz' --output-document 'representatives.contigs.fasta.gz'

awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' representatives.contigs.fasta.gz > $BENCHMARK_DIR/data/genomes/representatives.contigs.fasta.gz

wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466896.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466916.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466952.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466953.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466965.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2466996.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2467015.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2467039.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621010.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621033.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621107.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621155.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621229.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621247.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2621300.abund
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/SAMEA2622357.abund

cat *.abund | sed 's/\.fna//g' | cut -f1 | sort | uniq > list
mv *.abund $BENCHMARK_DIR/data/metagenomes/abund/

zgrep -A 1 -w -f list $BENCHMARK_DIR/data/genomes/representatives.contigs.fasta.gz > $BENCHMARK_DIR/data/genomes/repcontigs.fasta
rm -rf $BENCHMARK_DIR/data/genomes/representatives.contigs.fasta.gz

##################################################################################################################################################################
## Download simulated metagenomes

wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466896_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466896_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466916_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466916_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466952_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466952_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466953_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466953_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466965_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466965_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466996_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2466996_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2467015_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2467015_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2467039_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2467039_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621010_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621010_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621033_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621033_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621107_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621107_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621155_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621155_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621229_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621229_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621247_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621247_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621300_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2621300_2.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2622357_1.fastq.gz
wget https://zenodo.org/api/files/0f81e436-b108-435c-8339-9fbc9c4daae8/from_SAMEA2622357_2.fastq.gz

mv *.fastq.gz $BENCHMARK_DIR/data/metagenomes/fastq/

ls $BENCHMARK_DIR/data/metagenomes/fastq/*_1.fastq.gz > list1
ls $BENCHMARK_DIR/data/metagenomes/fastq/*_2.fastq.gz > list2
paste -d'\t' list1 list2 > readslist.txt; rm -rf list*

##################################################################################################################################################################
# Getting AMPs

echo "Searching for AMP in the Representative Contigs"

./FACS.sh -m c \
    --fasta $BENCHMARK_DIR/data/genomes/repcontigs.fasta \
    --outfolder $BENCHMARK_DIR/analysis/genomes/ \
    --outtag repcontigs \
    -t 3 \
    --block 100M \
    --log repcontigs.log > $BENCHMARK_DIR/analysis/genomes/bashlog.repcontigs.txt

while read a b;
do
	echo "Searching for AMP in the metagenome ${a/_1.fastq.gz/}"

    base_a=$(basename $a)
	./FACS.sh -m r \
        --fwd $a \
        --rev $b \
        --outfolder $BENCHMARK_DIR/analysis/metagenomes/ \
        --outtag ${a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${a/_1.fastq.gz/}.log \
        --clust 0 > $BENCHMARK_DIR/analysis/metagenomes/bashlog.${a/_1.fastq.gz/}.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenome ${a/_1.fastq.gz/}"
	./FACS.sh -m mr \
        --fwd $a \
        --rev $b \
        --ref $BENCHMARK_DIR/analysis/genomes/repcontigs.tsv.gz \
        --outfolder $BENCHMARK_DIR/analysis/abundances/ \
        --outtag ${a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${a/_1.fastq.gz/}.log > $BENCHMARK_DIR/analysis/abundances/bashlog.${a/_1.fastq.gz/}.txt

done < readslist.txt

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > $BENCHMARK_DIR/analysis/comparison/repcontigs.AMP.faa
python3 spurio.py -s 1 -e $(grep -c ">" $BENCHMARK_DIR/analysis/comparison/repcontigs.AMP.faa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q $BENCHMARK_DIR/analysis/comparison/repcontigs.AMP.faa -qt spurio_res
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > lista
grep -v -w -f lista $BENCHMARK_DIR/analysis/comparison/repcontigs.AMP.tsv.gx > AMP.tsv
mv AMP.tsv $BENCHMARK_DIR/analysis/comparison/repcontigs.nonspu.AMP.tsv

# Then all AMPs were pooled:

cd $BENCHMARK_DIR/analysis/metagenomes/
touch wow

for i in $(ls *gz);
do
	zcat $i | sed '/Access/d' > tmp
	sed -i "s/^/${i/.tsv.gz/}|/g" tmp
	cat tmp wow > tw; rm -rf tmp
	mv tw wow
done

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header
cat header wow | pigz --best > $BENCHMARK_DIR/analysis/metagenomes/AMP.tsv.gz
rm -rf header wow

zcat $BENCHMARK_DIR/analysis/metagenomes/AMP.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
python3 spurio.py -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q tmp.fa -qt spurio_res
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > lista
zgrep -v -w -f lista $BENCHMARK_DIR/analysis/metagenomes/AMP.tsv.gz > AMP.tsv
mv AMP.tsv  $BENCHMARK_DIR/analysis/metagenomes/AMP.nonspu.tsv

##################################################################################################################################################################
## Clustering AMPs

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.tsv.gz | cut -f2 | sed '1,1d' | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' > col1
zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.tsv.gz | cut -f1,2 | sed '1,1d' | awk '{print $2"\t"$1}' | sort -k1,1 |  awk -F'\t' -v OFS=';' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/;;/\t/g' > col2
sort -k1,1 col2 > tmp; mv tmp col2

for i in $($BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.tsv.gz | sed '1,1d' | cut -f2 | sort | uniq);
do
	zgrep "$i" $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.tsv.gz | cut -f2,3,4,5,6 | uniq >> col3
done

sort -k1,1 col3 > tmp; mv tmp col3
join col1 col2 | sed 's/ /\t/g' > tmp
join tmp col3 | sed 's/ /\t/g' > tmp2
rm -rf tmp col*

echo -e "Access\tPeptide\tCluster size\tCluster representatives\tAMP family\tAMP probability\tHemolytic peptide\tHemolytic probability" > header
awk '{print "SIM"NR"\t"$0}' tmp2 > tmp3
cat header tmp3 | pigz --best > $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.clstrs.tsv.gz
rm -rf tmp2 tmp3 header

##################################################################################################################################################################
# Comparison with prokaryotic genomes' AMP 

cd $BENCHMARK_DIR/analysis/comparison/

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.nonspu.AMP.clstrs.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa

blastp -db PAC.db -query tmp.fa -out SIM_clstrs_vsPAC.tsv\
	-evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0\
	-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

rm -rf tmp.fa

##################################################################################################################################################################
# Simulation of new metagenomes

grep ">" $BENCHMARK_DIR/data/genomes/repcontigs.fasta | sed 's/>//g' > headers_rep_genomes.txt
cat headers_rep_genomes.txt | awk -F "." '{print $1"."$2}' | sort | uniq > uniq_genomes
mkdir genomes_separated
zgrep -A 1 -f uniq_genomes repcontigs.fasta > contigs.test.fna
cat contigs.test.fna | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > genomes_genomeInfo.txt

for i in $(cat uniq_genomes);
do
        grep "$i" headers_rep_genomes.txt > tmp
        echo "Extracting sequence genome $i"
        grep -A 1 -w -f tmp contigs.test.fna > genomes_separated_files/$i.fna
        rm -rf tmp
done
rm -rf contigs.test.fna uniq_genomes headers_rep_genomes.txt

for i in $(ls $BENCHMARK_DIR/data/metagenomes/abund/);
do
	awk '{print $1"\t"int($2*100000000)}' $BENCHMARK_DIR/data/metagenomes/abund/$i | sort -k1,1 > $BENCHMARK_DIR/data/metagenomes/abund/${i/.abund/.AbundanceFile.txt}
done # number of reads

for i in $(ls $BENCHMARK_DIR/data/metagenomes/abund/*.AbundanceFile.txt);
do
	awk '{print $1}' $BENCHMARK_DIR/data/metagenomes/abund/$i > tmp # select list
	grep -f tmp genomes_genomeInfo.txt | sort -k1,1 > sample_genomeInfo.txt # select genomes sizes
	join ab/$i sample_genomeInfo.txt | awk '{print $1"\t"$2"\t"(150*$2)/$3}' > tmp
	while read a b c;
	do
		~/art_bin_MountRainier/art_illumina -ss HS25 -i sep_genomes/$a -p -l 150 -f $c -m 200 -s 10 -o ${a/.fna}
	done < tmp
	cat *1.fq | pigz --best > metagenomes_simulated/${i/.AbundanceFile.txt/_1.fastq.gz}
	cat *2.fq | pigz --best > metagenomes_simulated/${i/.AbundanceFile.txt/_2.fastq.gz}
	rm -rf tmp *.fq *.aln sample_genomeInfo.txt 
done
