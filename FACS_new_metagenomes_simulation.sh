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

##################################################################################################################################################################
##################################################################################################################################################################
