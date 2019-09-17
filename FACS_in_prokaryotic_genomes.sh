##################################################################################################################################################################
##################################################################################################################################################################
######################### Analysis of proteins from Complete genome assemblies and representative genomes from proGenomes database.
##################################################################################################################################################################
##################################################################################################################################################################
# To download genomes

mkdir refseq
mkdir refseq/archaea
mkdir refseq/bacteria
mkdir refseq/viral
mkdir refseq/progenomes

wget --header 'Host: progenomes.embl.de' --referer 'http://progenomes.embl.de/representatives.cgi' 'http://progenomes.embl.de/data/repGenomes/representatives.contigs.fasta.gz' --output-document 'representatives.contigs.fasta.gz'

##################################################################################################################################################################
## To download genome assemblies we used the ncbi-genome-download device (https://github.com/kblin/ncbi-genome-download)

ncbi-genome-download -F fasta --assembly-level complete bacteria -o refseq/bacteria --parallel 3
ncbi-genome-download -F fasta --assembly-level complete archaea -o refseq/archaea --parallel 3
ncbi-genome-download -F fasta --assembly-level complete viral -o refseq/viral --parallel 3

##################################################################################################################################################################
# Once all genomes downloaded
# Create a list with 3 columns: file, species, kingdom(archaea, bacteria, viral)

mv representatives.contigs.fasta.gz refseq/progenomes/
zgrep ">" refseq/progenomes/representatives.contigs.fasta.gz > col
sed 's/\..*//g' col > col2
for i in $(cat col2); do echo "bacteria" >> col3; done
paste -d'\t' col col2 col3 > progenomes.tab; rm -rf col col2 col3

ncbi-genome-download -F fasta --assembly-level complete viral --parallel 3 --dry-run | sed 's/ /_/g' > genomes_viral_list.txt
ncbi-genome-download -F fasta --assembly-level complete archaea --parallel 3 --dry-run | sed 's/ /_/g' > genomes_bacteria_list.txt
ncbi-genome-download -F fasta --assembly-level complete bacteria --parallel 3 --dry-run | sed 's/ /_/g' > genomes_archaea_list.txt

cp genomes_archaea_list.txt tmp; cut -f1 genomes_archaea_list.txt | sed '1,1d' | sort > t; mv t genomes_archaea_list.txt; for i in $(ls refseq/archaea/); do echo "$i" > tp; inc=`awk -F"_" '{print $1"_"$2}' tp`; name=`grep "$inc" genomes_archaea_list.txt | sed 's/ /_/g' | awk '{$1=""}1' | sed 's/\t//g'`; echo -e "$i\t$name\tarchaea" >> final_list; mv tmp genomes_archaea_list.txt

cp genomes_bacteria_list.txt tmp; cut -f1 genomes_bacteria_list.txt | sed '1,1d' | sort > t; mv t genomes_bacteria_list.txt; for i in $(ls refseq/bacteria/); do echo "$i" > tp; inc=`awk -F"_" '{print $1"_"$2}' tp`; name=`grep "$inc" genomes_bacteria_list.txt | sed 's/ /_/g' | awk '{$1=""}1' | sed 's/\t//g'`; echo -e "$i\t$name\tbacteria" >> final_list; mv tmp genomes_bacteria_list.txt

cp genomes_viral_list.txt tmp; cut -f1 genomes_viral_list.txt | sed '1,1d' | sort > t; mv t genomes_viral_list.txt; for i in $(ls refseq/viral/); do echo "$i" > tp; inc=`awk -F"_" '{print $1"_"$2}' tp`; name=`grep "$inc" genomes_viral_list.txt | sed 's/ /_/g' | awk '{$1=""}1' | sed 's/\t//g'`; echo -e "$i\t$name\tviral" >> final_list; mv tmp genomes_viral_list.txt

# Then run:

mkdir FACS
mkdir FACS/bacteria
mkdir FACS/archaea
mkdir FACS/viral
mkdir FACS/repcontigs

while read a b c
do

	echo -e "Doing $a genome from $b -- SK: $c"
	./FACS.sh -m c --fasta refseq/$a --outfolder FACS/$c --outtag ${a/.fasta.gz/} -t 3 --block 100M --log $b.log > FACS/$c/bashlog.$b.txt

done < final.list

./FACS.sh -m c --fasta refseq/progenomes/representatives.contigs.fasta.gz --outfolder FACS/progenomes/ --outtag proGenomes -t 3 --block 100M --log proGenomes.log > FACS/$c/bashlog.proGenomes.txt


##################################################################################################################################################################
# All results are then organized into:

mkdir AMP
mkdir BASHLOG
mkdir IDS
mkdir LOGs

cd FACS/archaea; mv *.logs LOGSs/; mv bashlog.* BASHLOG/; mv *ids.tsv.gz IDS/; mv *.tsv.gz AMP/; cd ../
cd FACS/bacteria; mv *.logs LOGSs/; mv bashlog.* BASHLOG/; mv *ids.tsv.gz IDS/; mv *.tsv.gz AMP/; cd ../
cd FACS/viral; mv *.logs LOGSs/; mv bashlog.* BASHLOG/; mv *ids.tsv.gz IDS/; mv *.tsv.gz AMP/; cd ../

##################################################################################################################################################################
# Then all AMPs were pooled:

mkdir all

cp FACS/bacteria/AMP/* all/
cp FACS/archaea/AMP/* all/
cp FACS/viral/AMP/* all/
cp FACS/progenomes/proGenomes.tsv.gz all/

cd all/

for i in $(ls );
do
	zcat $i | sed '1,1d' > tmp
	f=${i/.tsv.gz/}
	sed "s/smORF_/$f|smORF/g" tmp >> AMPs.tsv
	rm -rf tm
done

while read a b c d e f;
do
	echo "$a" | sed 's/|.*//g' > tmp
	fi=`cat tmp`
	rm -rf tmp
	val=`grep "$fi" ../final.list | cut -f3`
	echo -e "$a\t$b\t$c\t$d\t$e\t$f\tbacteria\t$val" >> tpe
done < AMPs.tsv; mv tpe AMPs.tsv

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability\tKingdom\tSpecies" > header

cat header AMPs.tsv > tmp; mv tmp AMPs.tsv

##################################################################################################################################################################
## AMPs from proGenomes were parsed manually using the numbers coding to species that were searched in the proGenomes database.
## These AMPs were then replicated in AMPs.TSV file accordingly to their original contigs present in the ids.tsv.gz files.
## All AMPs were analyzed manually to curate each cluster accordingly to the same sequence, a preliminary list of unique sequences
## obtained as clusters representatives was obtained as follows:

zcat AMPs.tsv.gz | cut -f2 | sed '1,1d' | sort | uniq -c | awk '{print $2"\t"$1}' > t
echo -e "Sequence\tCluster_size" > he; cat he t > preliminary_clusters; rm -rf he t

## Once all AMPs from assemblies and proGenomes database were in AMPs file, then:

pigz --best AMPs.tsv

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

zcat AMP.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
python3 spurio.py -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q tmp.fa -qt spurio_res
rm -rf tmp.fa
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > list
pigz -dc AMP.tsv.gz > tmp
grep -v -w -f list tmp > AMP.tsv
cut -f4 AMP.tsv | sed '1,1d' > nonspurious.scores
grep -w -f list tmp | cut -f4 > spurious.scores
echo "spurious" > t; cat t spurious.scores > tt; rm -rf t; mv tt spurious.scores
rm -rf AMP.tsv.gz tmp list
paste -d'\t' nonspurious.scores spurious.scores > scores_comparison; rm -rf *.scores

# To generate figure 2 it was used a R script with scores_comparison table:
R --vanilla --slave boxplot_chart_maker.R

##################################################################################################################################################################
##################################################################################################################################################################


