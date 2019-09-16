##################################################################################################################################################################
##################################################################################################################################################################
######################### Analysis of proteins from Complete genome assemblies and representative genomes from proGenomes database.
##################################################################################################################################################################
##################################################################################################################################################################
# To download metagenomes

mkdir metagenomes/
mkdir metagenomes/data
mkdir metagenoma/FACS
mkdir metagenoma/FACS/BASHLOG/
mkdir metagenoma/FACS/IDS/
mkdir metagenoma/FACS/AMP
mkdir metagenoma/FACS/LOG
cd metagenomes/data/

for i in $(cat qinetal2010.list); 
do 
	wget http://www.bork.embl.de/~arumugam/Qin_et_al_2010/$i
done

cd ../

##################################################################################################################################################################

while read a
do

	echo -e "Doing $a metagenome"
	./FACS.sh -m c --fasta data/$a --outfolder FACS/ --outtag ${a/.seq.fa.gz/}\
	-t 3 --block 100M --log ${a/.seq.fa.gz/}.log > FACS/BASHLOG/bashlog.${a/.seq.fa.gz/}.txt
	
	mv *.log FACS/LOG/
	mv *.ids.tsv.gz FACS/IDS/
	mv *.tsv.gz FACS/AMP/
	
done < qinetal2010.list

##################################################################################################################################################################
# Then all AMPs were pooled:

cd AMP/
touch wow

for i in $(ls *gz);
do
	zcat $i | sed '/Access/d' > tmp
	sed -i "s/^/${i/.tsv.gz/}|/g" tmp
	cat tmp wow > tw; rm -rf tmp
	mv tw wow
done

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header

cat header wow | pigz --best > qin2010_humangut.AMP.tsv.gz

rm -rf header wow

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

zcat qin2010_humangut.AMP.tsv.gz | awk '{print ">"$1"\n"$2}' | sed '1,2d' > tmp.fa
python3 spurio.py -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r /path/to/spurio/db/fullsource_filter.fa -q tmp.fa -qt spurio_res
rm -rf tmp.fa
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > list
zgrep -v -w -f list qin2010_humangut.AMP.tsv.gz > AMP.tsv
mv AMP.tsv qin2010_humangut.AMP.tsv
pigz --best qin2010_humangut.AMP.tsv

##################################################################################################################################################################
## Clustering AMPs

zcat qin2010_humangut.AMP.tsv.gz | cut -f2 | sed '1,1d' | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' > col1
zcat qin2010_humangut.AMP.tsv.gz | cut -f1,2 | sed '1,1d' | awk '{print $2"\t"$1}' | sort -k1,1 |  awk -F'\t' -v OFS=';' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/;;/\t/g' > col2
sort -k1,1 col2 > tmp; mv tmp col2

for i in $(zcat qin2010_humangut.AMP.tsv.gz | sed '1,1d' | cut -f2 | sort | uniq);
do
	zgrep "$i" qin2010_humangut.AMP.tsv.gz | cut -f2,3,4,5,6 | uniq >> col3
done

sort -k1,1 col3 > tmp; mv tmp col3
join col1 col2 | sed 's/ /\t/g' > tmp
join tmp col3 | sed 's/ /\t/g' > tmp2
rm -rf tmp col*

echo -e "Access\tPeptide\tCluster size\tCluster representatives\tAMP family\tAMP probability\tHemolytic peptide\tHemolytic probability" > header
awk '{print "QAC"NR"\t"$0}' tmp2 > tmp3
cat header tmp3 | pigz --best > qinetal2010_metagenomes.clstrs.tsv.gz
rm -rf tmp2 tmp3 header

##################################################################################################################################################################
# Comparison with prokaryotic genomes' AMPs
zcat AMPs.tsv.gz | sed '1,1d' | awk '{print ">"$1"\n"$2}' > tmp.fa
zcat qinetal2010_metagenomes.clstrs.tsv.gz | sed '1,1d' | awk '{print ">"$1"\n"$2}' > qin.fa

makeblastdb -in tmp.fa -dbtype prot -out PAC.db
rm -rf tmp.fa

blastp -db PAC.db -query qin.fa -out AMPs.QAC.tsv\
-evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0\
-outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

awk '$4 >= 70 && $13 <= 0.00001 && $15 >= 95' AMPs.QAC.tsv | sort -k1,1 -k14,14gr -k13,13g -k4,4gr | sort -u -k1,1 --merge > AMPs.QAC.parsed.tsv

##################################################################################################################################################################
##################################################################################################################################################################
