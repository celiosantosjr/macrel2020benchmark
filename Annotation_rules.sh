##################################################################################################################################################################
##################################################################################################################################################################
######################### Annotation and comparison
##################################################################################################################################################################
##################################################################################################################################################################
## To search hits from NCBI database...
## Sequences were uploaded in Blastp server and were searched using Blastp algorithm against nr database and results were downloaded from the webserver as
## a hits table (tsv).
## Searches involving other databases were made as bellow:
# To download databases
# 1. DRAMP

mkdir DRAMP

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926912/DRAMP_general_amps.xlsx' --output-document 'DRAMP_general_amps.xlsx'

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926900/DRAMP_patent_amps.xlsx' --output-document 'DRAMP_patent_amps.xlsx'

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926906/DRAMP_clinical_amps.xlsx' --output-document 'DRAMP_clinical_amps.xlsx'

## Manually, the database was cleanned from non-sequences or sequences containing non usual amino acids.
## This new file was saved in DRAMP.tsv

awk '{print ">"$1"\n"$2}' DRAMP.tsv > DRAMP.fa

makeblastdb -in DRAMP.fa -dbtype prot -out DRAMP.db

# 2. Patented proteins database

mkdir patented_proteins

cd patented proteins

aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/epop.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/jpop.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/kpop.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/nrnl1.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/nrnl2.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/nrpl1.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/nrpl2.gz
aria2c --file-allocation=none -c -x 10 -s 10 -d ./ ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/patent/uspopnr.gz

pigz -d *.gz

for i in $(ls *); do makeblastdb -in $i -dbtype prot -out ${i/.*/.pt.db}; done

# 3. CDD
## Install as recommended the pipeline rpsbproc
## ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/README

##################################################################################################################################################################
# Searching with DRAMP:

blastp -db DRAMP.db -query AMPs.fasta -out AMPs.dramp.tsv -evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

# Searching with patented proteins database:
for i in $(ls *.pt.db); 
do
	blastp -db DRAMP.db -query AMPs.fasta -out AMPs.${i/.pt.db/_pt}.tsv -evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
done
cat AMPs.*_pt.tsv > AMPs.patented.tsv

##################################################################################################################################################################
# Search with CDD

mkdir splitted/

parallel --pipe --block 1M --recstart ">" "cat > splitted/small-chunk{#}" < input.file >/dev/null 2>/dev/null

for i in $(ls splitted/);
do
	./rpsblast -query splitted/$i -db ./db/Cdd -outfmt 11 -out splitted/$i.aln
done

touch mo

for i in $(ls splitted/ | grep ".aln");
do
	blast_formatter -archive splitted/$i -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out splitted/${i/.aln/.out}
	awk '$11 <= 1e-5 && $3 >= 70 && $13 >= 95' ${i/.aln/.out} > tmp
	cat mo tmp > AMP_cdd.out
	mv AMP_cdd.out mo
	rm -rf tmp
done

cat splitted/*out | awk '$11 <= 0.00001' > AMP_cdd.out

cut -f2 AMP_cdd.out | sed 's/|/\t/g' | awk '{print $3}' > col

for i in $(cat col);
do
	va=`grep "$i" data/bitscore_specific.txt | awk '{print $3}'`
	fun=`grep "$i" data/cddid.tbl`
	echo $va >> val
	echo $fun >> fun
	unset va
	unset fun
done

sed -i 's/ /_/g' fun

paste -d '\t' AMP_cdd.out val fun | sed 's/ /\t/g' | awk '$12 >= $14' > AMP_cdd.tsv

rm -rf val fun col splitted/

##################################################################################################################################################################
# Parsing of Blast like results:

## 1. from NCBI database webserver
awk '$3 >= 70 && $4 <= 0.00001 && $6 >= 95' blastresults.tsv | sort -k1,1 -k5,5gr -k4,4g -k3,3gr | sort -u -k1,1 --merge > blastresults.parsed.tsv

## 2. using local processing with DRAMP and patented proteins database
awk '$4 >= 70 && $13 <= 0.00001 && $15 >= 95' AMPs.dramp.tsv | sort -k1,1 -k14,14gr -k13,13g -k4,4gr | sort -u -k1,1 --merge > AMPs.dramp.parsed.tsv
awk '$4 >= 70 && $13 <= 0.00001 && $15 >= 95' AMPs.patented.tsv | sort -k1,1 -k14,14gr -k13,13g -k4,4gr | sort -u -k1,1 --merge > AMPs.patented.parsed.tsv

## 3. CDD results
LANG=C; LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr AMP_cdd.tsv | sort -u -k1,1 --merge | awk '$13 >= 95 && $3 >= 70' > tmp; mv tmp AMP_cdd.tsv

##################################################################################################################################################################
##################################################################################################################################################################
