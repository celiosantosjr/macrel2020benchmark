##################################################################################################################################################################
##################################################################################################################################################################
######################### Annotation and comparison
##################################################################################################################################################################
##################################################################################################################################################################
## To search hits from NCBI database...
## Sequences were uploaded in Blastp server and were searched using Blastp algorithm against nr database and results were downloaded from the webserver as
## a hits table (tsv).
## Searches involving DRAMP database were made as bellow:

# To download database
# 1. DRAMP

mkdir DRAMP

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926912/DRAMP_general_amps.xlsx' --output-document 'DRAMP_general_amps.xlsx'

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926900/DRAMP_patent_amps.xlsx' --output-document 'DRAMP_patent_amps.xlsx'

wget 'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/15926906/DRAMP_clinical_amps.xlsx' --output-document 'DRAMP_clinical_amps.xlsx'

## Manually, the database was cleanned from non-sequences or sequences containing non usual amino acids.
## This new file was saved in DRAMP.tsv

awk '{print ">"$1"\n"$2}' DRAMP.tsv > DRAMP.fa

makeblastdb -in DRAMP.fa -dbtype prot -out DRAMP.db

##################################################################################################################################################################
# Searching with DRAMP:

blastp -db DRAMP.db -query AMPs.fasta -out AMPs.dramp.tsv -evalue 1e-5 -word_size 3 -qcov_hsp_perc 95.0 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

##################################################################################################################################################################
##################################################################################################################################################################
# Parsing of Blast like results:

## 1. from NCBI database webserver
awk '$3 >= 70 && $4 <= 0.00001 && $6 >= 95' blastresults.tsv | sort -k1,1 -k5,5gr -k4,4g -k3,3gr | sort -u -k1,1 --merge > blastresults.parsed.tsv

## 2. using local processing with DRAMP database
awk '$4 >= 70 && $13 <= 0.00001 && $15 >= 95' AMPs.dramp.tsv | sort -k1,1 -k14,14gr -k13,13g -k4,4gr | sort -u -k1,1 --merge > AMPs.dramp.parsed.tsv
awk '$4 >= 70 && $13 <= 0.00001 && $15 >= 95' AMPs.patented.tsv | sort -k1,1 -k14,14gr -k13,13g -k4,4gr | sort -u -k1,1 --merge > AMPs.patented.parsed.tsv
##################################################################################################################################################################
##################################################################################################################################################################
