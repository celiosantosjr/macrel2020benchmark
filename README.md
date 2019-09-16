# FACS Benchmarking (2019)

This repository includes code for benchmarking *FACS*.

This is a companion repository to:

> _FACS: FACS - Fast antimicrobial peptides screening in high throughput data_
> by Célio Dias Santos Júnior and Luis Pedro Coelho _in_ Bionformatics XX:YY (2019)
> [http://doi.org/XXX](http://doi.org/XXX)

It contains the rules to rebuild the benchmark proposed in the paper.

However, instead just running the code, we strongly recommend to read it.

Some steps depended on inputs obtained from manual curation, impossible to reproduce easily here.

## Order of scripts

The order to reproduce the entire paper benchmark is:

(1) FACS_in_prokaryotic_genomes.sh

(2) FACS_in_real_metagenomes.sh

(3) FACS_in_simulated_data.sh

(4) FACS_new_metagenomes_simulation.sh

(5) Annotation_rules.sh


## Third party softwares

In order to run all the codes, you will need:

- Spurio (https://bitbucket.org/bateman-group/spurio/src/master/)
- ArtMountRainier (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
- BlastAll+
- CDD (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd)
- ncbi-genome-download (https://github.com/kblin/ncbi-genome-download)
- pigz (https://zlib.net/pigz/)
- R v. 3.5+
- R package ggpubr (https://www.rdocumentation.org/packages/ggpubr/versions/0.2.3)
- R package dplyr (https://www.rdocumentation.org/packages/dplyr/versions/0.7.8)
