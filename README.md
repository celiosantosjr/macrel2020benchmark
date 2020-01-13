# MACREL Benchmarking (2019/20)

This repository includes code for benchmarking *MACREL*.

This is a companion repository to:

>   FACS: antimicrobial peptide screening in genomes and metagenomes
>   Celio Dias Santos-Junior, Shaojun Pan, Xing-Ming Zhao, Luis Pedro Coelho
>   bioRxiv 2019.12.17.880385; doi:
>   [https://doi.org/10.1101/2019.12.17.880385](https://doi.org/10.1101/2019.12.17.880385)

(The preprint still uses the old name of the tool, _FACS_ and will be updated
soon).

## Contents

It contains the rules to rebuild the benchmarks in the paper.

However, instead just running the code, we strongly recommend you read it.

Some steps depended on inputs obtained from manual curation

- To evaluate benchmarking results over tested AMP prediction models, please refer to the file "benchmark_AMP_models.xlsx".

- To reproduce benchmarking results over hemolytic peptides prediction model implemented in FACS, please follow the code bellow:

```
$ R --vanilla --slave hemolytic_peptides_model_benchmark.R
```

The other results showed in the MACREL benchmarking can be reproduced using the scripts in the following order:

(1) Benchmark.sh

(2) FACS_in_real_metagenomes.sh

(3) Annotation_rules.sh

### Third party softwares

In order to run all the codes, you will need besides MACREL:

- [Spurio](https://bitbucket.org/bateman-group/spurio/src/master/)
- [ArtMountRainier](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
- BlastAll+
- [pigz](https://zlib.net/pigz/)
- R v3.5+
- [samtools](http://samtools.sourceforge.net/)
- [bwa](https://github.com/lh3/bwa)
- [eXpress](https://pachterlab.github.io/eXpress/)
