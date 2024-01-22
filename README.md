# Maine-eDNA Metabarcoding Tutorial

This self-guided workflow tutorial will cover analysis of Illumina MiSeq data from initial reads to amplicon sequence variants (ASVs), taxonomy, and visualization. The workflow uses data generated by Maine-eDNA researchers for 18S using the [Comeau et al., 2011](https://doi.org/10.1371/journal.pone.0027492) primers E572F/E1009R for the V4 region. A containerized environment was created in CyVerse for software setup, which uses FastQC, Cutadapt, and R packages for analysis and visualization {Dada2, Phyloseq, GGplot2}.

> For more information on using CyVerse, see the [help documentation here](https://github.com/umaine-research/metabarcoding-cyverse/wiki).

### Suggested Publications

- [Bohmann et al. 2014](https://www.sciencedirect.com/science/article/pii/S016953471400086X?via%3Dihub) - Environmental DNA for wildlife biology and biodiversity monitoring. TREE, 29.
- [Schnell et al. 2015](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12402) - Tag jumps illuminated – reducing sequence‐to‐sample misidentifications in metabarcoding studies. Molecular Ecology Resources, 15, 1289-1303.
- [Callahan et al. 2016](https://www.nature.com/articles/nmeth.3869) - DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13, 581–583.
- [Ficetola et al. 2016](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12508) - How to limit false positives in environmental DNA and metabarcoding? Molecular Ecology, 16, 604-607.
- [Callahan et al. 2017](https://doi.org/10.1038/ismej.2017.119) - Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. ISME J, 11, 2639–2643.
- [Deiner et al. 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14350) - Environmental DNA metabarcoding: Transforming how we survey animal and plant communities. Molecular Ecology, 26, 5872-5895.
- [Alberdi et al. 2017](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12849) - Scrutinizing key steps for reliable metabarcoding of environmental samples. Methods in Ecology and Evolution, 9, 134-147.
- [Dickie et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12907) - Towards robust and repeatable sampling methods in eDNA-based studies. Molecular Ecology Resources, 18, 940-952.
- [Lamb et al. 2019](https://onlinelibrary.wiley.com/doi/10.1111/mec.14920) - How quantitative is metabarcoding: A meta‐analytical approach. Molecular Ecology, 28, 420-430. -Jamy et al. 2019. Long-read metabarcoding of the eukaryotic rDNA operon to phylogenetically and taxonomically resolve environmental diversity. Molecular Ecology Resources, 20, 429-443.
- [Zinger et al. 2019](https://onlinelibrary.wiley.com/doi/10.1111/mec.15060) - DNA metabarcoding—Need for robust experimental designs to draw sound ecological conclusions. Molecular Ecology, 28, 1857-1862
- [Antich et al. 2021](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04115-6) - To denoise or to cluster, that is not the question: optimizing pipelines for COI metabarcoding and metaphylogeography. BMC Bioinformatics, 22,177.

