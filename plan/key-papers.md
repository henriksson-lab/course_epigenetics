# Key papers: seed list per unit

Candidates for each unit's optional "Key papers" list (designprinciples §9).
Written from memory on 2026-11-17: **every entry must be checked (authors,
year, journal, DOI) before it goes into a lecture or exercise.** Tick `[x]`
when verified and add the DOI. Add, swap and drop freely; 3–6 per unit.

## Foundations

`chromatin-intro`
- [x] Luger et al. 1997, *Nature*: crystal structure of the nucleosome core particle. doi:10.1038/38444
- [x] Strahl & Allis 2000, *Nature*: the "histone code" hypothesis. doi:10.1038/47412
- [x] Bird 2007, *Nature*: "Perceptions of epigenetics", the definition debate. doi:10.1038/nature05913
- [x] Allis & Jenuwein 2016, *Nat Rev Genet*: the molecular hallmarks of epigenetic control (review). doi:10.1038/nrg.2016.59
- [x] Kornberg 1974, *Science*: chromatin structure, a repeating unit of histones and DNA. doi:10.1126/science.184.4139.868
- [x] Barski et al. 2007, *Cell*: high-resolution profiling of histone methylations. doi:10.1016/j.cell.2007.05.009

`gene-regulation`
- [ ] Banerji, Rusconi & Schaffner 1981, *Cell*: the first enhancer (SV40).
- [ ] Lambert et al. 2018, *Cell*: "The human transcription factors" census.
- [ ] Zaret & Carroll 2011, *Genes Dev*: pioneer transcription factors.

`genomic-intervals`
- [x] Quinlan & Hall 2010, *Bioinformatics*: BEDTools. doi:10.1093/bioinformatics/btq033
- [x] Lawrence et al. 2013, *PLoS Comput Biol*: GenomicRanges / Bioconductor infrastructure. doi:10.1371/journal.pcbi.1003118
- [x] Li et al. 2009, *Bioinformatics*: the SAM format and SAMtools. doi:10.1093/bioinformatics/btp352
- [x] Kent et al. 2010, *Bioinformatics*: BigWig and BigBed. doi:10.1093/bioinformatics/btq351
- [x] Frankish et al. 2019, *Nucleic Acids Res*: GENCODE. doi:10.1093/nar/gky955
- [x] ENCODE Project Consortium 2020, *Nature*: expanded encyclopaedias. doi:10.1038/s41586-020-2493-4

`enrichment-alignment`
- [ ] Langmead & Salzberg 2012, *Nat Methods*: Bowtie 2.
- [ ] Amemiya, Kundaje & Boyle 2019, *Sci Rep*: the ENCODE blacklist.
- [ ] Ramírez et al. 2016, *Nucleic Acids Res*: deepTools2.
- [ ] Robinson et al. 2011, *Nat Biotechnol*: IGV.

`epigenomics-design`
- [ ] ENCODE Project Consortium 2012, *Nature*: integrated encyclopedia of DNA elements.
- [ ] Landt et al. 2012, *Genome Res*: ChIP-seq guidelines and practices of ENCODE/modENCODE.
- [ ] Roadmap Epigenomics Consortium 2015, *Nature*: 111 reference human epigenomes.

## ATAC-seq

`atacseq-intro`
- [ ] Buenrostro et al. 2013, *Nat Methods*: ATAC-seq.
- [ ] Corces et al. 2017, *Nat Methods*: Omni-ATAC.
- [ ] Thurman et al. 2012, *Nature*: the accessible chromatin landscape (DNase-seq).

`atacseq-analysis`
- [ ] Yan et al. 2020, *Genome Biol*: ATAC-seq data analysis review.
- [ ] Zhang et al. 2008, *Genome Biol*: MACS.

`differential-accessibility`
- [ ] Love, Huber & Anders 2014, *Genome Biol*: DESeq2 (already known; the model reused).
- [ ] Lun & Smyth 2016, *Nucleic Acids Res*: csaw, normalisation pitfalls in differential binding.
- [ ] Reske, Wilson & Edwards 2020, *Epigenetics Chromatin*: normalisation choices in differential ATAC-seq.

`motifs-footprinting`
- [ ] Bailey et al. 2009, *Nucleic Acids Res*: MEME Suite.
- [ ] Heinz et al. 2010, *Mol Cell*: HOMER (and lineage-determining TFs).
- [ ] Schep et al. 2017, *Nat Methods*: chromVAR.
- [ ] Bentsen et al. 2020, *Nat Commun*: TOBIAS footprinting.

## ChIP-seq and successors

`chipseq-intro`
- [ ] Johnson et al. 2007, *Science*: genome-wide mapping of in vivo protein–DNA interactions (ChIP-seq).
- [ ] Barski et al. 2007, *Cell*: high-resolution profiling of histone methylations.
- [ ] Park 2009, *Nat Rev Genet*: ChIP-seq advantages and challenges (review).

`chipseq-peaks`
- [ ] Zhang et al. 2008, *Genome Biol*: MACS (Poisson local lambda).
- [ ] Li et al. 2011, *Ann Appl Stat*: IDR for replicate consistency.

`cutrun-cuttag`
- [ ] Skene & Henikoff 2017, *eLife*: CUT&RUN.
- [ ] Kaya-Okur et al. 2019, *Nat Commun*: CUT&Tag.
- [ ] Meers, Tenenbaum & Henikoff 2019, *Epigenetics Chromatin*: SEACR peak calling.

`chipseq-downstream`
- [ ] Yu, Wang & He 2015, *Bioinformatics*: ChIPseeker.
- [ ] Ross-Innes et al. 2012, *Nature*: differential ER binding and outcome in breast cancer (DiffBind).
- [ ] Orlando et al. 2014, *Cell Rep*: spike-in normalisation (ChIP-Rx).

`chromatin-states`
- [ ] Ernst & Kellis 2012, *Nat Methods*: ChromHMM.
- [ ] Ernst et al. 2011, *Nature*: chromatin-state dynamics in nine cell types.

## DNA methylation

`methylation-biology`
- [ ] Okano et al. 1999, *Cell*: DNMT3A/3B and de novo methylation.
- [ ] Tahiliani et al. 2009, *Science*: TET1 converts 5mC to 5hmC.
- [ ] Smith & Meissner 2013, *Nat Rev Genet*: DNA methylation in development (review).
- [ ] Law & Jacobsen 2010, *Nat Rev Genet*: methylation patterns in plants and animals.

`methylation-bisulfite`
- [ ] Frommer et al. 1992, *PNAS*: bisulfite sequencing.
- [ ] Lister et al. 2009, *Nature*: the first base-resolution human methylomes.
- [ ] Krueger & Andrews 2011, *Bioinformatics*: Bismark.
- [ ] Vaisvila et al. 2021, *Genome Res*: EM-seq.

`methylation-differential`
- [ ] Hansen, Langmead & Irizarry 2012, *Genome Biol*: BSmooth.
- [ ] Feng, Conneely & Wu 2014, *Nucleic Acids Res*: DSS (beta-binomial).
- [ ] Akalin et al. 2012, *Genome Biol*: methylKit.

`methylation-arrays`
- [ ] Aryee et al. 2014, *Bioinformatics*: minfi.
- [ ] Du et al. 2010, *BMC Bioinformatics*: beta vs M values.
- [ ] Houseman et al. 2012, *BMC Bioinformatics*: cell-type composition from methylation.
- [ ] Jaffe & Irizarry 2014, *Genome Biol*: cell-type heterogeneity as a confounder in EWAS.

`epigenetic-clocks`
- [ ] Horvath 2013, *Genome Biol*: the pan-tissue epigenetic clock.
- [ ] Hannum et al. 2013, *Mol Cell*: blood methylation clock.
- [ ] Lu et al. 2019, *Aging*: GrimAge.

`longread-epigenomics`
- [ ] Simpson et al. 2017, *Nat Methods*: detecting 5mC with nanopore sequencing.
- [ ] Stergachis et al. 2020, *Science*: Fiber-seq.

## 3D genome

`genome-3d-intro`
- [ ] Dekker et al. 2002, *Science*: chromosome conformation capture (3C).
- [ ] Dixon et al. 2012, *Nature*: topological domains (TADs).
- [ ] Nora et al. 2012, *Nature*: spatial partitioning of the X-inactivation centre.
- [ ] Fudenberg et al. 2016, *Cell Rep*: loop extrusion.
- [ ] Sanborn et al. 2015, *PNAS*: loop extrusion and CTCF orientation.

`hic-processing`
- [ ] Lieberman-Aiden et al. 2009, *Science*: Hi-C; A/B compartments.
- [ ] Imakaev et al. 2012, *Nat Methods*: iterative correction (ICE).
- [ ] Abdennur & Mirny 2020, *Bioinformatics*: cooler.

`hic-features`
- [ ] Rao et al. 2014, *Cell*: kilobase-resolution Hi-C and loops.
- [ ] Rao et al. 2017, *Cell*: cohesin loss eliminates loop domains.
- [ ] Hsieh et al. 2015, *Cell*: Micro-C.

## Integration and interpretation

`integration-rnaseq`
- [ ] Fulco et al. 2019, *Nat Genet*: activity-by-contact (ABC) model.
- [ ] Gasperini et al. 2019, *Cell*: CRISPRi screen of enhancers.

`noncoding-variants`
- [ ] Maurano et al. 2012, *Science*: disease variants concentrate in regulatory DNA.
- [ ] Finucane et al. 2015, *Nat Genet*: partitioned heritability by functional annotation.
- [ ] Degner et al. 2012, *Nature*: DNase I sensitivity QTLs.

`epigenetic-inheritance`
- [ ] Heard & Martienssen 2014, *Cell*: transgenerational epigenetic inheritance, myths and mechanisms.
- [ ] Cubas, Vincent & Coen 1999, *Nature*: a natural epimutation (Linaria).
- [ ] Morgan et al. 1999, *Nat Genet*: agouti viable yellow.

`cancer-epigenetics`
- [ ] Figueroa et al. 2010, *Cancer Cell*: IDH mutations impair TET2.
- [ ] Schwartzentruber et al. 2012, *Nature*: histone H3.3 mutations in glioblastoma.
- [ ] Capper et al. 2018, *Nature*: methylation-based CNS tumour classification.
- [ ] Flavahan, Gaskell & Bernstein 2017, *Science*: epigenetic plasticity and cancer (review).

## Machine learning

`ml-basics`
- [ ] Libbrecht & Noble 2015, *Nat Rev Genet*: machine learning in genetics and genomics (review).

`ml-prediction-tasks`
- [ ] Eraslan et al. 2019, *Nat Rev Genet*: deep learning for genomics (review).
- [ ] Whalen et al. 2022, *Nat Rev Genet*: navigating pitfalls of ML in genomics.

`ml-cnn-dna`
- [ ] Alipanahi et al. 2015, *Nat Biotechnol*: DeepBind.
- [ ] Zhou & Troyanskaya 2015, *Nat Methods*: DeepSEA.
- [ ] Kelley, Snoek & Rinn 2016, *Genome Res*: Basset.

`ml-interpretation-variants`
- [ ] Avsec et al. 2021, *Nat Genet*: BPNet.
- [ ] Avsec et al. 2021, *Nat Methods*: Enformer.
- [ ] Sasse et al. 2023, *Nat Genet*: sequence-based models fail to explain expression between individuals.
