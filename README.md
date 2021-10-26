# Variant analysis scripts and workflows for:

[A complete reference genome improves analysis of human genetic variation](https://www.biorxiv.org/content/10.1101/2021.07.12.452063v1)

Aganezov, S\*, Yan, SM\*, Soto, DC\*, Kirsche, M\*, Zarate, S\*, Avdeyev, P, Taylor, DJ, Shafin, K, Shumate, A, Xiao, C, Wagner, J, McDaniel, J, Olson, ND, Sauria, MEG, Vollger, MR, Rhie, A, Meredith, M, Martin, S, Lee, J, Koren, S, Rosenfeld, J, Paten, B, Layer, R, Chin, CS, Sedlazeck, FJ, Hansen, NF, Miller, DE, Phillippy, AM, Miga, KM, McCoy, RC†, Dennis, MY†, Zook, JW†, Schatz, MC† (2021) bioRxiv. [doi: https://doi.org/10.1101/2021.07.12.452063](https://www.biorxiv.org/content/10.1101/2021.07.12.452063v1)

## Primary analysis
- [1000genomes_metadata](1000genomes_metadata): how/where to download metadata for 1000 genomes samples
- [figures](figure): individual panels from Figures 2 and 3 along with related supplemental figures
- [population-variant-stats-per-sample](population-variant-stats-per-sample/): count of variants per sample in 1KGP across population/superpopulation as plotted in Figure 2 and related Supplemental Figures
- [svs.ipynb](svs.ipynb) Analysis of long read mapping and SV calls
- [superpopulation-af](superpopulation-af/): Allele frequency of variants in each 1KGP superpopulation as plotted in Figure 2
- [wdls](wdls): Workflows for processing the 1000 genomes samples in [AnVIL](http://anvilproject.org)


## Exploratory analysis
- [chromosome_repeats](chromosome_repeats): make plots of uniqueness/repetitiveness of each chromosome
- [chromosome_variants](chromosome_variants): make plots of variants per chromosome
- [local_ancestry](local_ancestry): make plots of number of variants in regions with different local ancestries
- [novel_region_af](novel_region_af): Allele frequency in novel regions (previously unresolved in GRCh38)
- [novel_regions](novel_regions): Notes on novel regions of T2T-CHM13 (previously unresolved by GRCh38)
- [population-variant-stats](population-variant-stats): Summary of variants per population
- [samtools_stats_chm13](samtools_stats_chm13): summary plots of mapping rates by population/superpopulation for 1000 genomes samples
