## WDLs for T2T Variants
This directory contains the WDL files used for large-scale short-read small-variant analysis.
### Data ingestion
- `wdls/download_aspera.wdl`: Downloads FASTQ files from the European Nucleotide Archive (ENA), given accession numbers
### Read alignment
- `wdls/t2t_alignment.wdl`: Given a reference FASTA file, sample name, paired-end FASTQ files, BWA index, and dedup distance (default = 100), performs alignment as described in [Aganezov, Yan, Soto, Kirsche, Zarate, *et al*. (2021)](https://doi.org/10.1101/2021.07.12.452063)
### Variant calling
- `wdls/haplotype_calling_chrom.wdl`: Given a reference FASTA (plus corresponding index and dict), a sample CRAM (plus corresponding index), the sample name, the sex of the sample, and whether or not to enable Spark (Spark functionality currently unavailable), performs GATK Haplotype Caller on a per-chromosome basis to generate per-chromosome VCF files for the sample as described in [Aganezov, Yan, Soto, Kirsche, Zarate, *et al*. (2021)](https://doi.org/10.1101/2021.07.12.452063)
- `wdls/joint_calling_interval.wdl`: Given a reference FASTA (plus corresponding index and dict), GenomicsDB tar file, interval, chromosome, start/end of interval, and margined start/end of interval, performs joint calling on the given interval across all samples included in the GenomicsDB tar file
- `wdls/t2t_concat_vcfs_chromosome.wdl`: Given an array of VCF files for different non-overlapping intervals, their corresponding index files, and the chromosome to which they all belong, concatenates the VCF files into a single chromosome-wide VCF file
- `wdls/t2t_genomics_db.wdl`: Given a GCP bucket, interval, chromosomes, margined start, and margined end, produces a tar file containing the genomics DB information for that interval across the samples specified in the sample map
- `wdls/variant_recalibration.wdl`: Performs variant recalibration on a VCF file with a number of databases as described in [Aganezov, Yan, Soto, Kirsche, Zarate, *et al*. (2021)](https://doi.org/10.1101/2021.07.12.452063)
### Downstream analysis
- `wdls/combine_mv_tables.wdl`: Given a file of concatenated Mendelian violation tables as produced by GATK VariantEval and a Mendelian violation table header, combines the tables into a summary Mendelian violation table with specific fields (# of variants, # of skipped variants, # of low-quality variants, # of Mendelian violations)
- `wdls/mendelian_analysis.wdl`: Given a reference FASTA (plus corresopnding index and dict), VCF file, sample list (one line per sample name), pedigree (in ped format), and population, performs an analysis of Mendelian violations using GATK VariantEval
- `wdls/local_ancestry_analysis.wdl`: Given a VCF file, windows BED file, ancestry BED file, population, and chromosome name, performs local ancestry analysis to identify which window corresponds to which ancestry/ancestries
### Miscellaneous utilities
- `wdls/bcftools_stats_region.wdl`: Given a VCF file in bgzip format, its corresponding index, and a BED file containing regions to examine, generates a stats file using `bcftools stats`
- `wdls/bcftools_stats_samples.wdl`: Given a multi-sample VCF file in bgzip format and a list of samples (one sample per line), generates a stats file using `bcftools stats`
- `wdls/bcftools_stats.wdl`: Given a VCF file in bgzip format, generates a stats file using `bcftools stats`
- `wdls/bcftools_view_subset.wdl`: Given a multi-sample VCF file, a list of samples (one sample per line), and a population name (for the output file name), generates a stats file using `bcftools stats`
- `wdls/bgzip_bcftools_index.wdl`: Given a VCF file, bgzips and indexes the file using bcftools
- `wdls/bgzip_tabix.wdl`: Given a VCF file, bgzips the file and indexes using tabix
- `wdls/concatenate_files.wdl`: Concatenates input files into a specified file name
- `wdls/filter_ac_variants.wdl`: Given a bgzipped VCF file, filters variants based on an allele count (AC); currently specifically looking for AC>1
- `wdls/filter_af_variants.wdl`: Given a bgzipped VCF file and list of samples, filters variants based on an allele frequency (AC); currently specifically looking for AF~=0.5
- `wdls/get_pass_records.wdl`: Given a bgzipped VCF file, filters the variants by whether or not they are annotated as PASS by GATK
- `wdls/gunzip.wdl`: Unzips a gzipped file
- `wdls/liftover.wdl`: Given a VCF file, chain file, and target reference FASTA (plus corresponding dict), lifts over VCF file from one reference to another
- `wdls/local_ancestry_analysis.wdl`: Given a CRAM file and regions BED file, subsets the CRAM file based on the regions given
- `wdls/subset_vcf_with_bed.wdl`: Given a VCF file in bgzip format, its corresponding index, and a BED file containing regions to examine, generates a stats file using `bcftools stats` as well as a VCF file corresponding to the region (and bgzipped version, plus an index)
- `wdls/tabix.wdl`: Performes `tabix` on a bgzipped VCF file