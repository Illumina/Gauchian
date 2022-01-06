# Gauchian: WGS-based GBA variant caller

Gauchian is a targeted variant caller for the GBA gene based on a whole-genome sequencing (WGS) BAM file. Gauchian uses a novel method to solve the problems caused by the high sequence similarity with the pseudogene paralog GBAP1 and is able to detect variants accurately in the Exons 9-11 homology region, such as large deletions or duplications between GBA and GBAP1, and GBAP1-like variants in GBA, including p.A495P, p.L483P, p.D448H, c.1263del, RecNciI, RecTL and c.1263del+RecTL. In addition to these challenging variants, Gauchian also calls known pathogenic or likely pathogenic GBA variants classified in ClinVar. Gauchian has been tested on Illumina WGS data with standard sequencing depth (>=30X). Gauchian does not work on targeted sequencing data. Please refer to our [preprint](https://www.medrxiv.org/content/10.1101/2021.11.12.21266253v1) for more details about the method.

## Installation

```bash
git clone https://github.com/Illumina/Gauchian
cd Gauchian
python3 setup.py install
```

## Running the program

```bash
gauchian --manifest MANIFEST_FILE \
         --genome [19/37/38] \
         --prefix OUTPUT_FILE_PREFIX \
         --outDir OUTPUT_DIRECTORY \
         --threads NUMBER_THREADS
```

The manifest is a text file in which each line should list the absolute path to an input WGS BAM/CRAM file. Full WGS BAM/CRAM files are recommended. If you would like to use a subsetted bamlet, please subset using region files in gauchian/data/GBA_region_*.bed.

For CRAM input, it’s suggested to provide the path to the reference fasta file with `--reference` in the command.

## Interpreting the output

The program produces a .tsv file in the directory specified by --outDir.
The fields are explained below:

| Fields in tsv                            | Explanation                                                                    |
|:-----------------------------------------|:-------------------------------------------------------------------------------|
| Sample                                   | Sample name                                                                    |
| is_biallelic(GBAP1-like_variant_exon9-11)| Whether the sample is called as biallelic for GBAP1-like variants in exon9-11  |
| is_carrier(GBAP1-like_variant_exon9-11)  | Whether the sample is called as a carrier for GBAP1-like variants in exon9-11  |
| CN(GBA+GBAP1)                            | Total copy number of GBA+GBAP1                                                 |
| deletion_breakpoint_in_GBA               | Whether the deletion breakpoint is in GBA gene if a deletion exists            |
| GBAP1-like_variant_exon9-11              | GBAP1-like variants called in exon9-11, two alleles separated by /             |
| other_unphased_variants                  | Other variants called (non-GBAP1-like variants or variants outside of exon9-11)|

A .json file is also produced that contains more information about each sample.

| Fields in json    | Explanation                                                                       |
|:------------------|:----------------------------------------------------------------------------------|
| Coverage_MAD      | Median absolute deviation of depth, measure of sample quality                     |
| Median_depth      | Sample median depth                                                               |
| deletion_CN       | CN of the unique region between GBA and GBAP1. This value plus 2 is the total CN  |
| deletion_CN_raw   | Raw normalized depth of the unique region between GBA and GBAP1                   |
| variant_raw_count | Supporting reads for each variant                                                 |
| snp_call          | GBA copy number call at GBA/GBAP1 differentiating sites                           |
| snp_raw           | Raw GBA copy number at GBA/GBAP1 differentiating sites                            |
| haplotypes        | Summary of haplotypes assembled across GBA/GBAP1 differentiating sites in Exon9-11|
