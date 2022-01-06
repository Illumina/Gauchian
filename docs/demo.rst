============
Program Demo
============

1. Download a WGS bam file::

    $ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239883/NA20815.final.cram
    $ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239883/NA20815.final.cram.crai

2. Make a manifest file::

    $ echo "$PWD/NA20815.final.cram" > manifest.txt

3. Run Gauchian, which takes about one to two minutes with single thread::

    $ gauchian -m manifest.txt -g 38 -o out -p testgba

4. Check the output file out/testgba.tsv

============= ========================================= ======================================= ============= ========================== =========================== =======================
Sample        is_biallelic(GBAP1-like_variant_exon9-11) is_carrier(GBAP1-like_variant_exon9-11) CN(GBA+GBAP1) deletion_breakpoint_in_GBA GBAP1-like_variant_exon9-11 other_unphased_variants
============= ========================================= ======================================= ============= ========================== =========================== =======================
NA20815.final False                                     True                                    4             N/A                        L483P/                      None
============= ========================================= ======================================= ============= ========================== =========================== =======================
