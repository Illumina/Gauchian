#!/usr/bin/env python3
#
# Gauchian: GBA variant caller
# Copyright 2021 Illumina, Inc.
# All rights reserved.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is licensed under the terms of the Polyform strict license
#
# ***As far as the law allows, the software comes as is, without
# any warranty or condition, and the licensor will not be liable
# to you for any damages arising out of these terms or the use
# or nature of the software, under any kind of legal claim.***
#
# You should have received a copy of the PolyForm Strict License 1.0.0
# along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
#
#


from pprint import pprint
from collections import namedtuple, Counter
from .snp_count import get_supporting_reads
from .copy_number_call import call_cn_snp
from .haplotype import (
    get_haplotypes_from_bam,
    extract_hap,
    filter_hap,
    call_hap_cn,
    get_block_read_count,
)


class Phasing:
    """
    A class for phasing haplotypes across differentiating sites.
    """
    var_call = namedtuple(
        "var_call",
        "total_cn is_biallelic is_carrier depth_per_site deletion_bp_in_gene \
            full_haplotypes variants_allele1 variants_allele2 variant_details"
    )
    assembled_haplotypes = namedtuple(
        "assembled_haplotypes", "full_haplotypes hap_5p hap_3p hap_mid hap3p_merged"
    )
    var_call_per_haplotype = namedtuple(
        "var_call_per_haplotype", "haplotype variant_name block_result is_variant"
    )
    block_read_support = namedtuple(
        "block_read_support", "read_count flankings read_count_ordered likelihood"
    )

    CN_likelihood_threshold = 0.7
    CN_likelihood_threshold_loose = 0.45
    CN_likelihood_threshold_stringent = 0.95

    def __init__(self):
        self.variant_sites = []
        self.trusted_hap_threshold = []
        self.trusted_hap_switch_point = []
        self.variant_names = []
        self.variant_names_in_one = []
        self.flanking_site = {}
        self.bamfile = None
        self.base_db = None
        self.total_cn = None
        self.num_sites = None
        self.depth_calls = None
        self.reference = None
        self.haplotype_per_read = None
        self.gene1_read_count = []
        self.gene2_read_count = []
        self.deletion_bp_in_gene = None
        self.variants_details = []
        self.is_carrier = False
        self.is_biallelic = False
        self.variants_called_allele1 = []
        self.variants_called_allele2 = []
        self.depth_base_cn = None
        self.candidate_gene_hap = []

    def set_par(self, bamfile, base_db, total_cn, reference):
        """
        Establish parameter values.
        Parameters:
            bamfile: bam file
            base_db: differentiating site info
            total_cn: summed CN of paralogs
            reference: reference file
        """
        self.bamfile = bamfile
        self.base_db = base_db
        self.total_cn = total_cn
        self.reference = reference
        self.num_sites = len(base_db.dsnp1)
        self.gene1_read_count, self.gene2_read_count = get_supporting_reads(
            self.bamfile,
            self.base_db.dsnp1,
            self.base_db.dsnp2,
            self.base_db.nchr,
            self.base_db.dindex
        )
        self.depth_calls = call_cn_snp(
            self.total_cn,
            self.gene1_read_count,
            self.gene2_read_count
            )
        self.haplotype_per_read = get_haplotypes_from_bam(
            self.bamfile,
            self.base_db,
            range(self.num_sites),
            reference=self.reference
            )
        if total_cn < 4:
            self.deletion_bp_in_gene = False

    def carrier_from_assembled_haps(self, var_sites_haps):
        """
        Test if a sample is a carrier based on assembled haplotypes.
        Parameters:
            var_sites_haps (list of str): haplotypes only at variant sites
        Returns:
            bool
        """
        # Limit to total_cn=3
        # All haplotypes assembled. All unique.
        if len(var_sites_haps) == self.total_cn and self.total_cn == 3:
            gene_cn = [0] * len(self.variant_sites)
            for hap in var_sites_haps:
                for i in range(len(hap)):
                    if hap[i] == "1":
                        gene_cn[i] += 1
            # CN(gene) is 1 in fully assembled haplotypes.
            if 1 in gene_cn:
                return True
        return False

    def get_variants(self, hap):
        """
        Get variant names given haplotype
        Parameters:
            hap (str): an assembed haplotype
        Returns:
            variants (list of str): names of variants carried by the haplotype
        """
        variants = []
        for i in range(len(hap)):
            if hap[i] == "2" and i in self.variant_sites:
                variant_index = self.variant_sites.index(i)
                variants.append(self.variant_names[variant_index])
        return variants

    def assess_assembled_haps(self, full_haplotypes):
        """
        Check if haplotypes have been completely assembled.
        Parameters:
            full_haplotypes (list of str): assembled haplotypes
        Returns:
            fully_assembled (bool)
            var_sites_haps (list of str): haplotypes only at variant sites
        """
        var_sites_haps = []
        fully_assembled = True
        for hap in full_haplotypes:
            if 'x' in hap:
                fully_assembled = False
            else:
                sites_of_interest = ''.join(
                    [hap[i] for i in self.variant_sites]
                    )
                var_sites_haps.append(sites_of_interest)
        return fully_assembled, var_sites_haps

    def filter_haps_get_var_sites(self, full_haplotypes, trust_hap=False):
        """
        From haplotypes, output variant sites (switching points)
        Parameters:
            full_haplotypes (list of str): assembled haplotypes
            trust_hap (bool): whether to trust trusted haplotypes
        Returns:
            var_index_in_haps (dict of str:list): index of candidate
            variant sites on each haplotype
        """
        var_index_in_haps = {}
        for hap in full_haplotypes:
            count_2 = hap.count("2")
            count_1 = self.depth_calls.count(1)
            count_0 = self.depth_calls.count(0)
            count_none = self.depth_calls.count(None)
            # this is the criteria for selecting candidate haplotypes,
            # comparing haplotypes against the depth calls
            if (
                count_2 <= count_1 + count_0 + count_none + 1 or
                (trust_hap and hap in self.trusted_hap)
            ):
                for var_index in range(self.num_sites):
                    if hap[var_index] == "2":
                        # find the first variant site past the switching point
                        if (var_index > 0 and hap[var_index - 1] == "1"):
                            n = var_index
                            while n < self.num_sites:
                                if n in self.variant_sites and hap[n] == "2":
                                    var_index_in_haps.setdefault(hap, [])
                                    if n not in var_index_in_haps[hap]:
                                        var_index_in_haps[hap].append(n)
                                    break
                                n += 1
                        if (
                            var_index + 1 < self.num_sites and
                            hap[var_index + 1] == "1"
                            ):
                            n = var_index
                            while n >= 0:
                                if n in self.variant_sites and hap[n] == "2":
                                    var_index_in_haps.setdefault(hap, [])
                                    if n not in var_index_in_haps[hap]:
                                        var_index_in_haps[hap].append(n)
                                    break
                                n -= 1
        return var_index_in_haps

    def get_depth_for_block(self, flanking_sites_per_var, var_index_in_flanking):
        """
        from single switching point, check if read counts
        support that there is only one copy of WT gene
        Parameters:
            flanking_sites_per_var (list of int): sites flanking a variant to check depth
            var_index_in_flanking (int): index of variant site in haplotypes
        Returns:
            named tuple (block_read_support) summarizing depth info at this site
        """
        switch_site_haps_count = extract_hap(
                self.haplotype_per_read, range(flanking_sites_per_var[0], flanking_sites_per_var[1] + 1)
            )
        switch_site_haps_count = filter_hap(switch_site_haps_count)
        if (
            self.total_cn == len(switch_site_haps_count) and
            min(switch_site_haps_count.values()) >= 2
        ):
            return(
                self.block_read_support(switch_site_haps_count, flanking_sites_per_var, None, 1)
                )
        if self.total_cn >= 2:
            block_read_count = get_block_read_count(switch_site_haps_count, var_index_in_flanking)
            if block_read_count is not None:
                if block_read_count == [0]:
                    return(
                        self.block_read_support(
                            switch_site_haps_count,
                            flanking_sites_per_var,
                            block_read_count,
                            1
                            )
                        )
                cn_call = call_hap_cn(block_read_count, self.total_cn)
                return(
                    self.block_read_support(
                        switch_site_haps_count,
                        flanking_sites_per_var,
                        block_read_count,
                        round(cn_call, 4)
                        )
                    )
        return(self.block_read_support(switch_site_haps_count, flanking_sites_per_var, None, None))

    def get_depth_for_blocks(self, good_var_indices, cn_threshold):
        """
        from a set of switching points, check if read counts
        support that there is only one copy of WT gene
        Parameters:
            good_var_indices (list of int): indices of variant sites on a haplotype
            cn_threshold (float): likelihood cutoff
        Returns:
            found_variant_in_short_hap (bool): whether a variant is found based on depth
            blocks (list of get_depth_for_block named tuple): summary of depth info
        """
        found_variant_in_short_hap = False
        blocks = []
        for var_index in good_var_indices:
            for flanking_sites_per_var in self.flanking_site[var_index]:
                var_index_in_flanking = None
                sites_range = list(range(flanking_sites_per_var[0], flanking_sites_per_var[1] + 1))
                if var_index in sites_range:
                    var_index_in_flanking = sites_range.index(var_index)
                block = self.get_depth_for_block(flanking_sites_per_var, var_index_in_flanking)
                blocks.append(block)
                if block.likelihood is not None and block.likelihood > cn_threshold:
                    found_variant_in_short_hap = True
        return found_variant_in_short_hap, blocks

    def compare_haplotype_against_depth(self, hap):
        """
        In rare cases there is more than one candidate haplotype.
        Pick the one that best explains the depth values
        Parameters:
            hap (list of str): a haplotype
        Returns:
            (float) ratio of matches compared to depth paterns
        """
        assert len(hap) == len(self.depth_calls)
        nmatch = 0
        nsites = len(hap)
        for i in range(nsites):
            if self.depth_calls[i] is not None:
                if (
                    (hap[i] == '2' and self.depth_calls[i] == 1) or
                    (hap[i] == '1' and self.depth_calls[i] >= 2)
                ):
                    nmatch += 1
        return nmatch/nsites

    def check_complete_loss(self, var_index_in_haps):
        """completely lose one gene, or one long conversion"""
        if var_index_in_haps == {}:
            gene1_total = sum(self.gene1_read_count)
            gene2_total = sum(self.gene2_read_count)
            combined_cn_call = call_cn_snp(
                self.total_cn,
                [gene1_total],
                [gene2_total],
                0.8
                )
            if combined_cn_call[0] is not None and combined_cn_call[0] == 0 and self.total_cn < 4:
                self.deletion_bp_in_gene = True
            if combined_cn_call[0] is not None and combined_cn_call[0] == 1:
                count_1 = self.depth_calls.count(1)
                count_none = self.depth_calls.count(None)
                if (
                    count_1 + count_none >= len(self.depth_calls) - 1 and
                    count_1 > count_none
                ):
                    self.is_carrier = True
                    self.variants_details.append(
                        self.var_call_per_haplotype(
                            None,
                            self.variant_names_in_one,
                            None,
                            True
                            )
                        )
                    if self.total_cn < 4:
                        self.deletion_bp_in_gene = True

    def assess_compound_het(self, full_haplotypes):
        """assess compound het"""
        fully_assembled, var_sites_haps = self.assess_assembled_haps(full_haplotypes)
        if fully_assembled is True:
            # Among all haplotypes,
            # there is not a copy that are WT in all positions
            if "1"*len(self.variant_sites) not in var_sites_haps:
                return True
        return False

    def assess_homozygous(self):
        """assess if there are homozygous sites"""
        homozygous_depth = False
        for i in self.variant_sites:
            if self.depth_calls[i] is not None and self.depth_calls[i] == 0:
                homozygous_depth = True
                self.variants_called_allele1.append(self.variant_names[self.variant_sites.index(i)])
                self.variants_called_allele2.append(self.variant_names[self.variant_sites.index(i)])
        return homozygous_depth

    def pick_best_variant(self, compound_het):
        """
        In rare scenarios where sample is not compound het but
        there are more than one variant haplotypes,
        pick one that best fits depth calls
        Parameters:
            compound_het (bool): whether sample is compound het
        """
        unique_vars = []
        unique_haps = []
        for var in self.variants_details:
            if var.is_variant:
                if var.variant_name not in unique_vars:
                    unique_vars.append(var.variant_name)
                if var.haplotype not in unique_haps:
                    unique_haps.append(var.haplotype)
        if len(unique_haps) > 1 and compound_het is False:
            hap_score = []
            var_list = []
            for var in self.variants_details:
                if var.is_variant:
                    var_list.append(var)
                    hap_score.append(self.compare_haplotype_against_depth(var.haplotype))
            best_vars = []
            max_score = max(hap_score)
            for i,score in enumerate(hap_score):
                if score == max_score:
                    best_vars.append(var_list[i].variant_name)
            best_vars = sorted(best_vars, key=lambda kv: len(kv))
            best_var = best_vars[0]
            variants_called_pick_one = []
            for var in self.variants_details:
                if var.variant_name == best_var:
                    variants_called_pick_one.append(var)
                    unique_vars = [var.variant_name]
                else:
                    variants_called_pick_one.append(var._replace(is_variant=False))
            self.variants_details = variants_called_pick_one
        unique_vars = sorted(unique_vars, key=lambda kv: len(kv))
        return unique_vars

    def get_candidate_gene_hap(self, full_haplotypes):
        """
        Get haplotypes that are (roughly) likely from the gene based
        on the number of gene bases
        """
        for haplotype in full_haplotypes:
            if haplotype.count("1") > len(haplotype)/2:
                self.candidate_gene_hap.append(haplotype)

    def get_base_cn(self, start, end):
        """
        Get (roughly) the number of gene haplotypes from depth calls
        Updates depth_base_cn
        Parameters:
            start (int): start position to count depth calls
            end (int): end position to count depth calls
        """
        clean_depth = []
        for i, site in enumerate(self.depth_calls):
            if start <= i <= end:
                if site is not None:
                    clean_depth.append(site)
        clean_depth_counter = sorted(
            Counter(clean_depth).items(), key=lambda kv: kv[1], reverse=True
        )
        if (
            clean_depth_counter[0] is not None and
            clean_depth_counter[0][1] >= len(self.depth_calls)/2
        ):
            self.depth_base_cn = clean_depth_counter[0][0]

    def check_deletion_bp_in_gene(self, var_index_in_haps):
        """
        In case of a deletion,
        check whether the deletion breakpoint is within the gene
        """
        return False

    def assess_hap(self, hap, var_index_in_haps, full_haplotypes):
        """
        Assess whether a haplotype is a variant haplotype that makes
        the sample a carrier i.e. having only one copy of the WT gene
        """
        return False

    def assemble_haplotypes(self, debug=False):
        """Assemble haplotypes that consist of differentiating sites"""
        return self.assembled_haplotypes([], None, None, None, None)

    def call_variants_from_haplotype(self, debug=False):
        """
        Call recombinant variants by phasing through differentiating sites
        """
        if debug is True:
            pprint(self.haplotype_per_read)

        # assemble haplotypes
        haps_assembled = self.assemble_haplotypes(debug=debug)
        full_haplotypes = haps_assembled.full_haplotypes
        self.get_candidate_gene_hap(full_haplotypes)

        homozygous_depth = self.assess_homozygous()
        compound_het = self.assess_compound_het(full_haplotypes)

        var_index_in_haps = self.filter_haps_get_var_sites(full_haplotypes)
        self.check_deletion_bp_in_gene(var_index_in_haps)
        # assess each haplotype to see if they makes the sample a carrier,
        # i.e. having one copy of the WT gene
        for hap in var_index_in_haps:
            self.assess_hap(hap, var_index_in_haps, full_haplotypes)

        # completely lose one gene, or one long conversion across all sites
        self.check_complete_loss(var_index_in_haps)

        unique_vars = self.pick_best_variant(compound_het)
        self.is_biallelic = homozygous_depth or (compound_het and len(unique_vars) >= 2)
        if self.is_biallelic:
            self.is_carrier = False
        if unique_vars != []:
            self.variants_called_allele1 += unique_vars[0]
            if len(unique_vars) > 1:
                self.variants_called_allele2 += unique_vars[1]

        return self.var_call(
            self.total_cn,
            self.is_biallelic,
            self.is_carrier,
            '_'.join([str(a) for a in self.depth_calls]),
            self.deletion_bp_in_gene,
            haps_assembled._asdict(),
            list(set(self.variants_called_allele1)),
            list(set(self.variants_called_allele2)),
            [a._asdict() for a in self.variants_details]
            )
