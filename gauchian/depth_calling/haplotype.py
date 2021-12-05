#!/usr/bin/env python3
#
# Gauchian: GBA variant caller
# Copyright 2021 Illumina, Inc.
# All rights reserved.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#


from scipy.stats import poisson
import re
import pysam
from collections import namedtuple
from .snp_count import reverse_complement

ERROR_RATE = 0.01
haplotype_groups = namedtuple(
    "haplotype_groups",
    "unique nonunique"
)
clip_start_pat = re.compile(r"^(\d+)S")
clip_end_pat = re.compile(r"(\d+)S$")


def get_NM(ltag):
    """Get NM tag for a read alignment"""
    for tag in ltag:
        if tag[0] == "NM":
            return tag[1]
    return None


def get_mismatches_from_softclip(read_alignment, ref, min_base_quality=13):
    """
    Get number of mismatches in soft-clipped sequences, considering base quality
    Sometimes a read is soft-clipped due to low base quality, so the size of the clip
    is not necessarily equal to the length of a mismatch.
    We want to keep such reads so the unclipped part is informative.
    """
    mismatch = 0
    if ref is None:
        return mismatch
    ref_pos = read_alignment.reference_start
    ref_end = read_alignment.reference_end
    cigar = read_alignment.cigarstring
    read_len = len(read_alignment.query_sequence)
    left_clip = sum([int(a) for a in clip_start_pat.findall(cigar)])
    right_clip = sum([int(a) for a in clip_end_pat.findall(cigar)])
    if left_clip > 0:
        for i in range(left_clip):
            if read_alignment.query_qualities[left_clip - (i+1)] >= min_base_quality:
                pos = ref_pos - (i+1)
                ref_base = ref.fetch(read_alignment.reference_name, pos, pos+1).upper()
                read_base = read_alignment.query_sequence[left_clip - (i+1)].upper()
                if ref_base != read_base:
                    mismatch += 1
    if right_clip > 0:
        for i in range(right_clip):
            if read_alignment.query_qualities[read_len - (right_clip - i)] >= min_base_quality:
                pos = ref_end + i
                ref_base = ref.fetch(read_alignment.reference_name, pos, pos+1).upper()
                read_base = read_alignment.query_sequence[read_len - (right_clip - i)].upper()
                if ref_base != read_base:
                    mismatch += 1
    return mismatch


def passing_read(pileupread, ref=None, long_read=False):
    """Return whether a read passes filter."""
    nm = get_NM(pileupread.alignment.tags)
    start_clip = sum([int(a) for a in clip_start_pat.findall(pileupread.alignment.cigarstring)])
    end_clip = sum([int(a) for a in clip_end_pat.findall(pileupread.alignment.cigarstring)])
    clip = start_clip + end_clip
    clip_mismatch = 0
    if ref is not None and clip > 15:
        clip_mismatch = get_mismatches_from_softclip(pileupread.alignment, ref)
    return (
        not pileupread.is_del
        and not pileupread.is_refskip
        and pileupread.alignment.is_secondary == 0
        and pileupread.alignment.is_supplementary == 0
        and pileupread.alignment.is_duplicate == 0
        and (nm is None or (nm is not None and nm <= 6) or long_read)
        and (
            (ref is None and clip <= 15)
            or (ref is not None and clip_mismatch <= 10)
        )
    )


def get_haplotypes_from_bam(bamfile_handle, base_db, target_positions, rc=None, reference=None):
    """
    Get haplotypes for each read at a set of target positions
    Parameters:
        bamfile_handle
        base_db (named tuple): differentiating site info
        target_positions (list of int): positions of interest
        rc (bool): reverse complement if any gene is on the reverse strand
        reference (str): reference fasta
    Returns:
        haplotype_per_read (dict of str : str): bases of each read at differentiating sites
    """
    dread = {}
    dread = get_bases_per_read(bamfile_handle, base_db, target_positions, rc=rc, reference=reference)
    base1, base2 = get_base1_base2(base_db, target_positions)
    haplotype_per_read = get_haplotype_per_read(dread, base1, base2)
    return haplotype_per_read


def get_haplotypes_from_bam_single_region(bamfile_handle, base_db, target_positions, rc=None, reference=None):
    """
    Get haplotypes for each read at a set of target positions, from only the gene region
    Parameters:
        bamfile_handle
        base_db (named tuple): differentiating site info
        target_positions (list of int): positions of interest
        rc (bool): reverse complement if any gene is on the reverse strand
        reference (str): reference fasta
    Returns:
        haplotype_per_read (dict of str : str): bases of each read at differentiating sites
    """
    dread = {}
    dread = get_bases_per_read(
        bamfile_handle, base_db, target_positions, region=0, min_mapq=10, rc=rc, reference=reference
    )
    base1, base2 = get_base1_base2(base_db, target_positions)
    haplotype_per_read = get_haplotype_per_read(dread, base1, base2)
    return haplotype_per_read


def get_bases_per_read(
    bamfile_handle, base_db, target_positions, region=None, min_mapq=0, rc=None, reference=None
):
    """
    Go through reads and get bases at sites of interest
    Paramters:
        bamfile_handle
        base_db (named tuple): differentiating site info
        target_positions (list of int): positions of interest
        region (int): analyze both regions (gene and pseudogene) or just one.
        min_mapq (int): mapq cutoff
        rc (bool): reverse complement if any gene is on the reverse strand
        reference (str): reference fasta
    Returns:
        dread (dict of str:list): bases at each site on individual reads,
        before converting to 1s and 2s
    """
    ref = None
    if reference is not None:
        ref = pysam.Fastafile(reference)
    dread = {}
    nchr = base_db.nchr
    dindex = base_db.dindex
    dsnps = [base_db.dsnp1, base_db.dsnp2]
    if region is not None:
        if region == 0:
            dsnps = [base_db.dsnp1]
        elif region == 1:
            dsnps = [base_db.dsnp2]
    for i in range(len(dsnps)):
        dsnp = dsnps[i]
        for snp_position_ori in dsnp:
            dsnp_index = dindex[snp_position_ori]
            snp_position = int(snp_position_ori.split("_")[0])
            if dsnp_index in target_positions:
                for pileupcolumn in bamfile_handle.pileup(
                    nchr,
                    snp_position - 1,
                    snp_position + 1,
                    truncate=True,
                    stepper="nofilter",
                    ignore_overlaps=False,
                    ignore_orphan=False,
                    min_base_quality=13,
                ):
                    site_position = pileupcolumn.pos + 1
                    if site_position == snp_position:
                        reg1_allele, reg2_allele = dsnp[snp_position_ori].split("_")
                        for read in pileupcolumn.pileups:
                            if (
                                passing_read(read, ref=ref)
                                and read.alignment.mapping_quality >= min_mapq
                            ):
                                read_name = read.alignment.query_name
                                read_seq = read.alignment.query_sequence
                                start_pos = read.query_position
                                end_pos = start_pos + min(
                                    len(reg1_allele), len(reg2_allele)
                                )
                                if end_pos < len(read_seq):
                                    hap = read_seq[start_pos:end_pos]
                                    if i == rc:
                                        hap = reverse_complement(hap)
                                    if read_name not in dread:
                                        dread.setdefault(read_name, {})
                                        for pos in target_positions:
                                            dread[read_name].setdefault(pos, None)
                                    if dread[read_name][dsnp_index] not in [None, hap]:
                                        dread[read_name][dsnp_index] = None
                                    dread[read_name][dsnp_index] = hap
    return dread


def get_haplotype_per_read(dread, base1, base2):
    """
    Translate bases into haplotypes, labeled as 1 or 2
    Parameters:
        dread (dict of str:list): bases at each site on individual reads
        base1 (list of str): gene1 bases
        base2 (list of str): gene2 bases
    Returns:
        haplotype_per_read (dict of str : str): bases of each read at differentiating sites
    """
    dread_hap = {}
    for read in dread:
        read_bases = dread[read]
        assert len(read_bases) == len(base1) == len(base2)
        pos_list = ["x"] * len(read_bases)
        for i, pos in enumerate(read_bases):
            base = read_bases[pos]
            if base is not None:
                for allele in base1[i].split(","):
                    if base == allele.upper():
                        pos_list[i] = "1"
                for allele in base2[i].split(","):
                    if base == allele.upper():
                        pos_list[i] = "2"
        dread_hap.setdefault(read, "".join(pos_list))
    return dread_hap


def get_base1_base2(base_db, target_positions):
    """
    Get expected bases corresponding to different alleles/paralogs
    at target positions
    Parameters:
        base_db (named tuple): differentiating site info
        target_positions (list of int): positions of interest
    Returns:
        list of gene1 bases, list of gene2 bases
    """
    base1 = [None] * len(base_db.dsnp1)
    base2 = [None] * len(base_db.dsnp1)
    dindex = base_db.dindex
    for pos in base_db.dsnp1:
        dsnp_index = dindex[pos]
        if dsnp_index in target_positions:
            index = int(pos.split("_")[1])
            allele1, allele2 = base_db.dsnp1[pos].split("_")
            base1[index] = allele1
            base2[index] = allele2
    return [base1[i] for i in target_positions], [base2[i] for i in target_positions]


def extract_hap(haplotype_per_read, positions_to_extract):
    """
    Extract haplotypes at given positions
    Paramters:
        haplotype_per_read (dict of str:str): bases of each read at differentiating sites
        positions_to_extract (list of int): positions of interest
    Returns:
        hap_count (dict of str:list): haplotypes and a list of supporting reads (as 1)
    """
    hap_count = {}
    for read_name in haplotype_per_read:
        hap = haplotype_per_read[read_name]
        hap_base = [hap[pos] for pos in positions_to_extract]
        if "x" not in hap_base:
            hap_count.setdefault("".join(hap_base), []).append(1)
    return hap_count


def group_haps(haplotype_per_read, hap_list):
    """
    Assign reads to partial haplotype blocks so haplotype blocks can be extended
    Parameters:
        haplotype_per_read (dict of str:str): bases of each read at differentiating sites
        hap_list (list of str): haplotypes to be extended
    Returns:
        (named tuple of dict of str:list): reads grouped with each haplotype,
        for uniquely and nonuniquely associated reads
    """
    matching_haplotype_groups = {}
    matching_haplotype_groups_nonunique = {}
    for read_name in haplotype_per_read:
        read_hap = haplotype_per_read[read_name]
        matching_haplotypes = []
        for haplotype_to_extend in hap_list:
            match, mismatch, extend = compare_two_haps(haplotype_to_extend, read_hap)
            if mismatch == 0 and match > 0 and extend >= 0:
                matching_haplotypes.append(haplotype_to_extend)
                matching_haplotype_groups_nonunique.setdefault(haplotype_to_extend, []).append(read_hap)
        # unique match
        if len(matching_haplotypes) == 1:
            matching_haplotype_groups.setdefault(matching_haplotypes[0], []).append(read_hap)
    return haplotype_groups(matching_haplotype_groups, matching_haplotype_groups_nonunique)


def extend_hap_5p(matching_haplotype_groups, min_read=None):
    """
    Extend a haplotype on 5p based on reads
        Parameters:
        matching_haplotype_groups (named tuple of dict of str:list):
        reads grouped with each haplotype to be extended, for uniquely
        and nonuniquely associated reads
        min_read (int): minimum read count for extension
    Returns:
        new_hap_list (list of str): left extended haplotypes
    """
    new_hap_list = []
    for haplotype_to_extend in matching_haplotype_groups.unique:
        if "x" in haplotype_to_extend:
            new_hap = haplotype_to_extend
            haps = matching_haplotype_groups.unique[haplotype_to_extend]
            pos_x = 0
            # find the first non-"x" position to start extending
            while pos_x < len(new_hap):
                if new_hap[pos_x] != "x":
                    break
                pos_x += 1
            if pos_x == 0:
                new_hap_list.append(haplotype_to_extend)
            else:
                extended = False
                pos = pos_x
                i = 0
                # going to the left site by site and
                # get consensus of the base from matching read haplotypes
                while pos > 0:
                    pos = pos - 1
                    i += 1
                    bases_at_site = [a[pos] for a in haps if a[pos] != "x"]
                    if bases_at_site != []:
                        uniq_bases = []
                        for uniq_base in list(set(bases_at_site)):
                            min_read_site = min_read
                            if min_read is None:
                                min_read_site = get_read_count_threshold(bases_at_site)
                            if bases_at_site.count(uniq_base) >= min_read_site:
                                uniq_bases.append(uniq_base)
                        for uniq_base in uniq_bases:
                            new_extended_haplotype = "x" * pos + uniq_base + new_hap[(pos + 1) :]
                            new_hap_list.append(new_extended_haplotype)
                            extended = True
                        break
                if extended is False and matching_haplotype_groups.nonunique:
                    pos = pos_x - 1
                    next_pos_nonuniques = [
                        a[pos] for a in matching_haplotype_groups.nonunique[haplotype_to_extend] if a[pos] != "x"
                        ]
                    if len(set(next_pos_nonuniques)) == 1 and len(next_pos_nonuniques) > 4:
                        new_extended_haplotype = "x" * pos + next_pos_nonuniques[0] + new_hap[(pos + 1) :]
                        new_hap_list.append(new_extended_haplotype)
                        extended = True
                # failed to extend
                if extended is False:
                    new_hap_list.append(haplotype_to_extend)
        else:
            new_hap_list.append(haplotype_to_extend)
    return new_hap_list


def extend_hap_3p(matching_haplotype_groups, min_read=None):
    """
    Extend a haplotype on 3p based on reads
    Parameters:
        matching_haplotype_groups (named tuple of dict of str:list):
        reads grouped with each haplotype to be extended, for uniquely
        and nonuniquely associated reads
        min_read (int): minimum read count for extension
    Returns:
        new_hap_list (list of str): right extended haplotypes
    """
    new_hap_list = []
    for haplotype_to_extend in matching_haplotype_groups.unique:
        if "x" in haplotype_to_extend:
            new_hap = haplotype_to_extend
            hap_len = len(new_hap)
            haps = matching_haplotype_groups.unique[haplotype_to_extend]
            pos_x = hap_len - 1
            # find the first non-"x" position from the 3p end to start extending
            while pos_x >= 0:
                if new_hap[pos_x] != "x":
                    break
                pos_x -= 1
            if pos_x == hap_len - 1:
                new_hap_list.append(haplotype_to_extend)
            else:
                extended = False
                pos = pos_x
                i = 0
                while pos < hap_len - 1:
                    pos = pos + 1
                    i += 1
                    bases_at_site = [a[pos] for a in haps if a[pos] != "x"]
                    if bases_at_site != []:
                        uniq_bases = []
                        for uniq_base in list(set(bases_at_site)):
                            min_read_site = min_read
                            if min_read is None:
                                min_read_site = get_read_count_threshold(bases_at_site)
                            if bases_at_site.count(uniq_base) >= min_read_site:
                                uniq_bases.append(uniq_base)
                        for uniq_base in uniq_bases:
                            new_extended_haplotype = new_hap[:pos] + uniq_base + "x" * (hap_len - 1 - pos)
                            new_hap_list.append(new_extended_haplotype)
                            extended = True
                        break
                if extended is False and matching_haplotype_groups.nonunique:
                    pos = pos_x + 1
                    next_pos_nonuniques = [
                        a[pos] for a in matching_haplotype_groups.nonunique[haplotype_to_extend] if a[pos] != "x"
                        ]
                    if len(set(next_pos_nonuniques)) == 1 and len(next_pos_nonuniques) > 4:
                        new_extended_haplotype = new_hap[:pos] + next_pos_nonuniques[0] + "x" * (hap_len - 1 - pos)
                        new_hap_list.append(new_extended_haplotype)
                        extended = True
                if extended is False:
                    new_hap_list.append(haplotype_to_extend)
        else:
            new_hap_list.append(haplotype_to_extend)
    return new_hap_list


def get_block_read_count(haps, var_index_in_flanking):
    """
    Return supporting reads for different haplotypes.
    Parameters:
        haps (dict of str: int): read counts of haplotypes
        var_index_in_flanking (int): index of variant site in haplotypes
    Returns:
        new_lcount (list of int): reordered read counts. first one should be WT gene
    """
    index_gene = None
    index_pseudo = None
    index_var = None
    wt_reads = 0
    lcount = []
    for i, hap in enumerate(haps):
        count = haps[hap]
        lcount.append(count)
        if hap[var_index_in_flanking] == '1':
            wt_reads += count
        if hap.count("1") == len(hap):
            index_gene = i
        elif hap.count("2") == len(hap):
            index_pseudo = i
        else:
            index_var = i
    # if no read carries WT gene base at the variant site
    if wt_reads == 0 and sum(lcount) > 10:
        return [0]
    new_lcount = None
    if index_gene is None:
        return None
    else:
        new_lcount = [lcount[index_gene]]
        if len(haps) == 2:
            new_lcount.append(lcount[1-index_gene])
        elif len(haps) == 3:
            if index_pseudo is not None:
                new_lcount.append(lcount[index_var])
                new_lcount.append(lcount[index_pseudo])
            else:
                for i in range(len(haps)):
                    if i != index_gene:
                        new_lcount.append(lcount[i])
        elif len(haps) == 4:
            for i in range(len(haps)):
                if i != index_gene:
                    new_lcount.append(lcount[i])
        else:
            return None
    return new_lcount


def get_possible_hap_cns(full_cn, haplotype_count):
    """
    Given the total copy number (3 or 4) and the number of haplotypes,
    Output all possible assignments of copy number per haplotype
    Parameters:
        full_cn (int): total cn of paralogs
        haplotype_count (int): number of haplotypes
    Returns:
        possible_hap_cns(list of tuples): each item is a tuple of
        size equal to haplotype_count, containing copy numbers of
        each haplotype
    """
    possible_hap_cns = set()
    if haplotype_count == 3:
        for i in range(1, full_cn - haplotype_count + 2):
            for j in range(1, full_cn - haplotype_count + 2):
                for k in range(1, full_cn - haplotype_count + 2):
                    if i+j+k == full_cn:
                        possible_hap_cns.add((i, j, k))
    if haplotype_count == 4:
        for i in range(1, full_cn - haplotype_count + 2):
            for j in range(1, full_cn - haplotype_count + 2):
                for k in range(1, full_cn - haplotype_count + 2):
                    for l in range(1, full_cn - haplotype_count + 2):
                        if i+j+k+l == full_cn:
                            possible_hap_cns.add((i, j, k, l))
    return possible_hap_cns


def call_hap_cn(lcount, full_cn):
    """
    Likelihood of having only one copy of WT gene, based on depth of haplotypes
    Parameters:
        lcount (list of int): read counts of haplotypes. first one should be WT gene
        full_cn (int): total cn of paralogs
    Returns:
        likelihood of the first haplotype having a cn of 1
    """
    nsum = sum(lcount)
    if full_cn == len(lcount):
        return 1
    if (
        nsum == 0 or
        full_cn < len(lcount) or
        0 in lcount or
        len(lcount) > 4 or
        full_cn >= 10
    ):
        return 0

    gene_count = lcount[0]
    # below are possible scenarios of copie numbers of WT, variant, pseudo
    if len(lcount) == 2:
        # 2 haplotypes, full_cn > 2
        prob = []
        for i in range(full_cn):
            if i > 0:
                prob_value = 1
                depthexpected = nsum * i / full_cn
                prob_value *= poisson.pmf(gene_count, depthexpected)
                depthexpected = nsum * (full_cn - i) / full_cn
                prob_value *= poisson.pmf(lcount[1], depthexpected)
                prob.append(prob_value)
        sum_prob = sum(prob)
        if sum_prob == 0:
            return 0
        post_prob = [a/sum_prob for a in prob]
        return post_prob[0]

    possible_hap_cns = get_possible_hap_cns(full_cn, len(lcount))
    prob = {}
    if len(lcount) == 3:
        # 3 haplotypes, full_cn > 3
        for possible_hap_cn in possible_hap_cns:
            prob_value = 1
            depthexpected = nsum * possible_hap_cn[0] / full_cn
            prob_value *= poisson.pmf(lcount[0], depthexpected)
            depthexpected = nsum * possible_hap_cn[1] / full_cn
            prob_value *= poisson.pmf(lcount[1], depthexpected)
            depthexpected = nsum * possible_hap_cn[2] / full_cn
            prob_value *= poisson.pmf(lcount[2], depthexpected)
            prob.setdefault(possible_hap_cn, prob_value)
    if len(lcount) == 4:
        # 4 haplotypes, full_cn > 4
        for possible_hap_cn in possible_hap_cns:
            prob_value = 1
            depthexpected = nsum * possible_hap_cn[0] / full_cn
            prob_value *= poisson.pmf(lcount[0], depthexpected)
            depthexpected = nsum * possible_hap_cn[1] / full_cn
            prob_value *= poisson.pmf(lcount[1], depthexpected)
            depthexpected = nsum * possible_hap_cn[2] / full_cn
            prob_value *= poisson.pmf(lcount[2], depthexpected)
            depthexpected = nsum * possible_hap_cn[3] / full_cn
            prob_value *= poisson.pmf(lcount[3], depthexpected)
            prob.setdefault(possible_hap_cn, prob_value)
    #print(lcount, full_cn, prob)
    sum_prob = sum(prob.values())
    if sum_prob == 0:
        return 0
    post_prob = []
    for possible_hap_cn in prob:
        if possible_hap_cn[0] == 1:
            post_prob.append(prob[possible_hap_cn]/sum_prob)
    if post_prob != []:
        return max(post_prob)
    return 0


def get_read_count_threshold(bases_at_site):
    """
    Get read count threshold depending on counts
    Parameters:
        bases_at_site (list of str): list of bases seen at a site
    Returns:
        min_count (int): minimum supporting reads to filter on
    """
    base_count = []
    for uniq_base in list(set(bases_at_site)):
        base_count.append(bases_at_site.count(uniq_base))
    min_count = 2
    if max(base_count) <= 4:
        min_count = 1
    return min_count


def join_subblocks(haplotype_per_read, lhap1, lhap2, pos1, pos2):
    """
    Merge two partially assembled haplotypes
    Parameters:
        haplotype_per_read (dict of str : str): bases of each read at differentiating sites
        lhap1 (list of str): haplotypes to be merged
        lhap2 (list of str): haplotypes to be merged
        pos1: looks for read support at one side before this position
        pos2: looks for read support at other side after this position
    Returns:
        full_hap (list of str): merged haplotypes
        dcount (dict of str : list): reads supporting each of the merged haplotypes
    """
    full_hap = []
    for sequence1 in lhap1:
        for sequence2 in lhap2:
            # Check if the pos1-pos2 region matches between the two sequences
            consensus = get_seq_consensus(sequence1, sequence2)
            if consensus is not False and consensus not in full_hap:
                    full_hap.append(consensus)
    dcount = {}
    for merged_hap in full_hap:
        for read_name in haplotype_per_read:
            read_hap = haplotype_per_read[read_name]
            match1, mismatch1, extend1 = compare_two_haps(read_hap[:pos1], merged_hap[:pos1])
            match2, mismatch2, extend2 = compare_two_haps(read_hap[pos2:], merged_hap[pos2:])
            if match1 > 0 and match2 > 0 and mismatch1 == 0 and mismatch2 == 0:
                dcount.setdefault(merged_hap, []).append(read_name)
    return full_hap, dcount


def filter_hap(haplotype_read_count, min_read=2):
    """
    Get haplotypes' supporting read count, filter by minimum read count
    Parameters:
        haplotype_read_count (dic of str : list): reads supporting each haplotype
        min_read (int): minimum supporting read count
    Returns:
        haplotype_read_count_filtered (dic of str : int): number of supporting reads
    """
    haplotype_read_count_filtered = {}
    for hap_item in haplotype_read_count:
        hap_count = sum(haplotype_read_count[hap_item])
        if hap_count >= min_read:
            haplotype_read_count_filtered.setdefault(hap_item, hap_count)
    return haplotype_read_count_filtered


def compare_two_haps(hap1, hap2):
    """
    Compare two haplotypes, and calculate number of bases for match/mismatch/extend
    Parameters:
        hap1 (str): sequence of haplotype 1
        hap2 (str): sequence of haplotype 2
    Returns:
        match (int): number of matching bases
        mismatch (int): number of mismatching bases
        extend (int): number of extending bases
    """
    assert len(hap1) == len(hap2)
    match = 0
    mismatch = 0
    extend = 0
    for i in range(len(hap1)):
        base1 = hap1[i]
        base2 = hap2[i]
        if "x" not in [base1, base2]:
            if base1 == base2:
                match += 1
            else:
                mismatch += 1
        if base1 == "x" and base2 != "x":
            extend += 1
    return match, mismatch, extend


def get_seq_consensus(hap1, hap2):
    """
    compare two haplotype blocks and produce consensus
    Parameters:
        hap1 (str): sequence of haplotype 1
        hap2 (str): sequence of haplotype 2
    Returns:
        consensus (str): consensus sequence
    """
    assert len(hap1) == len(hap2)
    if hap1.count("x") == len(hap1) or hap2.count("x") == len(hap2):
        return False
    consensus = ""
    match = 0
    for i in range(len(hap1)):
        base1 = hap1[i]
        base2 = hap2[i]
        if "x" not in [base1, base2]:
            if base1 != base2:
                return False
            match += 1
            consensus += base1
        elif base1 != "x":
            consensus += base1
        elif base2 != "x":
            consensus += base2
        else:
            consensus += "x"
    if match == 0:
        return False
    return consensus
