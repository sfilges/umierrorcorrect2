#!/usr/bin/env python3

import subprocess
import tempfile
from collections import Counter
from itertools import groupby
from math import log10

import pysam

from umierrorcorrect.core.constants import (
    COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS,
    DEFAULT_MAPPING_QUALITY,
    MAX_PHRED_SCORE,
    PHRED_TABLE_SIZE,
    READ_GROUP_TAG,
)

# Pre-computed lookup tables for phred score conversions (avoid repeated 10**x calculations)
# Phred scores typically range from 0-93 (ASCII 33-126)
_PHRED_TO_PROB = tuple(1.0 - (10.0 ** (-q / 10.0)) for q in range(PHRED_TABLE_SIZE))
_PHRED_TO_ERROR = tuple(10.0 ** (-q / 10.0) for q in range(PHRED_TABLE_SIZE))


def _phred_to_prob(phred: int) -> float:
    """Convert phred score to probability of correctness using lookup table."""
    if phred < PHRED_TABLE_SIZE:
        return _PHRED_TO_PROB[phred]
    return _PHRED_TO_PROB[PHRED_TABLE_SIZE - 1]  # Cap at max


def _phred_to_error(phred: int) -> float:
    """Convert phred score to error probability using lookup table."""
    if phred < PHRED_TABLE_SIZE:
        return _PHRED_TO_ERROR[phred]
    return _PHRED_TO_ERROR[PHRED_TABLE_SIZE - 1]  # Cap at max


class consensus_read:
    """Class for representing a consensus read, useful for writing to BAM"""

    def __init__(self, contig, regionid, position_start, name, count):
        self.contig = contig
        self.start_pos = position_start
        self._seq_parts = []  # Use list for O(1) append
        self._qual_parts = []  # Use list for O(1) append
        self._cigar_parts = []  # Use list for O(1) append
        self.indel_read = 0
        self.nmtag = 0
        self.is_split_read = False
        self.splits = []
        self.json = {}
        self.name = f"Consensus_read_{regionid}_{name}_Count={count}"
        self.count = count

    @property
    def seq(self):
        """Lazily join sequence parts."""
        return "".join(self._seq_parts)

    @property
    def qual(self):
        """Lazily join quality parts."""
        return "".join(self._qual_parts)

    @property
    def cigarstring(self):
        """Lazily join cigar parts."""
        return "".join(self._cigar_parts)

    def add_base(self, base, qual):
        self._seq_parts.append(base)
        self._qual_parts.append(qual)
        self._cigar_parts.append("0")  # todo mismatch

    def add_insertion(self, sequence):
        self._seq_parts.append(sequence)
        self._qual_parts.append("]" * len(sequence))
        self._cigar_parts.append("1" * len(sequence))
        self.nmtag += len(sequence)
        self.indel_read = 1

    def add_deletion(self, dellength):
        self._cigar_parts.append("2" * dellength)
        self.nmtag += int(dellength)
        self.indel_read = -1

    def get_cigar(self):
        cigarstring = self.cigarstring
        groups = groupby(cigarstring)
        cigar = tuple((int(label), sum(1 for _ in group)) for label, group in groups)
        return cigar

    def split_read(self, position1, position2):
        self.is_split_read = True
        if len(self.splits) == 0:
            self.splits.append((self.start_pos, position1))
            self.splits.append(position2)
        else:
            tmppos = self.splits[-1]
            if tmppos == position1:
                self.splits[-1] = position2
            else:
                self.splits[-1] = (tmppos, position1)
                self.splits.append(position2)

    def add_json_object(self, dictionary):
        self.json = dictionary

    def write_to_bam(self, f):
        if self.is_split_read is False:
            a = pysam.AlignedSegment()
            a.query_name = self.name
            a.query_sequence = self.seq
            a.flag = 0
            a.reference_id = f.references.index(self.contig)
            a.reference_start = self.start_pos
            a.mapping_quality = DEFAULT_MAPPING_QUALITY
            a.cigar = self.get_cigar()
            a.query_qualities = pysam.qualitystring_to_array(self.qual)
            a.tags = (("NM", self.nmtag), ("RG", READ_GROUP_TAG))
            f.write(a)
        else:
            j = 0
            for i, s in enumerate(self.splits):
                a = pysam.AlignedSegment()
                if isinstance(s, tuple):
                    start = s[0] - self.start_pos
                    end = s[1] - self.start_pos
                    a.reference_start = s[0]
                else:
                    start = s - self.start_pos
                    end = len(self.seq)
                    a.reference_start = s
                a.query_sequence = self.seq[start:end]
                if a.query_sequence:
                    parts = self.name.split("_Count=")
                    a.query_name = parts[0] + "_" + chr(j + 97) + "_Count=" + parts[1]
                    j += 1
                    a.flag = 0
                    a.reference_id = f.references.index(self.contig)
                    a.mapping_quality = DEFAULT_MAPPING_QUALITY
                    endc = end + self.cigarstring[start:end].count("2")  # add 1 for each deletion
                    groups = groupby(self.cigarstring[start:endc])
                    cigar = tuple((int(label), sum(1 for _ in group)) for label, group in groups)
                    a.cigar = cigar
                    a.query_qualities = pysam.qualitystring_to_array(self.qual)[start:end]
                    a.tags = (("NM", self.nmtag), ("RG", READ_GROUP_TAG))
                    f.write(a)
            # s=self.splits[-1]
            # a = pysam.AlignedSegment()
            # a.query_name = self.name + '_b'
            # start = self.splits[-1] - self.start_pos
            # a.query_sequence = self.seq[start:]
            # a.flag = 0
            # a.reference_id = f.references.index(self.contig)
            # a.reference_start = s
            # a.mapping_quality = 60
            # groups = groupby(self.cigarstring[start:])
            # cigar = tuple((int(label),
            #               sum(1 for _ in group)) for label, group in groups)
            # a.cigar = cigar
            # a.query_qualities = pysam.qualitystring_to_array(self.qual)[start:]
            # a.tags = (("NM", self.nmtag), ("RG", "L1"))
            # f.write(a)


def get_reference_sequence(fasta, chrx, start, stop):
    """Returns the fasta sequence of the reference for a given region"""
    chrx = str(chrx)
    ref = fasta.fetch(chrx, start, stop).upper()
    return ref


def get_most_common_allele(cons_pos):
    """Calculate the allele frequencies at one position and returns the
    allele ith the highest frequency."""
    cons_dict = {}
    for x, y in cons_pos.items():
        if x in "ID":
            for al in y:
                cons_dict[x + f"{al}"] = y[al]
        else:
            cons_dict[x] = len(y)
    cons_allele = max(cons_dict, key=cons_dict.get)  # get highest count allele
    cons_percent = (cons_dict[cons_allele] / sum(cons_dict.values())) * 100
    return (cons_allele, cons_percent)


def get_phred(character):
    """Get the numeric value of ASCII character associated with phred score (offset 33)"""
    value = ord(character) - 33
    return value


def get_ascii(value):
    """Get the ascii character for a given phred score (offset 33)"""
    ascii_letter = chr(value + 33)
    return ascii_letter


def calc_consensus(base, cons_pos):
    """Function for calculating the combined score for a base at a position.

    Uses pre-computed lookup tables for phred score conversions to avoid
    expensive 10**x calculations in the hot loop.
    """
    prod = 1.0
    for nucl in cons_pos:
        if nucl in base:
            for phred in cons_pos[nucl]:
                if phred != 0:
                    prod *= _phred_to_prob(phred)
        elif nucl in "ATCG":
            for phred in cons_pos[nucl]:
                if phred != 0:
                    prod *= _phred_to_error(phred)
    return prod


def calc_consensus_probabilities(cons_pos):
    """Function for calculating the probability of consensus at a given position
    Return the base with the highest probability"""
    p = {base: calc_consensus(base, cons_pos) for base in "ATCG"}
    denom = sum(p.values())
    if denom > 0:
        probs = {base: p[base] / denom for base in "ATCG"}
    else:
        probs = dict.fromkeys("ATCG", 0)
    cons_base = max(probs, key=probs.get)
    if probs[cons_base] == 1:
        cons_phred = MAX_PHRED_SCORE
    else:
        cons_phred = round(-10 * log10(1 - probs[cons_base]))
        if cons_phred > MAX_PHRED_SCORE:
            cons_phred = MAX_PHRED_SCORE
    return (cons_base, cons_phred)


def get_position_coverage(covpos):
    coverage = 0
    coverage = sum([len(covpos[x]) for x in covpos if x not in ["D", "I"]])
    if "D" in covpos:
        for numseqs in covpos["D"].values():
            coverage += numseqs
    return coverage


def _add_base_to_consensus(consensus, refpos, base, phred):
    """Add a base with quality score to the consensus dictionary."""
    if refpos not in consensus:
        consensus[refpos] = {}
    if base not in consensus[refpos]:
        consensus[refpos][base] = []
    consensus[refpos][base].append(phred)


def _process_read_without_indel(read, consensus):
    """Process a read without insertions or deletions."""
    sequence = read.seq
    qual = read.qual
    for qpos, refpos in read.get_aligned_pairs(matches_only=True):
        base = sequence[qpos]
        _add_base_to_consensus(consensus, refpos, base, get_phred(qual[qpos]))


def _process_read_with_indel(read, consensus):
    """Process a read that contains insertions or deletions."""
    positions = read.get_aligned_pairs(matches_only=True)
    q, ref = positions[0]
    q = q - 1
    ref = ref - 1

    for qpos, refpos in positions:
        sequence = read.seq
        qual = read.qual
        base = sequence[qpos]

        if qpos != q + 1:
            # insertion
            allele = read.seq[q + 1 : qpos]
            if refpos not in consensus:
                consensus[refpos] = {}
            if "I" not in consensus[refpos]:
                consensus[refpos]["I"] = {}
            if allele not in consensus[refpos]["I"]:
                consensus[refpos]["I"][allele] = 0
            consensus[refpos]["I"][allele] += 1
            _add_base_to_consensus(consensus, refpos, base, get_phred(qual[qpos]))

        elif refpos != ref + 1:
            # deletion
            dellength = refpos - (ref + 1)
            delpos = refpos - dellength
            if delpos not in consensus:
                consensus[delpos] = {}
            if "D" not in consensus[delpos]:
                consensus[delpos]["D"] = {}
            if dellength not in consensus[delpos]["D"]:
                consensus[delpos]["D"][dellength] = 0
            consensus[delpos]["D"][dellength] += 1
            _add_base_to_consensus(consensus, refpos, base, get_phred(qual[qpos]))

        else:
            _add_base_to_consensus(consensus, refpos, base, get_phred(qual[qpos]))

        q = qpos
        ref = refpos


def _get_consensus_base_and_qual(consensus_pos, poscov):
    """Calculate the consensus base and quality for a position."""
    if poscov < COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS:
        return calc_consensus_probabilities(consensus_pos)
    else:
        cons_base, _ = get_most_common_allele(consensus_pos)
        return cons_base, MAX_PHRED_SCORE


def _build_consensus_dict(group_seqs):
    """Build consensus dictionary from a group of reads."""
    consensus = {}
    for read in group_seqs:
        if read.cigarstring:
            if "I" not in read.cigarstring and "D" not in read.cigarstring:
                _process_read_without_indel(read, consensus)
            else:
                _process_read_with_indel(read, consensus)
    return consensus


def getConsensus3(group_seqs, contig, regionid, indel_freq_threshold, umi_info, consensus_freq_threshold, output_json):
    """Takes a list of pysam entries (rows in the BAM file) as input and generates a consensus sequence."""
    consensus = _build_consensus_dict(group_seqs)
    if len(consensus) > 0:
        # generate the consensus sequence
        consensus_sorted = sorted(consensus)
        consread = None
        add_consensus = True
        skippos = []  # if position is del
        prevpos = consensus_sorted[0] - 1
        for pos in sorted(consensus_sorted):
            if pos not in skippos:
                poscov = get_position_coverage(consensus[pos])
                if not consread:
                    consread = consensus_read(contig, regionid, pos, umi_info.centroid, umi_info.count)
                if not pos == prevpos + 1 and prevpos + 1 not in skippos:
                    for i in range(prevpos + 1, pos):
                        if i not in skippos:
                            consread.add_base("N", get_ascii(0))
                    consread.split_read(prevpos + 1, pos)
                if "I" in consensus[pos] and poscov >= 2:
                    # first add the insertion if it is in the majority of the reads, then add the base at the next position
                    cons_dict = consensus[pos]["I"]
                    cons_allele = max(cons_dict, key=cons_dict.get)
                    cons_num = cons_dict[cons_allele]
                    percent = (cons_num / poscov) * 100.0
                    if percent >= indel_freq_threshold:
                        sequence = cons_allele
                        consread.add_insertion(sequence)
                    del consensus[pos]["I"]
                    if poscov < COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS:
                        cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                    else:
                        cons_base, percent = get_most_common_allele(consensus[pos])
                        cons_qual = MAX_PHRED_SCORE
                    if not cons_base.startswith("D"):
                        consread.add_base(cons_base, get_ascii(cons_qual))

                elif "D" in consensus[pos] and poscov >= 2:
                    # add the deletions
                    a, percent = get_most_common_allele(consensus[pos])
                    if a.startswith("D"):
                        if percent >= indel_freq_threshold:
                            dellength = int(a.lstrip("D"))
                            # print(dellength)
                            consread.add_deletion(dellength)
                            if dellength > 1:
                                for i in range(1, dellength):
                                    skippos.append(pos + i)
                        else:
                            consread.add_base("N", get_ascii(0))
                            add_consensus = False
                    elif percent >= indel_freq_threshold:
                        if poscov < COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS:
                            cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                        else:
                            cons_base, percent = get_most_common_allele(consensus[pos])
                            cons_qual = MAX_PHRED_SCORE
                        consread.add_base(cons_base, get_ascii(cons_qual))
                    else:
                        consread.add_base("N", get_ascii(0))
                        add_consensus = False
                elif poscov >= 2:
                    # no indel
                    if poscov < COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS:
                        cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                    else:
                        cons_base, percent = get_most_common_allele(consensus[pos])
                        cons_qual = MAX_PHRED_SCORE
                    if consensus_freq_threshold:  # test if not None
                        if len(consensus[pos]) == 1:  # 100%
                            consread.add_base(cons_base, get_ascii(cons_qual))
                        else:
                            if cons_base in consensus[pos]:
                                percent = (len(consensus[pos][cons_base]) / len(group_seqs)) * 100.0
                                if percent >= consensus_freq_threshold:  # consensus frequency above threshold
                                    consread.add_base(cons_base, get_ascii(cons_qual))
                                else:
                                    consread.add_base("X", get_ascii(0))
                                    add_consensus = False
                    else:
                        consread.add_base(cons_base, get_ascii(cons_qual))
                else:
                    consread.add_base("N", get_ascii(0))
                    if consread.is_split_read and consread.splits[-1] == pos - 1:
                        consread.splits[-1] = pos  # extend gap
                    else:
                        consread.split_read(pos, pos + 1)
                # if umi_info.centroid == 'TCCTCACG':
                #        print(consread.start_pos,consread.splits)
                #        print(consread.seq)
                #        print(consread.qual)
                #        print(consread.cigarstring)
            prevpos = pos

        if add_consensus:
            if output_json:
                counts = Counter([x.seq for x in group_seqs])
                consread.add_json_object(counts)
            return consread
        else:
            return None
    else:
        return None


def getConsensusMostCommon(
    group_seqs, contig, regionid, indel_freq_threshold, umi_info, consensus_freq_threshold, output_json
):
    total_seqs = len(group_seqs)
    counts = Counter([x.seq for x in group_seqs])
    most_common_seq, n = counts.most_common(1)[0]
    percentage = (n / total_seqs) * 100
    if percentage >= consensus_freq_threshold:
        pos = min([x.pos for x in group_seqs])
        consread = consensus_read(contig, regionid, pos, umi_info.centroid, umi_info.count)
        if most_common_seq:
            for base in most_common_seq:
                consread.add_base(base, get_ascii(MAX_PHRED_SCORE))
            if output_json:
                consread.add_json_object(counts)
        return consread
    else:
        return None


def getConsensusMSA(
    group_seqs, contig, regionid, indel_frequency_threshold, umi_info, consensus_freq_threshold, output_json
):
    with tempfile.NamedTemporaryFile() as f:
        for a in group_seqs:
            name = a.query_name
            seq = a.seq
            if seq:
                f.write(b">" + bytes(name, "utf8") + b"\n" + bytes(seq, "utf8") + b"\n")
        f.seek(0)
        output = subprocess.check_output(["mafft", "--quiet", f.name])
        sequences = []
        printseq = False
        s = ""
        for line in output.decode().split("\n"):
            if line.startswith(">"):
                if printseq:
                    sequences.append(s)
                printseq = True
                s = ""
            else:
                line = line.rstrip()
                s += line
        sequences.append(s)
        consensus = {}
        for seq in sequences:
            for i in range(len(sequences[0])):
                base = seq[i]
                if i not in consensus:
                    consensus[i] = Counter()
                consensus[i][base] += 1
        pos = min([x.pos for x in group_seqs])
        consread = consensus_read(contig, regionid, pos, umi_info.centroid, umi_info.count)
        add_consensus = True
        n = len(sequences)
        for i in sorted(consensus):
            b = consensus[i].most_common(1)[0]
            fraction = (b[1] / n) * 100
            if b[0] not in "-":
                if fraction >= consensus_freq_threshold:
                    consread.add_base(b[0], get_ascii(MAX_PHRED_SCORE))
                else:
                    consread.add_base("N", get_ascii(0))
                    add_consensus = False
        if add_consensus:
            if output_json:
                counts = Counter([x.seq for x in group_seqs])
                consread.add_json_object(counts)
            return consread
        else:
            return None


def get_all_consensus(
    position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
):
    """Get the consensus sequences for all umis"""
    consensuses = {}
    for umi in position_matrix:
        consensuses[umi] = getConsensus3(
            position_matrix[umi],
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensuses


def get_all_consensus_most_common(
    position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
):
    consensus_seq = {}
    for umi in position_matrix:
        consensus_seq[umi] = getConsensusMostCommon(
            position_matrix[umi],
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensus_seq


def get_all_consensus_msa(
    position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
):
    consensus_seq = {}
    for umi in position_matrix:
        consensus_seq[umi] = getConsensusMSA(
            position_matrix[umi],
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensus_seq


def get_cons_dict(bamfilename, umis, contig, start, end, include_singletons):
    position_matrix = {}
    singleton_matrix = {}
    with pysam.AlignmentFile(bamfilename, "rb") as f:
        alignment = f.fetch(contig, start, end)
        for read in alignment:
            # Use rsplit with maxsplit=1 - more efficient than split(":")[-1]
            # as it only splits once from the right
            barcode = read.qname.rsplit(":", 1)[-1]
            pos = read.pos
            if start <= pos <= end:
                umi_info = umis.get(barcode)
                if umi_info is not None:
                    cluster = umi_info.centroid
                    cluster_size = umi_info.count
                    if cluster_size > 1:
                        if cluster not in position_matrix:
                            position_matrix[cluster] = []
                        position_matrix[cluster].append(read)
                    elif include_singletons:
                        if cluster not in singleton_matrix:
                            singleton_matrix[cluster] = read
    return (position_matrix, singleton_matrix)


def write_singleton_reads(singleton_matrix, region_id, g):
    for umi, read in singleton_matrix.items():
        read.query_name = f"Singleton_read_{region_id}_{umi}_Count=1"
        g.write(read)
