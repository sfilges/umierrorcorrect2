#!/usr/bin/env python3
"""Consensus sequence generation module.

Handles the generation of consensus sequences from clustered UMI reads,
including position-based consensus and MSA-based consensus.
"""

import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from itertools import groupby
from math import log10
from typing import TYPE_CHECKING, Any, Optional, Union

import pysam

from umierrorcorrect2.core.check_args import is_tool
from umierrorcorrect2.core.constants import (
    COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS,
    DEFAULT_MAPPING_QUALITY,
    MAX_PHRED_SCORE,
    PHRED_TABLE_SIZE,
    READ_GROUP_TAG,
)

if TYPE_CHECKING:
    from umierrorcorrect2.core.umi_cluster import umi_cluster

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


@dataclass
class ConsensusRead:
    """Class for representing a consensus read, useful for writing to BAM."""

    contig: str
    regionid: str
    start_pos: int
    umi_centroid: str
    umi_count: int
    name: str = field(init=False)
    count: int = field(init=False)
    _seq_parts: list[str] = field(default_factory=list, repr=False)
    _qual_parts: list[str] = field(default_factory=list, repr=False)
    _cigar_parts: list[str] = field(default_factory=list, repr=False)
    indel_read: int = 0
    nmtag: int = 0
    is_split_read: bool = False
    splits: list[Union[int, tuple[int, int]]] = field(default_factory=list)
    json: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Initialize derived fields."""
        self.count = self.umi_count
        self.name = f"Consensus_read_{self.regionid}_{self.umi_centroid}_Count={self.umi_count}"

    @property
    def seq(self) -> str:
        """Lazily join sequence parts."""
        return "".join(self._seq_parts)

    @property
    def qual(self) -> str:
        """Lazily join quality parts."""
        return "".join(self._qual_parts)

    @property
    def cigarstring(self) -> str:
        """Lazily join cigar parts."""
        return "".join(self._cigar_parts)

    def add_base(self, base: str, qual: str) -> None:
        """Add a single base and quality score to the read."""
        self._seq_parts.append(base)
        self._qual_parts.append(qual)
        self._cigar_parts.append("0")  # todo mismatch

    def add_insertion(self, sequence: str) -> None:
        """Add an insertion to the read."""
        length = len(sequence)
        self._seq_parts.append(sequence)
        self._qual_parts.append("]" * length)
        self._cigar_parts.append("1" * length)
        self.nmtag += length
        self.indel_read = 1

    def add_deletion(self, dellength: int) -> None:
        """Add a deletion to the read."""
        self._cigar_parts.append("2" * dellength)
        self.nmtag += int(dellength)
        self.indel_read = -1

    def get_cigar(self) -> tuple[tuple[int, int], ...]:
        """Generate CIGAR tuple for pysam."""
        cigarstring = self.cigarstring
        groups = groupby(cigarstring)
        cigar = tuple((int(label), sum(1 for _ in group)) for label, group in groups)
        return cigar

    def split_read(self, position1: int, position2: int) -> None:
        """Mark read as split and record split positions."""
        self.is_split_read = True
        if len(self.splits) == 0:
            self.splits.append((self.start_pos, position1))
            self.splits.append(position2)
        else:
            tmppos = self.splits[-1]
            if tmppos == position1:
                self.splits[-1] = position2
            else:
                self.splits[-1] = (tmppos, position1)  # type: ignore
                self.splits.append(position2)

    def add_json_object(self, dictionary: dict[str, Any]) -> None:
        """Add metadata dictionary."""
        self.json = dictionary

    def write_to_bam(self, f: pysam.AlignmentFile) -> None:
        """Write the consensus read to a BAM file."""
        if not self.is_split_read:
            a = pysam.AlignedSegment()
            a.query_name = self.name
            a.query_sequence = self.seq
            a.flag = 0
            a.reference_id = f.references.index(self.contig)
            a.reference_start = self.start_pos
            a.mapping_quality = DEFAULT_MAPPING_QUALITY
            a.cigar = self.get_cigar()  # type: ignore
            a.query_qualities = pysam.qualitystring_to_array(self.qual)  # type: ignore
            a.tags = (("NM", self.nmtag), ("RG", READ_GROUP_TAG))  # type: ignore
            f.write(a)
        else:
            self._write_split_read_to_bam(f)

    def _write_split_read_to_bam(self, f: pysam.AlignmentFile) -> None:
        """Helper to write split reads."""
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
                a.query_name = f"{parts[0]}_{chr(j + 97)}_Count={parts[1]}"
                j += 1
                a.flag = 0
                a.reference_id = f.references.index(self.contig)
                a.mapping_quality = DEFAULT_MAPPING_QUALITY

                # Calculate CIGAR for this segment
                endc = end + self.cigarstring[start:end].count("2")  # add 1 for each deletion
                groups = groupby(self.cigarstring[start:endc])
                cigar = tuple((int(label), sum(1 for _ in group)) for label, group in groups)

                a.cigar = cigar  # type: ignore
                a.query_qualities = pysam.qualitystring_to_array(self.qual)[start:end]  # type: ignore
                a.tags = (("NM", self.nmtag), ("RG", READ_GROUP_TAG))  # type: ignore
                f.write(a)


def get_reference_sequence(fasta: pysam.FastaFile, chrx: str, start: int, stop: int) -> str:
    """Returns the fasta sequence of the reference for a given region."""
    chrx = str(chrx)
    ref = fasta.fetch(chrx, start, stop).upper()
    return ref


def get_most_common_allele(cons_pos: dict[str, Any]) -> tuple[str, float]:
    """Calculate the allele with the highest frequency at a position.

    Args:
        cons_pos: Dictionary of base counts/qualities at a position.

    Returns:
        Tuple of (most_common_allele, percentage).
    """
    cons_dict = {}
    for x, y in cons_pos.items():
        if x in "ID":
            # Insertions/Deletions are nested dicts
            for al in y:
                cons_dict[f"{x}{al}"] = y[al]
        else:
            # Bases are lists of qualities
            cons_dict[x] = len(y)

    if not cons_dict:
        return ("N", 0.0)

    cons_allele = max(cons_dict, key=cons_dict.get)  # type: ignore
    total = sum(cons_dict.values())
    cons_percent = (cons_dict[cons_allele] / total) * 100 if total > 0 else 0.0
    return (cons_allele, cons_percent)


def get_phred(character: str) -> int:
    """Get the numeric value of ASCII character associated with phred score (offset 33)."""
    return ord(character) - 33


def get_ascii(value: int) -> str:
    """Get the ascii character for a given phred score (offset 33)."""
    return chr(value + 33)


def calc_consensus(base: str, cons_pos: dict[str, list[int]]) -> float:
    """Calculate the combined probability score for a base at a position."""
    prod = 1.0
    for nucl in cons_pos:
        # If the read base matches the candidate consensus base
        if nucl in base:
            for phred in cons_pos[nucl]:
                if phred != 0:
                    prod *= _phred_to_prob(phred)
        # If the read base is standard ATCG but does not match candidate
        elif nucl in "ATCG":
            for phred in cons_pos[nucl]:
                if phred != 0:
                    prod *= _phred_to_error(phred)
    return prod


def calc_consensus_probabilities(cons_pos: dict[str, list[int]]) -> tuple[str, int]:
    """Calculate the probability of consensus at a given position.

    Returns:
        Tuple of (consensus_base, consensus_phred_score).
    """
    p = {base: calc_consensus(base, cons_pos) for base in "ATCG"}
    denom = sum(p.values())

    if denom > 0:
        probs = {base: p[base] / denom for base in "ATCG"}
    else:
        probs = dict.fromkeys("ATCG", 0.0)

    cons_base = max(probs, key=probs.get)  # type: ignore
    prob_val = probs[cons_base]

    if prob_val >= 1.0:  # Floating point tolerance handling implied
        cons_phred = MAX_PHRED_SCORE
    else:
        try:
            cons_phred = round(-10 * log10(1 - prob_val))
        except (ValueError, ZeroDivisionError):
            cons_phred = MAX_PHRED_SCORE

        if cons_phred > MAX_PHRED_SCORE:
            cons_phred = MAX_PHRED_SCORE

    return (cons_base, cons_phred)


def get_position_coverage(covpos: dict[str, Any]) -> int:
    """Calculate total coverage at a position, including deletions."""
    coverage = sum([len(covpos[x]) for x in covpos if x not in ["D", "I"]])
    if "D" in covpos:
        for numseqs in covpos["D"].values():
            coverage += numseqs
    return coverage


def _add_base_to_consensus(consensus: dict[int, dict], refpos: int, base: str, phred: int) -> None:
    """Add a base with quality score to the consensus dictionary."""
    if refpos not in consensus:
        consensus[refpos] = {}
    if base not in consensus[refpos]:
        consensus[refpos][base] = []
    consensus[refpos][base].append(phred)


def _process_read_without_indel(read: pysam.AlignedSegment, consensus: dict[int, dict]) -> None:
    """Process a read without insertions or deletions."""
    sequence = read.query_sequence
    qual = read.query_qualities
    if sequence is None or qual is None:
        return

    for qpos, refpos in read.get_aligned_pairs(matches_only=True):
        base = sequence[qpos]
        # Pysam qualities are already integers, no need for ord() - 33 if using query_qualities
        # But legacy code used query_alignment_qualities or string...
        # Wait, previous code used read.qual which was likely a string?
        # Checking old code: read.qual is property of consensus_read, but here 'read' is pysam AlignedSegment?
        # Yes, from _build_consensus_dict -> group_seqs -> pysam entries.
        # pysam.AlignedSegment.qual is NOT a property. query_qualities is array('B'). query_qualities is string (deprecated?) or query_quality_string.
        # Let's stick to using what works for pysam.
        # Original code: get_phred(qual[qpos]).
        # If qual is array, no need for get_phred.
        # Let's check how 'read' is passed. It comes from f.fetch().
        _add_base_to_consensus(consensus, refpos, base, qual[qpos])


def _process_read_with_indel(read: pysam.AlignedSegment, consensus: dict[int, dict]) -> None:
    """Process a read that contains insertions or deletions."""
    positions = read.get_aligned_pairs(matches_only=True)
    if not positions:
        return

    q_prev, ref_prev = positions[0]
    q_prev -= 1
    ref_prev -= 1

    sequence = read.query_sequence
    qual = read.query_qualities
    if sequence is None or qual is None:
        return

    for qpos, refpos in positions:
        base = sequence[qpos]
        q_score = qual[qpos]

        if qpos != q_prev + 1:
            # Insertion
            allele = sequence[q_prev + 1 : qpos]
            if refpos not in consensus:
                consensus[refpos] = {}
            if "I" not in consensus[refpos]:
                consensus[refpos]["I"] = {}
            if allele not in consensus[refpos]["I"]:
                consensus[refpos]["I"][allele] = 0
            consensus[refpos]["I"][allele] += 1
            _add_base_to_consensus(consensus, refpos, base, q_score)

        elif refpos != ref_prev + 1:
            # Deletion
            dellength = refpos - (ref_prev + 1)
            delpos = refpos - dellength
            if delpos not in consensus:
                consensus[delpos] = {}
            if "D" not in consensus[delpos]:
                consensus[delpos]["D"] = {}
            if dellength not in consensus[delpos]["D"]:
                consensus[delpos]["D"][dellength] = 0
            consensus[delpos]["D"][dellength] += 1
            _add_base_to_consensus(consensus, refpos, base, q_score)

        else:
            _add_base_to_consensus(consensus, refpos, base, q_score)

        q_prev = qpos
        ref_prev = refpos


def _build_consensus_dict(group_seqs: list[pysam.AlignedSegment]) -> dict[int, dict]:
    """Build consensus dictionary from a group of reads."""
    consensus: dict[int, dict] = {}
    for read in group_seqs:
        if read.cigarstring:
            if "I" not in read.cigarstring and "D" not in read.cigarstring:
                _process_read_without_indel(read, consensus)
            else:
                _process_read_with_indel(read, consensus)
    return consensus


def _get_consensus_base_simple(consensus_pos: dict, poscov: int) -> tuple[str, int]:
    """Calculate consensus base (simple majority or prob based on coverage)."""
    if poscov < COVERAGE_THRESHOLD_FOR_SIMPLE_CONSENSUS:
        return calc_consensus_probabilities(consensus_pos)
    else:
        cons_base, _ = get_most_common_allele(consensus_pos)
        return cons_base, MAX_PHRED_SCORE


def get_consensus_position_based(
    group_seqs: list[pysam.AlignedSegment],
    contig: str,
    regionid: str,
    indel_freq_threshold: float,
    umi_info: "umi_cluster",
    consensus_freq_threshold: float,
    output_json: bool,
) -> Optional[ConsensusRead]:
    """Generate a consensus sequence from a group of reads (position-based)."""
    consensus = _build_consensus_dict(group_seqs)
    if not consensus:
        return None

    consensus_sorted = sorted(consensus)
    consread = None
    add_consensus = True
    deletion_skip_positions: list[int] = []
    prevpos = consensus_sorted[0] - 1

    for pos in consensus_sorted:
        if pos in deletion_skip_positions:
            prevpos = pos
            continue

        poscov = get_position_coverage(consensus[pos])

        # Initialize consensus read if first position
        if not consread:
            consread = ConsensusRead(contig, regionid, pos, umi_info.centroid, umi_info.count)

        # Handle gaps (Ns or Splits)
        if pos != prevpos + 1 and (prevpos + 1) not in deletion_skip_positions:
            for i in range(prevpos + 1, pos):
                if i not in deletion_skip_positions:
                    consread.add_base("N", get_ascii(0))
            consread.split_read(prevpos + 1, pos)

        # Logic for Insertion
        if "I" in consensus[pos] and poscov >= 2:
            cons_dict = consensus[pos]["I"]
            cons_allele = max(cons_dict, key=cons_dict.get)  # type: ignore
            cons_num = cons_dict[cons_allele]
            percent = (cons_num / poscov) * 100.0

            if percent >= indel_freq_threshold:
                consread.add_insertion(cons_allele)

            # Remove insertion info to proceed with base calculation
            # Note: This follows original logic where base at insertion site is also processed
            del consensus[pos]["I"]

            cons_base, cons_qual = _get_consensus_base_simple(consensus[pos], poscov)

            if not cons_base.startswith("D"):
                consread.add_base(cons_base, get_ascii(cons_qual))

        # Logic for Deletion
        elif "D" in consensus[pos] and poscov >= 2:
            allele, percent = get_most_common_allele(consensus[pos])

            if allele.startswith("D"):
                if percent >= indel_freq_threshold:
                    dellength = int(allele.lstrip("D"))
                    consread.add_deletion(dellength)
                    if dellength > 1:
                        # Mark positions to skip
                        for i in range(1, dellength):
                            deletion_skip_positions.append(pos + i)
                else:
                    consread.add_base("N", get_ascii(0))
                    add_consensus = False
            elif percent >= indel_freq_threshold:
                # Majority is NOT deletion
                cons_base, cons_qual = _get_consensus_base_simple(consensus[pos], poscov)
                consread.add_base(cons_base, get_ascii(cons_qual))
            else:
                # Ambiguous
                consread.add_base("N", get_ascii(0))
                add_consensus = False

        # Logic for Match/Mismatch
        elif poscov >= 2:
            cons_base, cons_qual = _get_consensus_base_simple(consensus[pos], poscov)

            if consensus_freq_threshold:
                # Check frequency threshold
                if len(consensus[pos]) == 1:
                    consread.add_base(cons_base, get_ascii(cons_qual))
                elif cons_base in consensus[pos]:
                    count_base = len(consensus[pos][cons_base])
                    percent = (count_base / len(group_seqs)) * 100.0
                    if percent >= consensus_freq_threshold:
                        consread.add_base(cons_base, get_ascii(cons_qual))
                    else:
                        consread.add_base("X", get_ascii(0))
                        add_consensus = False
                else:
                    consread.add_base("X", get_ascii(0))
                    add_consensus = False
            else:
                consread.add_base(cons_base, get_ascii(cons_qual))

        else:
            # Low coverage (< 2)
            consread.add_base("N", get_ascii(0))
            if consread.is_split_read and consread.splits[-1] == pos - 1:
                # extend gap
                if isinstance(consread.splits[-1], int):
                    consread.splits[-1] = pos
                # Logic for extending split read in original code was:
                # if consread.splits[-1] == pos - 1: consread.splits[-1] = pos
            else:
                consread.split_read(pos, pos + 1)

        prevpos = pos

    if add_consensus and consread:
        if output_json:
            counts = Counter([x.query_sequence for x in group_seqs if x.query_sequence])
            consread.add_json_object(dict(counts))  # type: ignore
        return consread

    return None


def get_consensus_most_common(
    group_seqs: list[pysam.AlignedSegment],
    contig: str,
    regionid: str,
    indel_freq_threshold: float,
    umi_info: "umi_cluster",
    consensus_freq_threshold: float,
    output_json: bool,
) -> Optional[ConsensusRead]:
    """Generate consensus by taking the most common full sequence."""
    total_seqs = len(group_seqs)
    # Filter None sequences
    seqs = [x.query_sequence for x in group_seqs if x.query_sequence]
    if not seqs:
        return None

    counts = Counter(seqs)
    most_common_seq, n = counts.most_common(1)[0]
    percentage = (n / total_seqs) * 100

    if percentage >= consensus_freq_threshold:
        pos = min([x.reference_start for x in group_seqs if x.reference_start is not None])
        consread = ConsensusRead(contig, regionid, pos, umi_info.centroid, umi_info.count)

        for base in most_common_seq:
            consread.add_base(base, get_ascii(MAX_PHRED_SCORE))

        if output_json:
            consread.add_json_object(dict(counts))
        return consread

    return None


def get_consensus_msa(
    group_seqs: list[pysam.AlignedSegment],
    contig: str,
    regionid: str,
    indel_frequency_threshold: float,
    umi_info: "umi_cluster",
    consensus_freq_threshold: float,
    output_json: bool,
) -> Optional[ConsensusRead]:
    """Generate consensus using Multiple Sequence Alignment (MAFFT)."""

    # Dependency check
    if not is_tool("mafft"):
        raise RuntimeError("MAFFT is not installed or not in PATH. Please install MAFFT to use MSA consensus.")

    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        for a in group_seqs:
            name = a.query_name  # type: ignore
            seq = a.query_sequence
            if seq:
                f.write(f">{name}\n{seq}\n")
        f.flush()

        try:
            # Use text mode for subprocess to avoid bytes handling issues
            output = subprocess.check_output(["mafft", "--quiet", f.name], text=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            # Log warning or re-raise? For now, re-raise with context
            raise RuntimeError(f"MAFFT execution failed: {e.stderr}") from e

    sequences: list[str] = []
    current_seq_parts: list[str] = []

    for line in output.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_seq_parts:
                sequences.append("".join(current_seq_parts))
                current_seq_parts = []
        else:
            current_seq_parts.append(line)
    if current_seq_parts:
        sequences.append("".join(current_seq_parts))

    if not sequences:
        return None

    consensus: dict[int, Counter] = {}
    seq_len = len(sequences[0])

    for seq in sequences:
        # Check if MSA produced equal length sequences (it should)
        if len(seq) != seq_len:
            continue
        for i, base in enumerate(seq):
            if i not in consensus:
                consensus[i] = Counter()
            consensus[i][base] += 1

    pos = min([x.reference_start for x in group_seqs if x.reference_start is not None])
    consread = ConsensusRead(contig, regionid, pos, umi_info.centroid, umi_info.count)
    add_consensus = True
    n = len(sequences)

    for i in sorted(consensus):
        b = consensus[i].most_common(1)[0]
        fraction = (b[1] / n) * 100
        base_char = b[0]

        if base_char != "-":
            if fraction >= consensus_freq_threshold:
                consread.add_base(base_char, get_ascii(MAX_PHRED_SCORE))
            else:
                consread.add_base("N", get_ascii(0))
                add_consensus = False

    if add_consensus:
        if output_json:
            counts = Counter([x.query_sequence for x in group_seqs if x.query_sequence])
            consread.add_json_object(dict(counts))
        return consread

    return None


def get_all_consensus(
    position_matrix: dict[str, list[pysam.AlignedSegment]],
    umis: dict[str, "umi_cluster"],
    contig: str,
    regionid: str,
    indel_frequency_cutoff: float,
    consensus_frequency_cutoff: float,
    output_json: bool,
) -> dict[str, Optional[ConsensusRead]]:
    """Get the consensus sequences for all umis (position-based)."""
    consensuses = {}
    for umi, reads in position_matrix.items():
        consensuses[umi] = get_consensus_position_based(
            reads,
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensuses


def get_all_consensus_most_common(
    position_matrix: dict[str, list[pysam.AlignedSegment]],
    umis: dict[str, "umi_cluster"],
    contig: str,
    regionid: str,
    indel_frequency_cutoff: float,
    consensus_frequency_cutoff: float,
    output_json: bool,
) -> dict[str, Optional[ConsensusRead]]:
    """Get the consensus sequences for all umis (most common sequence)."""
    consensus_seq = {}
    for umi, reads in position_matrix.items():
        consensus_seq[umi] = get_consensus_most_common(
            reads,
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensus_seq


def get_all_consensus_msa(
    position_matrix: dict[str, list[pysam.AlignedSegment]],
    umis: dict[str, "umi_cluster"],
    contig: str,
    regionid: str,
    indel_frequency_cutoff: float,
    consensus_frequency_cutoff: float,
    output_json: bool,
) -> dict[str, Optional[ConsensusRead]]:
    """Get the consensus sequences for all umis (MSA based)."""
    consensus_seq = {}
    for umi, reads in position_matrix.items():
        consensus_seq[umi] = get_consensus_msa(
            reads,
            contig,
            regionid,
            indel_frequency_cutoff,
            umis[umi],
            consensus_frequency_cutoff,
            output_json,
        )
    return consensus_seq


def get_cons_dict(
    bamfilename: str, umis: dict[str, "umi_cluster"], contig: str, start: int, end: int, include_singletons: bool
) -> tuple[dict[str, list[pysam.AlignedSegment]], dict[str, pysam.AlignedSegment]]:
    """Read BAM file and group reads by UMI cluster."""
    position_matrix: dict[str, list[pysam.AlignedSegment]] = {}
    singleton_matrix: dict[str, pysam.AlignedSegment] = {}

    with pysam.AlignmentFile(bamfilename, "rb") as f:
        alignment = f.fetch(contig, start, end)
        for read in alignment:
            # Use rsplit with maxsplit=1 - more efficient than split(":")[-1]
            barcode = read.qname.rsplit(":", 1)[-1]  # type: ignore
            pos = read.reference_start

            if pos is not None and start <= pos <= end:
                umi_info = umis.get(barcode)
                if umi_info is not None:
                    cluster = umi_info.centroid
                    cluster_size = umi_info.count

                    if cluster_size > 1:
                        if cluster not in position_matrix:
                            position_matrix[cluster] = []
                        position_matrix[cluster].append(read)
                    elif include_singletons:
                        # For singletons, just keep one read (logic from original code)
                        if cluster not in singleton_matrix:
                            singleton_matrix[cluster] = read

    return (position_matrix, singleton_matrix)


def write_singleton_reads(
    singleton_matrix: dict[str, pysam.AlignedSegment], region_id: str, g: pysam.AlignmentFile
) -> None:
    """Write singleton reads directly to the output BAM."""
    for umi, read in singleton_matrix.items():
        read.query_name = f"Singleton_read_{region_id}_{umi}_Count=1"  # type: ignore
        g.write(read)
