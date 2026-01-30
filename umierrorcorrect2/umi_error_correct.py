#!/usr/bin/env python
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any

import pysam

from umierrorcorrect2.core.consensus import (
    get_all_consensus,
    get_all_consensus_most_common,
    get_all_consensus_msa,
    get_cons_dict,
    get_reference_sequence,
    write_singleton_reads,
)
from umierrorcorrect2.core.constants import DEFAULT_FAMILY_SIZES
from umierrorcorrect2.core.get_cons_info import calc_major_nonref_allele_frequency, get_cons_info, write_consensus
from umierrorcorrect2.core.get_regions_from_bed import merge_regions, read_bed, sort_regions
from umierrorcorrect2.core.group import read_bam_from_bed, read_bam_from_tag, readBam
from umierrorcorrect2.core.logging_config import get_logger
from umierrorcorrect2.core.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from umierrorcorrect2.models.models import UMIErrorCorrectConfig

logger = get_logger(__name__)

# Column indices for consensus file parsing
CONS_FILE_POS_COL = 2  # Column index for position
CONS_FILE_FSIZE_COL = 13  # Column index for family size
CONS_FILE_ALLELE_START = 5  # Start of allele count columns
CONS_FILE_ALLELE_END = 15  # End of allele count columns (exclusive)


def write_to_json(cons_read: Any) -> dict[str, Any]:
    """Convert a consensus read to a JSON-serializable dictionary."""
    outdict: dict[str, Any] = {}
    outdict["Name"] = cons_read.name
    outdict["Consensus"] = cons_read.seq
    outdict["Members"] = dict(cons_read.json)
    return outdict


def cluster_consensus_worker(args: tuple) -> None:
    """Run UMI clustering and consensus read generation on one region."""
    (
        umi_dict,
        samplename,
        tmpfilename,
        regionid,
        contig,
        start,
        end,
        edit_distance_threshold,
        bamfilename,
        include_singletons,
        annotations,
        fasta,
        consensus_method,
        indel_frequency_cutoff,
        consensus_frequency_cutoff,
        output_json,
    ) = args  # extract args

    indel_frequency_cutoff = float(indel_frequency_cutoff)
    consensus_frequency_cutoff = float(consensus_frequency_cutoff)
    # UMI clustering
    adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
    clusters = get_connected_components(umi_dict, adj_matrix)
    umis = merge_clusters(umi_dict, clusters)

    # Consensus sequence generation
    position_matrix, singleton_matrix = get_cons_dict(
        bamfilename, umis, contig, start, end, True
    )  # include_singletons=True
    if consensus_method == "position":
        consensus_seq = get_all_consensus(
            position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
        )
    elif consensus_method == "most_common":
        consensus_seq = get_all_consensus_most_common(
            position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
        )
    elif consensus_method == "msa":
        consensus_seq = get_all_consensus_msa(
            position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff, output_json
        )

    outfilename = tmpfilename

    # Write consensus reads (and singletons) to a BAM file
    num_cons = 0
    with pysam.AlignmentFile(bamfilename, "rb") as f, pysam.AlignmentFile(outfilename, "wb", template=f) as g:  # type: ignore
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
                num_cons += 1
        if include_singletons:
            write_singleton_reads(singleton_matrix, regionid, g)

    if output_json:
        json_filename = outfilename.replace(".bam", ".json")
        outputlist = []
        for cons_read in consensus_seq.values():
            if cons_read:
                outputlist.append(write_to_json(cons_read))
        json_object = json.dumps(outputlist)
        with Path(json_filename).open("w") as f:
            f.write(json_object)

    # Generate info for cons file
    cons = get_cons_info(consensus_seq, singleton_matrix)
    consfilename = outfilename.rstrip(".bam") + ".cons"
    if len(cons) > 0:
        startpos = min(list(cons.keys()))  # take the rightmost coordinate as start
        endpos = max(list(cons.keys())) + 1  # take the leftmost coordinate as end
        with pysam.FastaFile(fasta) as fasta_file:
            ref_seq = get_reference_sequence(fasta_file, contig, startpos, endpos)
        with Path(consfilename).open("w") as g:
            write_consensus(g, cons, ref_seq, startpos, contig, annotations, samplename, False)
    else:  # empty file
        Path(consfilename).touch()


def update_bam_header(bamfile: str, samplename: str) -> dict[str, Any]:
    """Update BAM header with sample name."""
    with pysam.AlignmentFile(bamfile, "rb") as f:
        new_header = f.header.copy().to_dict()
        template = {"ID": "L1", "SM": samplename, "LB": samplename, "PL": "ILLUMINA"}

    new_header["RG"] = [template]
    return new_header


def merge_bams(output_path: str | Path, original_bamfile: str, bamfilelist: list[str], sample_name: str) -> None:
    """Merge all BAM files in bamfilelist, and remove temporary files."""
    new_header = update_bam_header(original_bamfile, sample_name)
    output_bam = Path(output_path) / f"{sample_name}_consensus_reads.bam"
    g = pysam.AlignmentFile(str(output_bam), "wb", header=new_header)
    for filename in bamfilelist:
        with pysam.AlignmentFile(filename, "rb") as f1:
            for line in f1:
                g.write(line)
    g.close()

    for filename in bamfilelist:
        Path(filename).unlink()


def merge_cons(output_path: str | Path, consfilelist: list[str], sample_name: str) -> None:
    """Merge all cons files in consfilelist and remove temporary files."""
    output_tsv = Path(output_path) / f"{sample_name}_cons.tsv"
    with output_tsv.open("w") as g:
        g.write(
            "Sample Name\tContig\tPosition\tName\tReference\tA\tC\tG\tT\tI\tD\tN\tCoverage\tConsensus group size\tMax Non-ref Allele Count\tMax Non-ref Allele Frequency\tMax Non-ref Allele\n"
        )
        for filename in consfilelist:
            with Path(filename).open() as f:
                for line in f:
                    g.write(line)

    for filename in consfilelist:
        Path(filename).unlink()


def check_duplicate_positions(cons_file: str) -> dict[str, list[str]]:
    """Check for duplicate positions in consensus file.

    Uses proper temporary files and avoids shell=True for security.
    """
    chrlist = []

    # Create temporary files in a secure way
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp1:
        tmp1_path = tmp1.name
        with Path(cons_file).open() as f:
            f.readline()
            for line in f:
                parts = line.split("\t")
                if parts[1] not in chrlist:
                    chrlist.append(parts[1])
                if parts[13] == "0":
                    tmp1.write(" ".join(parts[1:3]) + "\n")

    # Use proper subprocess chaining without shell=True
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp2:
        tmp2_path = tmp2.name

    try:
        # Find sort and uniq commands
        sort_cmd = shutil.which("sort") or "sort"
        uniq_cmd = shutil.which("uniq") or "uniq"

        with Path(tmp2_path).open("w") as outfile:
            # Chain: sort tmp1 | uniq -d > tmp2
            p1 = subprocess.Popen([sort_cmd, tmp1_path], stdout=subprocess.PIPE)
            p2 = subprocess.Popen([uniq_cmd, "-d"], stdin=p1.stdout, stdout=outfile)
            if p1.stdout:
                p1.stdout.close()  # Allow p1 to receive SIGPIPE if p2 exits
            p2.communicate()

        duppos = {}
        for chrx in chrlist:
            duppos[chrx] = []

        with Path(tmp2_path).open() as f:
            for line in f:
                line = line.rstrip()
                parts = line.split()
                if len(parts) >= 2:
                    chrx = parts[0]
                    pos = parts[1]
                    duppos[chrx].append(pos)

        return duppos

    finally:
        # Clean up temporary files
        tmp1 = Path(tmp1_path)
        tmp2 = Path(tmp2_path)
        if tmp1.exists():
            tmp1.unlink()
        if tmp2.exists():
            tmp2.unlink()


def sum_lists(*args: tuple[int, ...]) -> list[int]:
    """Sum lists elementwise."""
    return list(map(sum, zip(*args)))


def merge_duplicate_positions(args: tuple[str, list[str], str]) -> None:
    """Merge duplicate positions in consensus file."""
    chrx, duppos, cons_file = args
    dupcons = {}
    with Path(cons_file).open() as f:
        line = f.readline()
        for line in f:
            parts = line.split("\t")
            pos = parts[CONS_FILE_POS_COL]
            contig = parts[CONS_FILE_POS_COL - 1]
            if contig == chrx and pos in duppos:
                fsize = parts[CONS_FILE_FSIZE_COL]
                if pos not in dupcons:
                    dupcons[pos] = {}
                if fsize not in dupcons[pos]:
                    dupcons[pos][fsize] = []
                dupcons[pos][fsize].append(line)
    newpos = {}
    for pos in dupcons:
        newpos[pos] = {}
        for fsize in dupcons[pos]:
            newpos[pos][fsize] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for s in dupcons[pos][fsize]:
                parts = s.split("\t")
                newpos[pos][fsize] = [
                    int(x) + int(y)
                    for x, y in zip(newpos[pos][fsize], parts[CONS_FILE_ALLELE_START:CONS_FILE_ALLELE_END])
                ]
    with Path(cons_file).open() as f, Path(cons_file + "_new" + chrx).open("w") as g:
        line = f.readline()
        g.write(line)
        positions = []
        fsizes = [str(x) for x in DEFAULT_FAMILY_SIZES]
        for line in f:
            parts = line.split("\t")
            pos = parts[CONS_FILE_POS_COL]
            contig = parts[CONS_FILE_POS_COL - 1]
            if contig == chrx:
                if pos not in newpos:
                    g.write(line)
                else:
                    if pos not in positions:
                        for fsize in fsizes:
                            if fsize in newpos[pos]:
                                tmp = newpos[pos][fsize]
                                newlist = [str(x) for x in tmp]
                                consdict = {
                                    "A": tmp[0],
                                    "C": tmp[1],
                                    "G": tmp[2],
                                    "T": tmp[3],
                                    "I": tmp[4],
                                    "D": tmp[5],
                                    "N": tmp[6],
                                }
                                # tot = sum(consdict.values())
                                tot = sum(consdict[key] for key in consdict if key != "I")
                                if tot > 0:
                                    refbase = parts[4]
                                    nonrefcons = {key: consdict[key] for key in consdict if key != refbase}
                                    mna, freq, count = calc_major_nonref_allele_frequency(nonrefcons, parts[4], tot)
                                    # frac=(newpos[pos][fsize][9]/newpos[pos][fsize][7])*1.0
                                    g.write(
                                        "\t".join(parts[0:5])
                                        + "\t"
                                        + "\t".join(newlist[0:8])
                                        + "\t"
                                        + fsize
                                        + "\t"
                                        + str(count)
                                        + "\t"
                                        + str(freq)
                                        + "\t"
                                        + mna
                                        + "\n"
                                    )
                        positions.append(pos)


def merge_duplicate_positions_all_chromosomes(duppos: dict[str, list[str]], cons_file: str, num_cpus: int) -> None:
    """Merge duplicate positions across all chromosomes."""
    argvec = []
    # print(duppos)
    for chrx in duppos:
        tmpargs = (chrx, duppos[chrx], cons_file)
        argvec.append(tmpargs)
    p = Pool(int(num_cpus))
    p.map(merge_duplicate_positions, argvec)
    merge_tmp_cons_files(duppos.keys(), cons_file)
    cons_file2 = Path(cons_file + "2")
    if cons_file2.is_file():
        Path(cons_file).unlink()
        cons_file2.rename(cons_file)


def merge_tmp_cons_files(chrlist: Iterable[str], cons_file: str) -> None:
    """Merge all temporary consensus files into a single file."""
    try:
        chrlist_sorted = sorted(chrlist, key=int)
    except ValueError:
        chrlist_sorted = sorted(chrlist)
    tmpfilelist = [cons_file + "_new" + str(x) for x in chrlist_sorted]
    with Path(cons_file + "2").open("w") as g:
        i = 0
        for filename in tmpfilelist:
            with Path(filename).open() as f:
                if i > 0:
                    f.readline()
                for line in f:
                    g.write(line)
                i += 1
    for filename in tmpfilelist:
        Path(filename).unlink()


def index_bam_file(filename: str, num_threads: int = 1) -> None:
    """Index the consensus reads BAM file."""
    pysam.sort("-@", str(num_threads), filename, "-o", filename + ".sorted", catch_stdout=False)
    Path(filename + ".sorted").rename(filename)
    pysam.index(filename, catch_stdout=False)


def split_into_chunks(umi_dict: dict[str, int], clusters: list[list[str]]) -> list[dict[str, int]]:
    """Split a region into chunks if it contains more than 100000 raw reads.

    Keeps all barcodes in the same cluster together in the same chunk.
    """
    n = 0
    newdicts = []
    newdict = {}
    for c in clusters:
        for j in c:
            count = umi_dict[j]
            newdict[j] = count
            n += count
        if n > 100000:
            newdicts.append(newdict)
            newdict = {}
            n = 0
    newdicts.append(newdict)  # add remaining entries
    return newdicts


def cluster_umis_all_regions(
    regions: dict[str, dict[int | str, dict[str, int]]],
    ends: dict[str, dict[int | str, int]],
    edit_distance_threshold: int,
    samplename: str,
    bamfilename: str,
    output_path: str | Path,
    include_singletons: bool,
    fasta: str | Path,
    bedregions: dict[str, list],
    num_cpus: int,
    consensus_method: str,
    indel_frequency_cutoff: float,
    consensus_frequency_cutoff: float,
    outputjson: bool = False,
    region_from_tag: bool = False,
    starts: dict[str, dict[int | str, int]] | None = None,
) -> list[str]:
    """Run UMI clustering and error correction using num_cpus threads.

    Processes one region on each thread.
    """
    if starts is None:
        starts = {}
    argvec = []
    bamfilelist = []
    i = 0
    j = 0

    for contig in regions:
        for pos in regions[contig]:
            annotations = bedregions.get(contig, [])

            if region_from_tag:
                i = pos
                posx = int(starts[contig][pos])
                j = 0
            else:
                posx = int(pos)
            tmpfilename = f"{output_path}/tmp_{i}.bam"  # noqa: S108
            numreads = sum(regions[contig][pos].values())
            if numreads > 100000:  # split in chunks
                umi_dict = regions[contig][pos]
                adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
                clusters = get_connected_components(umi_dict, adj_matrix)
                newdicts = split_into_chunks(umi_dict, clusters)
                for x in newdicts:
                    tmpfilename = f"{output_path}/tmp_{i}.bam"  # noqa: S108
                    argvec.append(
                        (
                            x,
                            samplename,
                            tmpfilename,
                            i,
                            contig,
                            posx,
                            int(ends[contig][pos]),
                            int(edit_distance_threshold),
                            bamfilename,
                            include_singletons,
                            annotations,
                            fasta,
                            consensus_method,
                            indel_frequency_cutoff,
                            consensus_frequency_cutoff,
                            outputjson,
                        )
                    )
                    bamfilelist.append(f"{output_path}/tmp_{i}.bam")  # noqa: S108
                    if not region_from_tag:
                        i += 1
                    else:
                        i = i + "_" + str(j)
                        j += 1
            else:
                argvec.append(
                    (
                        regions[contig][pos],
                        samplename,
                        tmpfilename,
                        i,
                        contig,
                        posx,
                        int(ends[contig][pos]),
                        int(edit_distance_threshold),
                        bamfilename,
                        include_singletons,
                        annotations,
                        fasta,
                        consensus_method,
                        indel_frequency_cutoff,
                        consensus_frequency_cutoff,
                        outputjson,
                    )
                )
                bamfilelist.append(f"{output_path}/tmp_{i}.bam")  # noqa: S108
                if not region_from_tag:
                    i += 1

    p = Pool(int(num_cpus))

    p.map(cluster_consensus_worker, argvec)
    return bamfilelist


def cluster_umis_on_position(
    bamfilename: str, position_threshold: int, group_method: str, bedfilename: str | None = None
) -> tuple:
    """Cluster UMIs by genomic position."""
    position_threshold = int(position_threshold)
    if group_method == "fromBed":
        regions, ends = read_bam_from_bed(bamfilename, bedfilename, position_threshold)
    elif group_method == "fromTag":
        regions, starts, ends = read_bam_from_tag(bamfilename)
    else:
        regions, ends = readBam(bamfilename, position_threshold)

    if group_method == "fromTag":
        return (regions, ends, starts)
    else:
        return (regions, ends)


def run_umi_errorcorrect(config: UMIErrorCorrectConfig) -> None:
    """Run UMI clustering and consensus read generation (error correction).

    Args:
        config: Configuration object containing all parameters for error correction.
    """
    logger.info("Starting UMI clustering")

    # Derive group_method from config flags
    if config.regions_from_bed:
        group_method = "fromBed"
    elif config.regions_from_tag:
        group_method = "fromTag"
    else:
        group_method = "automatic"

    # Validate consensus method
    consensus_method = config.consensus_method.lower()
    if consensus_method.startswith("pos"):
        consensus_method = "position"
    elif consensus_method.startswith("most"):
        consensus_method = "most_common"
    elif consensus_method.startswith("msa") or consensus_method.startswith("multiple"):
        consensus_method = "msa"
    else:
        raise ValueError("Please choose consensus method between 'position', 'most_common', 'msa'")

    logger.info(f"Group by position method: {group_method}")
    logger.info(f"Consensus method: {consensus_method}")

    # Config validator handles bam_file and sample_name auto-detection
    output_path = config.output_path
    bam_file = str(config.bam_file) if config.bam_file else None
    sample_name = config.sample_name
    bed_file = str(config.bed_file) if config.bed_file else None

    if not bam_file:
        raise ValueError("bam_file must be provided or auto-detected")
    if not sample_name:
        raise ValueError("sample_name must be provided or auto-detected")

    # Cluster UMIs on position
    if group_method == "fromTag":
        regions, ends, starts = cluster_umis_on_position(bam_file, config.position_threshold, group_method, bed_file)
    else:
        regions, ends = cluster_umis_on_position(bam_file, config.position_threshold, group_method, bed_file)
        starts = None

    nregions = sum(len(regions[chrx]) for chrx in regions)
    logger.info(f"Number of regions: {nregions}")

    edit_distance_threshold = config.edit_distance_threshold
    num_cpus = config.num_threads if config.num_threads else cpu_count()
    logger.info("Starting Consensus sequence generation")
    logger.info(f"Starting {num_cpus} threads")

    fasta = str(config.reference_file)

    # Load BED regions if provided
    if bed_file:
        bedregions = read_bed(bed_file)
        bedregions = sort_regions(bedregions)
        if group_method == "fromBed":
            bedregions = merge_regions(bedregions, 0)
    else:
        bedregions = {}

    # Run clustering on all regions
    if group_method == "fromTag":
        bamfilelist = cluster_umis_all_regions(
            regions,
            ends,
            edit_distance_threshold,
            sample_name,
            bam_file,
            str(output_path),
            config.include_singletons,
            fasta,
            bedregions,
            num_cpus,
            consensus_method,
            config.indel_frequency_threshold,
            config.consensus_frequency_threshold,
            config.output_json,
            config.regions_from_tag,
            starts,
        )
    else:
        bamfilelist = cluster_umis_all_regions(
            regions,
            ends,
            edit_distance_threshold,
            sample_name,
            bam_file,
            str(output_path),
            config.include_singletons,
            fasta,
            bedregions,
            num_cpus,
            consensus_method,
            config.indel_frequency_threshold,
            config.consensus_frequency_threshold,
            config.output_json,
        )

    # Merge and index BAM files
    merge_bams(output_path, bam_file, bamfilelist, sample_name)
    consensus_bam = output_path / f"{sample_name}_consensus_reads.bam"
    index_bam_file(str(consensus_bam), num_cpus)

    # Merge consensus files
    consfilelist = [x.rstrip(".bam") + ".cons" for x in bamfilelist]
    merge_cons(output_path, consfilelist, sample_name)
    cons_file = str(output_path / f"{sample_name}_cons.tsv")

    # Optionally remove large files
    if config.remove_large_files:
        bam_path = output_path / Path(bam_file).name
        if bam_path.exists():
            bam_path.unlink()

    # Handle duplicate positions in cons file
    duppos = check_duplicate_positions(cons_file)
    if any(duppos.values()):
        merge_duplicate_positions_all_chromosomes(duppos, cons_file, num_cpus)

    # Generate stats file from consensus BAM (single source of truth)
    from umierrorcorrect2.get_consensus_statistics import get_stat, write_stats_file

    stats = get_stat(consensus_bam, bed_file)
    write_stats_file(stats, output_path, sample_name)

    logger.info(
        f"Consensus generation complete, output written to {output_path}/{sample_name}_consensus_reads.bam, "
        f"{output_path}/{sample_name}_cons.tsv"
    )


def main(config: UMIErrorCorrectConfig) -> None:
    """Main entry point for UMI error correction."""
    run_umi_errorcorrect(config)
