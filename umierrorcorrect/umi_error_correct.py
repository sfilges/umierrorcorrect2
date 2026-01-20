#!/usr/bin/env python
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pysam

from umierrorcorrect.core.consensus import (
    get_all_consensus,
    get_all_consensus_most_common,
    get_all_consensus_msa,
    get_cons_dict,
    get_reference_sequence,
    write_singleton_reads,
)
from umierrorcorrect.core.get_cons_info import calc_major_nonref_allele_frequency, get_cons_info, write_consensus
from umierrorcorrect.core.get_regions_from_bed import get_overlap, merge_regions, read_bed, sort_regions
from umierrorcorrect.core.group import read_bam_from_bed, read_bam_from_tag, readBam
from umierrorcorrect.core.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from umierrorcorrect.core.utils import check_output_directory, get_sample_name


def write_to_json(cons_read):
    outdict = {}
    outdict["Name"] = cons_read.name
    outdict["Consensus"] = cons_read.seq
    outdict["Members"] = dict(cons_read.json)
    return outdict


def cluster_consensus_worker(args):
    """Run UMI clustering and consensus read generation on one region"""
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
    with pysam.AlignmentFile(bamfilename, "rb") as f, pysam.AlignmentFile(outfilename, "wb", template=f) as g:
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
    statfilename = outfilename.rstrip(".bam") + ".hist"
    if len(cons) > 0:
        startpos = min(list(cons.keys()))  # take the rightmost coordinate as start
        endpos = max(list(cons.keys())) + 1  # take the leftmost coordinate as end
        with pysam.FastaFile(fasta) as f:
            ref_seq = get_reference_sequence(f, contig, startpos, endpos)
        with Path(consfilename).open("w") as g:
            write_consensus(g, cons, ref_seq, startpos, contig, annotations, samplename, False)
    else:  # empty file
        Path(consfilename).touch()
    # Write to hist/stat file
    if len(cons) > 0:
        with Path(statfilename).open("w") as g2:
            name = get_overlap(annotations, start, endpos)
            regionname = f"{contig}:{start}-{endpos}"
            g2.write(
                "\t".join(
                    [
                        str(regionid),
                        regionname,
                        name,
                        "consensus_reads: " + str(num_cons),
                        "singletons: " + str(len(singleton_matrix)),
                    ]
                )
                + "\n"
            )
    else:  # empty file
        Path(statfilename).touch()


def update_bam_header(bamfile, samplename):
    with pysam.AlignmentFile(bamfile, "rb") as f:
        new_header = f.header.copy().to_dict()
        template = {"ID": "L1", "SM": samplename, "LB": samplename, "PL": "ILLUMINA"}

    new_header["RG"] = [template]
    return new_header


def merge_bams(output_path, original_bamfile, bamfilelist, sample_name):
    """Merge all BAM files for in bamfilelist, and remove temporary files"""
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


def merge_cons(output_path, consfilelist, sample_name):
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


def check_duplicate_positions(cons_file):
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


def sum_lists(*args):
    return list(map(sum, zip(*args)))


def merge_duplicate_positions(args):
    chrx, duppos, cons_file = args
    dupcons = {}
    a = 13
    b = 2
    with Path(cons_file).open() as f:
        line = f.readline()
        for line in f:
            parts = line.split("\t")
            pos = parts[b]
            contig = parts[b - 1]
            if contig == chrx and pos in duppos:
                fsize = parts[a]
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
                newpos[pos][fsize] = [int(a) + int(b) for a, b in zip(newpos[pos][fsize], parts[5:15])]
    with Path(cons_file).open() as f, Path(cons_file + "_new" + chrx).open("w") as g:
        line = f.readline()
        g.write(line)
        positions = []
        fsizes = ["0", "1", "2", "3", "4", "5", "7", "10", "20", "30"]
        for line in f:
            parts = line.split("\t")
            pos = parts[b]
            contig = parts[b - 1]
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

    # os.remove(cons_file)
    # os.rename(cons_file+'_new',cons_file)


def merge_duplicate_positions_all_chromosomes(duppos, cons_file, num_cpus):
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


def merge_tmp_cons_files(chrlist, cons_file):
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


def merge_stat(output_path, statfilelist, sample_name):
    """Merge all stat files in statfilelist and remove temporary files."""
    output_hist = Path(output_path) / f"{sample_name}.hist"
    with output_hist.open("w") as g:
        for filename in statfilelist:
            with Path(filename).open() as f:
                for line in f:
                    g.write(line)

    for filename in statfilelist:
        Path(filename).unlink()


def merge_duplicate_stat(output_path, samplename):
    convert = lambda text: int(text) if text.isdigit() else text.lower()  # noqa: E731
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]  # noqa: E731
    histfile = Path(output_path) / f"{samplename}.hist"
    regions = {}
    # histfile
    with histfile.open() as f:
        for line in f:
            line = line.rstrip()
            parts = line.split("\t")
            idx = parts[0]
            chrx = parts[1].split(":")[0]
            if chrx not in regions:
                regions[chrx] = {}
            pos = parts[1].split(":")[1].split("-")[0]
            end = parts[1].split(":")[1].split("-")[1]
            name = parts[2]
            numcons = parts[3].split()[1]
            numsing = parts[4].split()[1]
            if pos not in regions[chrx]:
                regions[chrx][pos] = (idx, int(pos), int(end), name, int(numcons), int(numsing))
            else:
                tmp = regions[chrx][pos]
                if "-" in tmp[0]:
                    newid = tmp[0].split("-")[0] + "-" + idx
                else:
                    newid = tmp[0] + "-" + idx
                if int(end) > int(tmp[2]):
                    newend = end
                else:
                    newend = tmp[2]
                newnumcons = tmp[4] + int(numcons)
                newnumsing = tmp[5] + int(numsing)
                regions[chrx][pos] = (newid, tmp[1], newend, name, newnumcons, newnumsing)
    histfile2 = Path(str(histfile) + "2")
    with histfile2.open("w") as g:
        for chrx in sorted(regions, key=alphanum_key):
            for pos in regions[chrx]:
                tmp = regions[chrx][pos]
                g.write(
                    f"{tmp[0]}\t{str(chrx)}:{tmp[1]}-{tmp[2]}\t{tmp[3]}\tconsensus_reads: {tmp[4]}\tsingletons: {tmp[5]}\n"
                )
    histfile2.rename(histfile)


def index_bam_file(filename, num_threads=1):
    """Index the consensus reads bam file"""
    pysam.sort("-@", str(num_threads), filename, "-o", filename + ".sorted", catch_stdout=False)
    Path(filename + ".sorted").rename(filename)
    pysam.index(filename, catch_stdout=False)


def split_into_chunks(umi_dict, clusters):
    """If one region contains more than 100000 raw reads, split in chunks of 100000.
    keep all barcodes in same cluster in the same chunk.
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

    # n = 0
    # i = 0
    # newdicts = []
    # b=list(dictname.values())
    # a=iter(dictname.items())
    # for count in b:
    #    n += count
    #    i += 1
    #    if n > 100000:
    #        newdicts.append(dict(islice(a,i))) #add (more than) 100000 raw reads to newdicts (from 0 to index i)
    #        n = 0
    #        i = 0
    # newdicts.append(dict(islice(a,i))) #add remaining entries
    # return(newdicts)


def cluster_umis_all_regions(
    regions,
    ends,
    edit_distance_threshold,
    samplename,
    bamfilename,
    output_path,
    include_singletons,
    fasta,
    bedregions,
    num_cpus,
    consensus_method,
    indel_frequency_cutoff,
    consensus_frequency_cutoff,
    outputjson=False,
    region_from_tag=False,
    starts=None,
):
    """Function for running UMI clustering and error correction using num_cpus threads,
    i.e. one region on each thread."""
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


def cluster_umis_on_position(bamfilename, position_threshold, group_method, bedfilename=None):
    """Function for 0cluster umis on position"""
    # position_threshold = 20
    # group_method='fromBed'
    # group_method='automatic'
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


def run_umi_errorcorrect(args):
    """Run the umi clustering and consensus read generation (error correction)"""
    logging.info("Starting UMI clustering")
    args.output_path = check_output_directory(args.output_path)

    if args.regions_from_bed:
        group_method = "fromBed"
    elif args.regions_from_tag:
        group_method = "fromTag"
    else:
        group_method = "automatic"
    args.consensus_method = args.consensus_method.lower()
    if args.consensus_method.startswith("pos"):
        consensus_method = "position"
    elif args.consensus_method.startswith("most"):
        consensus_method = "most_common"
    elif args.consensus_method.startswith("msa") or args.consensus_method.startswith("multiple"):
        consensus_method = "msa"
    else:
        print("Please choose consensus method between 'position','most_common','MSA'")
        sys.exit(1)
    logging.info(f"Group by position method: {group_method}")
    logging.info(f"Consensus method: {consensus_method}")
    output_path = Path(args.output_path)
    if not args.bam_file and args.sample_name:
        # see if it is possible to guess bam file from previous step.
        testname = output_path / f"{args.sample_name}.sorted.bam"
        if testname.is_file():
            args.bam_file = str(testname)
    if not args.bam_file:
        bamfile = list(output_path.glob("*sorted.bam"))
        if len(bamfile) > 1:
            print(
                "Too many sorted.bam files in output folder, please specify which sample to run with -s (sample name) or -b (path to bam file)."
            )
            sys.exit(1)
        else:
            args.bam_file = str(bamfile[0])
    if not args.sample_name:
        args.sample_name = get_sample_name(args.bam_file)
    if group_method == "fromTag":
        regions, ends, starts = cluster_umis_on_position(
            args.bam_file, args.position_threshold, group_method, args.bed_file
        )
    else:
        regions, ends = cluster_umis_on_position(args.bam_file, args.position_threshold, group_method, args.bed_file)

    # print(regions)
    nregions = 0
    for chrx in regions:
        nregions += len(regions[chrx])
    logging.info(f"Number of regions, {nregions}")

    edit_distance_threshold = args.edit_distance_threshold
    num_cpus = int(args.num_threads) if args.num_threads else int(cpu_count())
    logging.info("Starting Consensus sequence generation")
    logging.info(f"Starting {num_cpus} threads")
    fasta = args.reference_file
    if args.bed_file:
        bedregions = read_bed(args.bed_file)
        bedregions = sort_regions(bedregions)
        if group_method == "fromBed":
            bedregions = merge_regions(bedregions, 0)
    else:
        bedregions = {}
    # print(bedregions)
    if group_method == "fromTag":
        bamfilelist = cluster_umis_all_regions(
            regions,
            ends,
            edit_distance_threshold,
            args.sample_name,
            args.bam_file,
            args.output_path,
            args.include_singletons,
            fasta,
            bedregions,
            num_cpus,
            consensus_method,
            args.indel_frequency_threshold,
            args.consensus_frequency_threshold,
            args.output_json,
            args.regions_from_tag,
            starts,
        )
    else:
        bamfilelist = cluster_umis_all_regions(
            regions,
            ends,
            edit_distance_threshold,
            args.sample_name,
            args.bam_file,
            args.output_path,
            args.include_singletons,
            fasta,
            bedregions,
            num_cpus,
            consensus_method,
            args.indel_frequency_threshold,
            args.consensus_frequency_threshold,
            args.output_json,
        )
    merge_bams(args.output_path, args.bam_file, bamfilelist, args.sample_name)
    consensus_bam = output_path / f"{args.sample_name}_consensus_reads.bam"
    index_bam_file(str(consensus_bam), num_cpus)
    consfilelist = [x.rstrip(".bam") + ".cons" for x in bamfilelist]
    merge_cons(args.output_path, consfilelist, args.sample_name)
    cons_file = str(output_path / f"{args.sample_name}_cons.tsv")
    if args.remove_large_files:
        (output_path / args.bam_file).unlink()

    statfilelist = [x.rstrip(".bam") + ".hist" for x in bamfilelist]
    merge_stat(args.output_path, statfilelist, args.sample_name)
    duppos = check_duplicate_positions(cons_file)
    if any(duppos):
        merge_duplicate_positions_all_chromosomes(duppos, cons_file, num_cpus)
    merge_duplicate_stat(args.output_path, args.sample_name)
    logging.info(
        f"Consensus generation complete, output written to {args.output_path}/{args.sample_name}_consensus_reads.bam, "
        f"{args.output_path}/{args.sample_name}_cons.tsv"
    )


def main(args):
    run_umi_errorcorrect(args)
