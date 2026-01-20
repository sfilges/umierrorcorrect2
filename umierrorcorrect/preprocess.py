#!/usr/bin/env python3
"""
UMI error correct, preprocess.py - remove UMI and append to the header of Fastq sequences.
================================

:Author: Tobias Osterlund

Purpose
-------

Preprocess the fastq files by removing the unique molecular index and add it to the header of the fastq entry.

"""

import logging
import subprocess
import sys
from pathlib import Path

from umierrorcorrect.core.check_args import check_args_fastq
from umierrorcorrect.core.logging_config import log_subprocess_stderr
from umierrorcorrect.core.read_fastq_records import read_fastq, read_fastq_paired_end
from umierrorcorrect.core.utils import check_output_directory


def generate_random_dir(tmpdir):
    """Generate a directory for storing temporary files, using a timestamp."""
    import datetime

    newtmpdir = Path(tmpdir) / f"r{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}"
    newtmpdir = check_output_directory(str(newtmpdir))
    return newtmpdir


def run_unpigz(filename, tmpdir, num_threads, program):
    """Unzip the fastq.gz files using parallel gzip (pigz)."""
    input_path = Path(filename)
    outfilename = Path(tmpdir) / input_path.name.removesuffix(".gz")
    if program == "pigz":
        command = ["unpigz", "-p", num_threads, "-c", filename]
    elif program == "gzip":
        command = ["gunzip", "-c", filename]
    with outfilename.open("w") as g:
        p = subprocess.Popen(command, stdout=g, stderr=subprocess.PIPE)
        _, stderr = p.communicate()
        log_subprocess_stderr(stderr, program)
        p.wait()
    return str(outfilename)


def run_gunzip(filename, tmpdir):
    """Unzip the fastq.gz files using parallel gzip (pigz)."""
    input_path = Path(filename)
    outfilename = Path(tmpdir) / input_path.name.removesuffix(".gz")
    command = ["gunzip", "-c", filename]
    with outfilename.open("w") as g:
        p = subprocess.Popen(command, stdout=g, stderr=subprocess.PIPE)
        _, stderr = p.communicate()
        log_subprocess_stderr(stderr, "gunzip")
        p.wait()
    return str(outfilename)


def run_pigz(filename, num_threads, program):
    """Zip the new fastq files with parallel gzip (pigz)."""
    if program == "pigz":
        command = ["pigz", "-p", num_threads, filename]
    elif program == "gzip":
        command = ["gzip", filename]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = p.communicate()
    log_subprocess_stderr(stderr, program)
    p.wait()


def preprocess_se(infilename, outfilename, barcode_length, spacer_length):
    """Run the preprocessing for single end data (one fastq file)."""
    with Path(infilename).open() as f, Path(outfilename).open("w") as g:
        read_start = barcode_length + spacer_length
        nseqs = 0
        for name, seq, qual in read_fastq(f):
            nseqs += 1
            barcode = seq[:barcode_length]
            # g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            parts = name.split()
            newname = ":".join([parts[0], barcode]) + " " + parts[-1]
            g.write("\n".join([newname, seq[read_start:], "+", qual[read_start:]]) + "\n")
    return nseqs


def preprocess_pe(r1file, r2file, outfile1, outfile2, barcode_length, spacer_length, dual_index):
    """Run the preprocessing for paired end data (two fastq files)."""
    read_start = barcode_length + spacer_length
    with (
        Path(r1file).open() as f1,
        Path(r2file).open() as f2,
        Path(outfile1).open("w") as g1,
        Path(outfile2).open("w") as g2,
    ):
        nseqs = 0
        for name1, seq1, qual1, name2, seq2, qual2 in read_fastq_paired_end(f1, f2):
            nseqs += 1
            if dual_index:
                barcode = seq1[:barcode_length] + seq2[:barcode_length]
            else:
                barcode = seq1[:barcode_length]
            parts1 = name1.split()
            parts2 = name2.split()
            newname1 = ":".join([parts1[0], barcode]) + " " + parts1[-1]
            newname2 = ":".join([parts2[0], barcode]) + " " + parts2[-1]
            g1.write("\n".join([newname1, seq1[read_start:], "+", qual1[read_start:]]) + "\n")
            if dual_index:
                g2.write("\n".join([newname2, seq2[read_start:], "+", qual2[read_start:]]) + "\n")
            else:
                g2.write("\n".join([newname2, seq2, "+", qual2]) + "\n")
    return 2 * nseqs


def run_preprocessing(args):
    """Start preprocessing."""
    logging.info(f"Start preprocessing of sample {args.sample_name}")

    if args.tmpdir:
        newtmpdir = generate_random_dir(args.tmpdir)
    else:
        newtmpdir = generate_random_dir(args.output_path)
    # args.chunksize=int(args.chunksize)
    # Unzip the fastq.gz files
    if not args.read1.endswith("gz"):
        r1file = args.read1
        removerfiles = False
        if args.mode == "paired":
            r2file = args.read2
    else:
        removerfiles = True
        if args.mode == "paired":
            r1file = run_unpigz(args.read1, newtmpdir, args.num_threads, args.gziptool)
            r2file = run_unpigz(args.read2, newtmpdir, args.num_threads, args.gziptool)
        else:
            r1file = run_unpigz(args.read1, newtmpdir, args.num_threads, args.gziptool)

    logging.info(f"Writing output files to {args.output_path}")
    if args.adapter_trimming is True:
        if args.adapter_sequence.lower() == "illumina":
            adapter = "AGATCGGAAGAGC"
        elif args.adapter_sequence.lower() == "nextera":
            adapter = "CTGTCTCTTATA"
        elif args.adapter_sequence.lower() == "small-rna":
            adapter = "ATGGAATTCTCG"
        else:
            adapter = args.adapter_sequence.upper()

        output_path = Path(args.output_path)
        if args.mode == "single":
            outfilename = str(output_path / f"{args.sample_name}_trimmed.fastq")
            command = ["cutadapt", "-a", adapter, "-o", outfilename, "-O", "3", "-m", "20", r1file]
        else:
            outfile1 = str(output_path / f"{args.sample_name}_R1_trimmed.fastq")
            outfile2 = str(output_path / f"{args.sample_name}_R2_trimmed.fastq")
            command = [
                "cutadapt",
                "-a",
                adapter,
                "-A",
                adapter,
                "-o",
                outfile1,
                "-p",
                outfile2,
                "-O",
                "3",
                "-m",
                "20",
                r1file,
                r2file,
            ]
        logging.info(f"Performing adapter trimming using cutadapt with adapter sequence {adapter}")
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, stderr = p.communicate()
        log_subprocess_stderr(stderr, "cutadapt")
        p.wait()
        if args.mode == "single":
            Path(r1file).unlink()
            r1file = outfilename
        else:
            Path(r1file).unlink()
            Path(r2file).unlink()
            r1file = outfile1
            r2file = outfile2

    output_path = Path(args.output_path)
    if args.mode == "single":
        outfilename = str(output_path / f"{args.sample_name}_umis_in_header.fastq")
        nseqs = preprocess_se(r1file, outfilename, args.umi_length, args.spacer_length)
        run_pigz(outfilename, args.num_threads, args.gziptool)
        Path(r1file).unlink()
        Path(newtmpdir).rmdir()
        fastqfiles = [f"{outfilename}.gz"]
    else:
        if args.reverse_index:
            # switch forward and reverse read
            r1filetmp = r1file
            r1file = r2file
            r2file = r1filetmp
            outfile1 = str(output_path / f"{args.sample_name}_R2_umis_in_header.fastq")
            outfile2 = str(output_path / f"{args.sample_name}_R1_umis_in_header.fastq")
        else:
            # r1file=args.read1
            # r2file=args.read2
            outfile1 = str(output_path / f"{args.sample_name}_R1_umis_in_header.fastq")
            outfile2 = str(output_path / f"{args.sample_name}_R2_umis_in_header.fastq")
        nseqs = preprocess_pe(r1file, r2file, outfile1, outfile2, args.umi_length, args.spacer_length, args.dual_index)
        run_pigz(outfile1, args.num_threads, args.gziptool)
        run_pigz(outfile2, args.num_threads, args.gziptool)
        r1path = Path(r1file)
        r2path = Path(r2file)
        if removerfiles is True and r1path.is_file():
            r1path.unlink()
        if removerfiles is True and r2path.is_file():
            r2path.unlink()
        Path(newtmpdir).rmdir()
        fastqfiles = [outfile1 + ".gz", outfile2 + ".gz"]
    logging.info("Finished preprocessing")
    return (fastqfiles, nseqs)


def main(args):
    try:
        args = check_args_fastq(args)  # check if combination of arguments are correct
    except ValueError as e:
        print(e)
        sys.exit(1)
    fastqfiles, nseqs = run_preprocessing(args)
    print(nseqs)
