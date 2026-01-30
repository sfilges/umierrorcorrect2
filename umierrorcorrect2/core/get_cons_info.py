from collections import Counter

from umierrorcorrect2.core.constants import DEFAULT_FAMILY_SIZES, SINGLETON_FAMILY_SIZES
from umierrorcorrect2.core.get_regions_from_bed import get_all_annotations


def get_cons_info(consensus_seq, singletons, fsizes=DEFAULT_FAMILY_SIZES):
    """loop through the consensus reads to collapse alleles for each position"""
    cons = {}
    for consensus_read in consensus_seq.values():
        if consensus_read:
            pos = consensus_read.start_pos
            count = consensus_read.count
            if consensus_read.indel_read == 0:
                for base in consensus_read.seq:
                    if pos not in cons:
                        cons[pos] = {}
                    for fsize in fsizes:
                        if fsize == 0:
                            if fsize not in cons[pos]:
                                cons[pos][fsize] = Counter()
                            if base not in "N":
                                cons[pos][fsize][base] += count
                        elif count >= fsize:
                            if fsize not in cons[pos]:
                                cons[pos][fsize] = Counter()
                            if base not in "N":
                                cons[pos][fsize][base] += 1
                    pos += 1
            else:
                i = 0
                skipbase = []
                cigar = consensus_read.cigarstring
                for base in consensus_read.seq:
                    c = cigar[i]
                    if i not in skipbase:
                        if pos not in cons:
                            cons[pos] = {}
                        if c == "0":  # match or mismatch
                            for fsize in fsizes:
                                if fsize == 0:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    if not base == "N":
                                        cons[pos][fsize][base] += count
                                elif count >= fsize:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    if not base == "N":
                                        cons[pos][fsize][base] += 1
                            pos += 1
                            i += 1
                        elif c == "1":  # insertion
                            if cigar[i + 1] == "1":  # if next base also have an insertion
                                j = 0
                                insertion = True
                                while insertion:
                                    j += 1
                                    if cigar[i + j] == "1":
                                        skipbase.append(i + j)
                                    else:
                                        insertion = False
                                # insstring=consensus_read.seq[i-1:i+j-1]
                            # else:
                            #    insstring=consensus_read.seq[i-1]
                            for fsize in fsizes:
                                if fsize == 0:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    cons[pos][fsize]["I"] += count
                                elif count >= fsize:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    cons[pos][fsize]["I"] += 1
                            i += 1
                        elif c == "2":  # deletion
                            for fsize in fsizes:
                                if fsize == 0:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    cons[pos][fsize]["D"] += count

                                elif count >= fsize:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    cons[pos][fsize]["D"] += 1
                            deletion = False
                            if cigar[i + 1] == "2":
                                deletion = True
                                while deletion:
                                    i += 1
                                    c = cigar[i]
                                    pos += 1
                                    # if pos not in cons:
                                    #    cons[pos] = {}
                                    # for fsize in fsizes:
                                    #    if fsize == 0:
                                    #        if fsize not in cons[pos]:
                                    #            cons[pos][fsize] = Counter()
                                    #        cons[pos][fsize]['D'] += count

                                    #    elif count >= fsize:
                                    #        if fsize not in cons[pos]:
                                    #            cons[pos][fsize] = Counter()
                                    #        cons[pos][fsize]['D'] += 1
                                    if cigar[i + 1] == "2":
                                        deletion = True
                                    else:
                                        deletion = False
                            pos += 1
                            if pos not in cons:
                                cons[pos] = {}
                            for fsize in fsizes:
                                if fsize == 0:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    if base not in "N":
                                        cons[pos][fsize][base] += count

                                elif count >= fsize:
                                    if fsize not in cons[pos]:
                                        cons[pos][fsize] = Counter()
                                    if base not in "N":
                                        cons[pos][fsize][base] += 1
                            pos += 1
                            i += 1
                    else:
                        i += 1

    for read in singletons.values():
        if "I" not in read.cigarstring and "D" not in read.cigarstring:
            sequence = read.query_sequence
            for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                base = sequence[qpos]
                if refpos not in cons:
                    cons[refpos] = {}
                for fsize in SINGLETON_FAMILY_SIZES:
                    if fsize not in cons[refpos]:
                        cons[refpos][fsize] = Counter()
                    cons[refpos][fsize][base] += 1

        else:  # insertion or deletion
            positions = read.get_aligned_pairs(matches_only=True)
            q, ref = positions[0]
            q = q - 1
            ref = ref - 1
            for qpos, refpos in positions:
                if not qpos == q + 1:
                    # allele = read.query_sequence[q+1:qpos]
                    # inspos=refpos-1
                    if refpos not in cons:
                        cons[refpos] = {}
                    for fsize in SINGLETON_FAMILY_SIZES:
                        if fsize not in cons[refpos]:
                            cons[refpos][fsize] = Counter()
                        cons[refpos][fsize]["I"] += 1
                        # cons[refpos][fsize][base] += 1
                elif not refpos == ref + 1:
                    # deletion
                    dellength = refpos - (ref + 1)
                    delpos = refpos - dellength
                    if delpos not in cons:
                        cons[delpos] = {}
                    for fsize in SINGLETON_FAMILY_SIZES:
                        if fsize not in cons[delpos]:
                            cons[delpos][fsize] = Counter()
                        cons[delpos][fsize]["D"] += 1
                sequence = read.query_sequence
                base = sequence[qpos]
                if refpos not in cons:
                    cons[refpos] = {}
                for fsize in SINGLETON_FAMILY_SIZES:
                    if fsize not in cons[refpos]:
                        cons[refpos][fsize] = Counter()
                    cons[refpos][fsize][base] += 1
                else:
                    sequence = read.query_sequence
                    # qual = read.qual
                    base = sequence[qpos]
                    if refpos not in cons:
                        cons[refpos] = {}
                    for fsize in SINGLETON_FAMILY_SIZES:
                        if fsize not in cons[refpos]:
                            cons[refpos][fsize] = Counter()
                        cons[refpos][fsize][base] += 1
                q = qpos
                ref = refpos

    return cons


def calc_major_nonref_allele_frequency(cons, ref, tot):
    comp = cons
    allele = max(comp, key=comp.get)
    count = cons[allele]

    if tot > 0:
        frac = 1.0 * (cons[allele] / tot)
    else:
        frac = 0

    if frac > 0:
        return (allele, frac, count)
    else:
        return ("", 0, 0)


def write_consensus(f, cons, ref_seq, start, contig, annotation, samplename, only_target_regions):
    bases = ["A", "C", "G", "T", "I", "D", "N"]
    # print(list(cons.keys())[0],list(cons.keys())[-1],start,len(ref_seq))

    for pos in sorted(cons):
        annotation_pos = get_all_annotations(annotation, pos + 1)
        if not (annotation_pos == "" and only_target_regions):
            # if len(ref_seq)<(pos-start+1):
            #     print("error",contig,start,ref_seq)
            refbase = ref_seq[pos - start]

            for fsize in cons[pos]:
                line = []
                line.append(samplename)
                line.append(contig)
                line.append(str(pos + 1))
                line.append(annotation_pos)
                line.append(refbase)
                consline = cons[pos][fsize]
                tot = sum(consline[key] for key in consline if key != "I")
                nonrefcons = {key: consline[key] for key in consline if key != refbase}
                if len(nonrefcons) > 0:
                    mna, freq, count = calc_major_nonref_allele_frequency(nonrefcons, refbase, tot)
                else:
                    mna = ""
                    freq = 0
                    count = 0
                for base in bases:
                    if base in cons[pos][fsize]:
                        line.append(str(cons[pos][fsize][base]))
                    else:
                        line.append(str(0))
                line.append(str(tot))
                line.append(str(fsize))
                line.append(str(count))
                line.append(str(freq))
                line.append(mna)
                f.write("\t".join(line) + "\n")
