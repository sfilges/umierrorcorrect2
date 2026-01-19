#!/usr/bin/env python3
# from collections import Counter
# import time
import itertools

# import sys
import pickle
from collections.abc import Generator

from umierrorcorrect.src.constants import SUBSTRING_OPTIMIZATION_THRESHOLD


class umi_cluster:
    """Container for UMI cluster information."""

    def __init__(self, name: str, count: int) -> None:
        self.centroid: str = name
        self.count: int = count

    def add_count(self, newcount: int) -> None:
        self.count += newcount

    def change_centroid(self, newname: str) -> None:
        self.centroid = newname


def hamming_distance(a: str, b: str) -> int:
    """Returns the Hamming distance between two strings of equal length."""
    try:
        assert len(a) == len(b)
        return sum(i != j for i, j in zip(a, b))
    except AssertionError as error:
        print(f"Barcode lengths are not equal for {a}. {b}")
        raise error


def create_substring_matrix(
    barcodedict: dict[str, int], edit_distance_threshold: int
) -> list[dict[str, list[str]]]:
    """Divide each barcode in two or three substrings of (approximately) equal length."""
    edit_distance_threshold = int(edit_distance_threshold)
    umi_length = len(list(barcodedict.keys())[0])
    if edit_distance_threshold <= 1:
        s = round(umi_length // 2)
        substr_dict1 = {}
        substr_dict2 = {}
        for barcode in barcodedict:
            sub1 = barcode[:s]
            sub2 = barcode[s:]
            if sub1 not in substr_dict1:
                substr_dict1[sub1] = []
            if sub2 not in substr_dict2:
                substr_dict2[sub2] = []
            substr_dict1[sub1].append(barcode)
            substr_dict2[sub2].append(barcode)
        return [substr_dict1, substr_dict2]
    if edit_distance_threshold == 2:
        s = round(umi_length // 3)
        substr_dict1 = {}
        substr_dict2 = {}
        substr_dict3 = {}
        for barcode in barcodedict:
            sub1 = barcode[:s]
            sub2 = barcode[s : 2 * s]
            sub3 = barcode[2 * s :]
            if sub1 not in substr_dict1:
                substr_dict1[sub1] = []
            if sub2 not in substr_dict2:
                substr_dict2[sub2] = []
            if sub3 not in substr_dict3:
                substr_dict3[sub3] = []
            substr_dict1[sub1].append(barcode)
            substr_dict2[sub2].append(barcode)
            substr_dict3[sub3].append(barcode)
        return [substr_dict1, substr_dict2, substr_dict3]


def get_adj_matrix_from_substring(
    barcodedict: dict[str, int], substrdictlist: list[dict[str, list[str]]]
) -> Generator[tuple[str, str], None, None]:
    """A generator that generates combinations to test for edit distance."""
    umi_length = len(list(barcodedict.keys())[0])
    if len(substrdictlist) == 2:
        substr_dict1, substr_dict2 = substrdictlist
        s = round(umi_length // 2)
        for barcode in barcodedict:
            neighbors = set()
            sub1 = barcode[:s]
            neighbors = neighbors.union(substr_dict1[sub1])
            sub2 = barcode[s:]
            neighbors = neighbors.union(substr_dict2[sub2])
            neighbors.remove(barcode)
            for neighbor in neighbors:
                yield barcode, neighbor
    if len(substrdictlist) == 3:
        substr_dict1, substr_dict2, substr_dict3 = substrdictlist
        s = round(umi_length // 3)
        for barcode in barcodedict:
            neighbors = set()
            sub1 = barcode[:s]
            neighbors = neighbors.union(substr_dict1[sub1])
            sub2 = barcode[s : 2 * s]
            neighbors = neighbors.union(substr_dict2[sub2])
            sub3 = barcode[2 * s :]
            neighbors = neighbors.union(substr_dict3[sub3])
            neighbors.remove(barcode)
            for neighbor in neighbors:
                # comb.append( (barcode, neighbor))
                yield barcode, neighbor
    # return(comb)


def cluster_barcodes(
    barcodedict: dict[str, int], edit_distance_threshold: int
) -> dict[str, list[str]]:
    """Cluster barcodes by edit distance."""
    edit_distance_threshold = int(edit_distance_threshold)
    adj_matrix: dict[str, list[str]] = {}
    if len(barcodedict) > SUBSTRING_OPTIMIZATION_THRESHOLD:
        # compare substrings for speedup
        substring_matrix = create_substring_matrix(barcodedict, edit_distance_threshold)
        comb = get_adj_matrix_from_substring(barcodedict, substring_matrix)
    else:
        comb = itertools.combinations(barcodedict.keys(), 2)
    # print(comb)
    # comb=itertools.combinations(barcodedict.keys(),2)
    for a, b in comb:
        if hamming_distance(a, b) <= edit_distance_threshold:
            if barcodedict[a] >= barcodedict[b]:
                if a not in adj_matrix:
                    adj_matrix[a] = []
                if b not in adj_matrix[a]:
                    adj_matrix[a].append(b)
            else:
                if b not in adj_matrix:
                    adj_matrix[b] = []
                if a not in adj_matrix[b]:
                    adj_matrix[b].append(a)
    return adj_matrix


def get_connected_components(
    barcodedict: dict[str, int], adj_matrix: dict[str, list[str]]
) -> list[list[str]]:
    """Get connected components from the adjacency matrix."""
    clusters: list[list[str]] = []
    added: list[str] = []
    umi_sorted = sorted(barcodedict, key=lambda x: barcodedict[x], reverse=True)  # sort umis by counts, reversed
    for umi in umi_sorted:
        if umi not in added:
            if umi in adj_matrix:
                cluster = []
                cluster.append(umi)
                added.append(umi)
                for neighbor in adj_matrix[umi]:
                    if neighbor not in added:
                        cluster.append(neighbor)
                        added.append(neighbor)
                clusters.append(cluster)
            else:
                cluster = [umi]
                clusters.append(cluster)
                added.append(umi)
    # print(len(added))
    # print(len(barcodedict))
    return clusters


def merge_clusters(
    barcodedict: dict[str, int], clusters: list[list[str]]
) -> dict[str, umi_cluster]:
    """Merge UMI clusters and return dictionary mapping barcodes to cluster info."""
    umis: dict[str, umi_cluster] = {}
    # add all umis separately
    for name, count in barcodedict.items():
        umis[name] = umi_cluster(name, count)
    for cluster in clusters:
        # merge umi counts for clusters larger than 1
        if len(cluster) > 1:
            # first item in the list is the centroid
            centroid = cluster[0]
            neighbors = cluster[1:]
            for neighbor in neighbors:
                umis[centroid].add_count(umis[neighbor].count)
            for neighbor in neighbors:
                umis[neighbor] = umis[centroid]
    # print(umis)
    return umis


def main():
    # print(hamming_distance('ATTAA','ACTAA'))
    # print(hamming_distance('ATAA','ACTAA'))
    edit_distance_threshold = 1
    with open("/home/xsteto/tmp/umierrorcorrect/test.pickle", "rb") as f:
        regions = pickle.load(f)
        # adj_matrix=cluster_barcodes(regions['2'][33141588])
        umis = regions["17"][7577495]
        adj_matrix = cluster_barcodes(regions["17"][7577495], edit_distance_threshold)
        clusters = get_connected_components(regions["17"][7577495], adj_matrix)
        # print(clusters)
        print(len(clusters))
        n = 0
        for cl in clusters:
            # if 'ACACTCGTAGTA' in cl:
            if "CATGGCGAGCCT" in cl:
                print(cl)
                for c in cl:
                    print(c, regions["17"][7577495][c])
                    print(umis[c])
            for c in cl:
                if "CATGGCGAGCCT" in c:
                    print(cl)
                n += 1
        print(n)
        # ifor a in adj_matrix:
        #    print (a+'\t'+' '.join(adj_matrix[a]))
        # print(adj_matrix)
        # umis=merge_clusters(regions['2'][33141588],adj_matrix)
        umis = merge_clusters(regions["17"][7577495], clusters)
        # print(umis)
        # for umi in umis:
        #    print(umi,umis[umi].centroid,umis[umi].count)
    with open("/home/xsteto/umierrorcorrect/umi.pickle", "wb") as g:
        pickle.dump(umis, g)


if __name__ == "__main__":
    main()
