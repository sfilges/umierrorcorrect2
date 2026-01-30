#!/usr/bin/env python3
import itertools
from collections.abc import Generator
from dataclasses import dataclass
from operator import ne

from umierrorcorrect2.core.constants import SUBSTRING_OPTIMIZATION_THRESHOLD
from umierrorcorrect2.core.logging_config import get_logger

logger = get_logger(__name__)

# Try to import numba for JIT-compiled hamming distance
try:
    from numba import njit

    _HAS_NUMBA = True
except Exception:  # Catch ImportError and other initialization errors (e.g. from llvmlite)
    _HAS_NUMBA = False


@dataclass
class umi_cluster:
    """Container for UMI cluster information."""

    centroid: str
    count: int

    def add_count(self, newcount: int) -> None:
        self.count += newcount

    def change_centroid(self, newname: str) -> None:
        self.centroid = newname


# ----------------------------------------------
# Hamming distance implementations
# ----------------------------------------------


def _hamming_original(a, b):
    """
    Returns the Hamming distance between two strings of equal length

    Original implementation from umierrorcorrect2, now deprecated.
    """
    try:
        assert len(a) == len(b)
        return sum(i != j for i, j in zip(a, b))
    except AssertionError:
        logger.warning(f"Barcode lengths are not equal for {a} and {b}.")
        raise


def _hamming_native(a: str, b: str) -> int:
    """Optimized native Python hamming distance using map."""
    return sum(map(ne, a, b))


if _HAS_NUMBA:

    @njit(cache=True)
    def _hamming_numba_core(a: bytes, b: bytes) -> int:
        """Numba JIT-compiled hamming distance for byte sequences."""
        count = 0
        for i in range(len(a)):
            if a[i] != b[i]:
                count += 1
        return count

    def _hamming_numba(a: str, b: str) -> int:
        """Numba-accelerated hamming distance wrapper."""
        return _hamming_numba_core(a.encode(), b.encode())

    hamming_distance = _hamming_numba
    _HAMMING_MSG = "Using Numba-accelerated Hamming distance."
else:
    hamming_distance = _hamming_native
    _HAMMING_MSG = "Using native Python Hamming distance."

# ----------------------------------------------
# UMI clustering functions
# ----------------------------------------------


# TODO: Do these implementations match [UMI-tools](https://github.com/CGATOxford/UMI-tools/tree/77e186c6b51917136fe5b4faf3f01ead7eb75aff)?
# TODO: Can this be improved? Or better to call UMI-tools directly?
def create_substring_matrix(barcodedict: dict[str, int], edit_distance_threshold: int) -> list[dict[str, list[str]]]:
    """Divide each barcode in two or three substrings of (approximately) equal length."""
    edit_distance_threshold = int(edit_distance_threshold)
    umi_length = len(list(barcodedict.keys())[0])
    if edit_distance_threshold <= 1:
        s = round(umi_length // 2)
        substr_dict1: dict[str, list[str]] = {}
        substr_dict2: dict[str, list[str]] = {}
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
        substr_dict1: dict[str, list[str]] = {}
        substr_dict2: dict[str, list[str]] = {}
        substr_dict3: dict[str, list[str]] = {}
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
    return []  # Default case if threshold > 2 or invalid


def get_adj_matrix_from_substring(
    barcodedict: dict[str, int], substrdictlist: list[dict[str, list[str]]]
) -> Generator[tuple[str, str], None, None]:
    """A generator that generates combinations to test for edit distance."""
    umi_length = len(list(barcodedict.keys())[0])
    if len(substrdictlist) == 2:
        substr_dict1, substr_dict2 = substrdictlist
        s = round(umi_length // 2)
        for barcode in barcodedict:
            neighbors: set[str] = set()
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
            neighbors: set[str] = set()
            sub1 = barcode[:s]
            neighbors = neighbors.union(substr_dict1[sub1])
            sub2 = barcode[s : 2 * s]
            neighbors = neighbors.union(substr_dict2[sub2])
            sub3 = barcode[2 * s :]
            neighbors = neighbors.union(substr_dict3[sub3])
            neighbors.remove(barcode)
            for neighbor in neighbors:
                yield barcode, neighbor


def cluster_barcodes(barcodedict: dict[str, int], edit_distance_threshold: int) -> dict[str, list[str]]:
    """Cluster barcodes by edit distance."""
    logger.debug(_HAMMING_MSG)
    edit_distance_threshold = int(edit_distance_threshold)
    adj_matrix: dict[str, list[str]] = {}
    if len(barcodedict) > SUBSTRING_OPTIMIZATION_THRESHOLD:
        # compare substrings for speedup
        substring_matrix = create_substring_matrix(barcodedict, edit_distance_threshold)
        comb = get_adj_matrix_from_substring(barcodedict, substring_matrix)
    else:
        comb = itertools.combinations(barcodedict.keys(), 2)
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


def get_connected_components(barcodedict: dict[str, int], adj_matrix: dict[str, list[str]]) -> list[list[str]]:
    """Get connected components from the adjacency matrix."""
    clusters: list[list[str]] = []
    added: set[str] = set()
    umi_sorted = sorted(barcodedict, key=lambda x: barcodedict[x], reverse=True)  # sort umis by counts, reversed
    for umi in umi_sorted:
        if umi not in added:
            if umi in adj_matrix:
                cluster = [umi]
                added.add(umi)
                for neighbor in adj_matrix[umi]:
                    if neighbor not in added:
                        cluster.append(neighbor)
                        added.add(neighbor)
                clusters.append(cluster)
            else:
                cluster = [umi]
                clusters.append(cluster)
                added.add(umi)
    return clusters


def merge_clusters(barcodedict: dict[str, int], clusters: list[list[str]]) -> dict[str, umi_cluster]:
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
    return umis
