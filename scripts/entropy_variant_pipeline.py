#!/usr/bin/env python3
"""
Entropy-augmented variant calling prototype.

This script supplements a standard germline pipeline by layering
entropy-aware locus discovery, lightweight local assembly, and
custom scoring logic on top of baseline VCF calls.  It produces
baseline, augmented, and novel-only VCFs plus a JSON summary of
supporting metrics.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from statistics import mean
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pysam


# --------------------------------------------------------------------------------------
# Data containers
# --------------------------------------------------------------------------------------


@dataclass
class WindowMetrics:
    chrom: str
    start: int
    end: int
    total_reads: int
    evaluated_bases: int
    mismatch_density: float
    indel_density: float
    softclip_rate: float
    base_entropy: float
    cigar_entropy: float
    alignment_entropy: float
    local_entropy_score: float
    error_rate: float
    entropy_z: float = 0.0
    error_z: float = 0.0


@dataclass
class VariantCandidate:
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
    source: str  # 'baseline' or 'assembly'
    assembly_support: float = 0.0
    contig_id: Optional[str] = None
    window: Optional[WindowMetrics] = None
    allele_balance: float = 0.0
    base_q: float = 0.0
    map_q: float = 0.0
    total_depth: int = 0
    alt_depth: int = 0
    score: float = 0.0
    classification: str = ""


class RunningStats:
    """Welford running statistics."""

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0

    def update(self, value: float) -> None:
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        delta2 = value - self.mean
        self.m2 += delta * delta2

    @property
    def variance(self) -> float:
        if self.n < 2:
            return 0.0
        return self.m2 / (self.n - 1)

    @property
    def std(self) -> float:
        return math.sqrt(self.variance)


class IntervalIndex:
    """Simple interval lookups without external deps."""

    def __init__(self) -> None:
        self._intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        self._starts: Dict[str, List[int]] = {}

    def add(self, chrom: str, start: int, end: int) -> None:
        self._intervals[chrom].append((start, end))

    def finalize(self) -> None:
        for chrom, entries in self._intervals.items():
            entries.sort()
            self._starts[chrom] = [s for s, _ in entries]

    def contains(self, chrom: str, pos0: int) -> bool:
        """pos0 is 0-based."""
        entries = self._intervals.get(chrom)
        if not entries:
            return False
        starts = self._starts[chrom]
        import bisect

        idx = bisect.bisect_right(starts, pos0) - 1
        if idx < 0:
            return False
        start, end = entries[idx]
        return start <= pos0 < end


# --------------------------------------------------------------------------------------
# Helper utilities
# --------------------------------------------------------------------------------------


CIGAR_OPS = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
}


def shannon_entropy(counts: Iterable[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    entropy = 0.0
    for value in counts:
        if value == 0:
            continue
        p = value / total
        entropy -= p * math.log2(p)
    return entropy


def ti_tv_ratio(variants: Sequence[VariantCandidate]) -> float:
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    ti = 0
    tv = 0
    for variant in variants:
        if len(variant.ref) == len(variant.alt) == 1:
            pair = (variant.ref.upper(), variant.alt.upper())
            if pair[0] == pair[1]:
                continue
            if pair in transitions:
                ti += 1
            else:
                tv += 1
    if tv == 0:
        return float("inf") if ti > 0 else 0.0
    return ti / tv


def indel_length_distribution(variants: Sequence[VariantCandidate]) -> Dict[int, int]:
    distribution: Dict[int, int] = defaultdict(int)
    for variant in variants:
        if len(variant.ref) != len(variant.alt):
            delta = abs(len(variant.ref) - len(variant.alt))
            distribution[delta] += 1
    return dict(sorted(distribution.items()))


def variant_type(variant: VariantCandidate) -> str:
    if len(variant.ref) == len(variant.alt):
        if len(variant.ref) == 1:
            return "SNV"
        return "MNV"
    if len(variant.ref) < len(variant.alt):
        if len(variant.alt) - len(variant.ref) >= 20:
            return "micro_insertion"
        return "insertion"
    if len(variant.ref) - len(variant.alt) >= 20:
        return "micro_deletion"
    return "deletion"


def clamp(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


# --------------------------------------------------------------------------------------
# Complexity scanner
# --------------------------------------------------------------------------------------


class ComplexityScanner:
    def __init__(
        self,
        bam_path: str,
        reference_path: str,
        window_size: int,
        step: int,
        min_reads: int,
        min_mapq: int,
    ) -> None:
        self.bam_path = bam_path
        self.reference_path = reference_path
        self.window_size = window_size
        self.step = step
        self.min_reads = min_reads
        self.min_mapq = min_mapq
        self._contig_lengths: Dict[str, int] = {}
        self._window_cache: Dict[Tuple[str, int], WindowMetrics] = {}
        self.covered_bases = 0

        with pysam.AlignmentFile(
            self.bam_path, "rb", reference_filename=self.reference_path
        ) as bam:
            for contig in bam.header.references:
                self._contig_lengths[contig] = bam.get_reference_length(contig)

    # -- window iteration ----------------------------------------------------------------
    def _window_iter(self) -> Iterator[WindowMetrics]:
        with pysam.AlignmentFile(
            self.bam_path, "rb", reference_filename=self.reference_path
        ) as bam:
            for chrom, length in self._contig_lengths.items():
                for start in range(0, length, self.step):
                    end = min(start + self.window_size, length)
                    metrics = self._compute_metrics(bam, chrom, start, end)
                    if metrics.total_reads < self.min_reads:
                        continue
                    yield metrics

    def _compute_metrics(
        self, bam: pysam.AlignmentFile, chrom: str, start: int, end: int
    ) -> WindowMetrics:
        mismatch = 0.0
        indel = 0.0
        softclip = 0.0
        total_bases = 0
        read_count = 0
        base_counts: Counter[str] = Counter()
        cigar_counts: Counter[str] = Counter()

        window_len = max(end - start, 1)
        for read in bam.fetch(contig=chrom, start=start, stop=end):
            if not self._usable_read(read):
                continue
            overlap = read.get_overlap(start, end)
            if overlap <= 0:
                continue
            read_count += 1
            total_bases += overlap
            aln_len = max(read.query_alignment_length or 1, 1)
            frac = overlap / aln_len
            nm = read.get_tag("NM") if read.has_tag("NM") else 0
            indel_bases = sum(
                length for op, length in (read.cigartuples or []) if op in (1, 2)
            )
            softclip_bases = sum(
                length for op, length in (read.cigartuples or []) if op == 4
            )
            mismatch += max(nm - indel_bases, 0) * frac
            indel += indel_bases * frac
            softclip += softclip_bases * frac
            for op, length in (read.cigartuples or []):
                cigar_counts[CIGAR_OPS.get(op, "?")] += length
            self._update_base_counts(read, base_counts, start, end)

        mismatch_density = mismatch / max(total_bases, 1)
        indel_density = indel / max(total_bases, 1)
        softclip_rate = softclip / max(total_bases, 1)
        base_entropy = shannon_entropy(base_counts.values())
        cigar_entropy = shannon_entropy(cigar_counts.values())
        # matches approximate
        match_prob = clamp(1.0 - mismatch_density - indel_density, 0.0, 1.0)
        alignment_entropy = shannon_entropy(
            [
                int(match_prob * 1000),
                int(mismatch_density * 1000),
                int(indel_density * 1000),
            ]
        )
        local_entropy_score = (
            0.5 * base_entropy + 0.3 * cigar_entropy + 0.2 * alignment_entropy
        )
        error_rate = mismatch_density + indel_density

        metrics = WindowMetrics(
            chrom=chrom,
            start=start,
            end=end,
            total_reads=read_count,
            evaluated_bases=window_len,
            mismatch_density=mismatch_density,
            indel_density=indel_density,
            softclip_rate=softclip_rate,
            base_entropy=base_entropy,
            cigar_entropy=cigar_entropy,
            alignment_entropy=alignment_entropy,
            local_entropy_score=local_entropy_score,
            error_rate=error_rate,
        )
        return metrics

    def _update_base_counts(
        self,
        read: pysam.AlignedSegment,
        base_counts: Counter[str],
        start: int,
        end: int,
    ) -> None:
        if read.query_sequence is None:
            return
        seq = read.query_sequence
        if not seq:
            return
        quals = read.query_qualities
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if ref_pos is None:
                continue
            if start <= ref_pos < end and query_pos is not None:
                base = seq[query_pos]
                if base == "N":
                    continue
                base_counts[base.upper()] += 1

    def _usable_read(self, read: pysam.AlignedSegment) -> bool:
        if read.is_unmapped:
            return False
        if read.is_secondary or read.is_supplementary:
            return False
        if read.is_duplicate:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        return True

    # -- public API ---------------------------------------------------------------------

    def detect_complex_windows(
        self,
        entropy_z: float,
        error_z: float,
    ) -> Tuple[List[WindowMetrics], Dict[str, float]]:
        entropy_stats = RunningStats()
        error_stats = RunningStats()
        covered_bases = 0
        first_pass: List[WindowMetrics] = []
        for metrics in self._window_iter():
            entropy_stats.update(metrics.local_entropy_score)
            error_stats.update(metrics.error_rate)
            covered_bases += metrics.evaluated_bases
            first_pass.append(metrics)

        self.covered_bases = covered_bases
        entropy_threshold = (
            entropy_stats.mean + entropy_z * entropy_stats.std
            if entropy_stats.n
            else 0.0
        )
        error_threshold = (
            error_stats.mean + error_z * error_stats.std if error_stats.n else 0.0
        )

        flagged: List[WindowMetrics] = []
        for metrics in first_pass:
            metrics.entropy_z = (
                (metrics.local_entropy_score - entropy_stats.mean)
                / entropy_stats.std
                if entropy_stats.std > 0
                else 0.0
            )
            metrics.error_z = (
                (metrics.error_rate - error_stats.mean) / error_stats.std
                if error_stats.std > 0
                else 0.0
            )
            is_complex = (
                metrics.local_entropy_score >= entropy_threshold
                or metrics.error_rate >= error_threshold
            )
            if is_complex:
                flagged.append(metrics)
            cache_key = (metrics.chrom, metrics.start)
            if cache_key not in self._window_cache:
                self._window_cache[cache_key] = metrics

        stats = {
            "entropy_mean": entropy_stats.mean,
            "entropy_std": entropy_stats.std,
            "entropy_threshold": entropy_threshold,
            "error_mean": error_stats.mean,
            "error_std": error_stats.std,
            "error_threshold": error_threshold,
            "covered_bases": covered_bases,
        }
        return flagged, stats

    def get_window_for_position(self, chrom: str, pos1: int) -> WindowMetrics:
        start = max(((pos1 - 1) // self.step) * self.step, 0)
        cache_key = (chrom, start)
        if cache_key in self._window_cache:
            return self._window_cache[cache_key]

        with pysam.AlignmentFile(
            self.bam_path, "rb", reference_filename=self.reference_path
        ) as bam:
            end = min(start + self.window_size, self._contig_lengths[chrom])
            metrics = self._compute_metrics(bam, chrom, start, end)
            self._window_cache[cache_key] = metrics
            return metrics


# --------------------------------------------------------------------------------------
# Local assembler (naive De Bruijn graph)
# --------------------------------------------------------------------------------------


class DeBruijnGraph:
    def __init__(self, k: int) -> None:
        self.k = k
        self.graph: Dict[str, Counter[str]] = defaultdict(Counter)
        self.indegree: Counter[str] = Counter()

    def add_sequence(self, seq: str) -> None:
        seq = seq.upper()
        if len(seq) < self.k:
            return
        for i in range(len(seq) - self.k + 1):
            kmer = seq[i : i + self.k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix][suffix] += 1
            self.indegree[suffix] += 1

    def assemble(self, max_contigs: int = 25) -> List[Tuple[str, float]]:
        contigs: List[Tuple[str, float]] = []
        visited_edges = set()
        nodes = set(self.graph.keys()) | set(self.indegree.keys())
        for node in nodes:
            out_neighbors = self.graph.get(node)
            indeg = self.indegree.get(node, 0)
            outdeg = sum(out_neighbors.values()) if out_neighbors else 0
            if outdeg == 0:
                continue
            if indeg != 1 or (out_neighbors and len(out_neighbors) != 1):
                contig = self._walk(node, visited_edges)
                if contig:
                    contigs.append(contig)
            if len(contigs) >= max_contigs:
                break

        # handle isolated cycles
        if len(contigs) < max_contigs:
            for node in nodes:
                if all(((node, nxt) in visited_edges) for nxt in self.graph.get(node, [])):
                    continue
                contig = self._walk(node, visited_edges)
                if contig:
                    contigs.append(contig)
                if len(contigs) >= max_contigs:
                    break
        return contigs

    def _walk(
        self, start: str, visited_edges: set
    ) -> Optional[Tuple[str, float]]:
        current = start
        sequence = current
        support = float("inf")
        while True:
            neighbors = self.graph.get(current)
            if not neighbors:
                break
            next_node, weight = max(neighbors.items(), key=lambda item: item[1])
            edge = (current, next_node)
            if edge in visited_edges:
                break
            visited_edges.add(edge)
            support = min(support, float(weight))
            sequence += next_node[-1]
            indeg = self.indegree.get(next_node, 0)
            out_neighbors = self.graph.get(next_node)
            if not out_neighbors or indeg != 1 or len(out_neighbors) != 1:
                current = next_node
                continue
            current = next_node
        if len(sequence) <= self.k:
            return None
        return sequence, support if support != float("inf") else 0.0


# --------------------------------------------------------------------------------------
# Variant pipeline
# --------------------------------------------------------------------------------------


class EntropyVariantPipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.sample_ids: List[str] = []
        self.baseline_variants: Dict[Tuple[str, int, str, str], VariantCandidate] = {}
        self.assembly_variants: Dict[Tuple[str, int, str, str], VariantCandidate] = {}
        self.flagged_windows: List[WindowMetrics] = []
        self.window_index = IntervalIndex()
        self.stats: Dict[str, float] = {}
        self.scanner = ComplexityScanner(
            bam_path=args.bam,
            reference_path=args.reference,
            window_size=args.window_size,
            step=args.window_step,
            min_reads=args.min_window_reads,
            min_mapq=args.min_mapq,
        )
        self.reference = pysam.FastaFile(args.reference)

    # -- pipeline -----------------------------------------------------------------------
    def run(self) -> None:
        self._load_baseline_variants()
        self.flagged_windows, self.stats = self.scanner.detect_complex_windows(
            entropy_z=self.args.entropy_z, error_z=self.args.error_z
        )
        for window in self.flagged_windows:
            self.window_index.add(window.chrom, window.start, window.end)
        self.window_index.finalize()
        self._discover_variants_from_windows()
        self._score_variants()
        self._write_outputs()
        self._write_summary()

    # -- baseline -----------------------------------------------------------------------
    def _load_baseline_variants(self) -> None:
        vcf = pysam.VariantFile(self.args.baseline_vcf)
        self.sample_ids = list(vcf.header.samples)
        for record in vcf:
            for alt in record.alts or []:
                key = (record.chrom, record.pos, record.ref, alt)
                candidate = VariantCandidate(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alt=alt,
                    source="baseline",
                    classification="standard",
                )
                self.baseline_variants[key] = candidate
        vcf.close()

    # -- candidate discovery ------------------------------------------------------------
    def _discover_variants_from_windows(self) -> None:
        for window in self.flagged_windows:
            reads = self._extract_read_sequences(window)
            if not reads:
                continue
            graph = DeBruijnGraph(k=self.args.kmer)
            for seq in reads:
                graph.add_sequence(seq)
            contigs = graph.assemble(max_contigs=self.args.max_contigs)
            ref_seq = self.reference.fetch(window.chrom, window.start, window.end)
            for idx, (contig, support) in enumerate(contigs):
                variants = self._derive_variants_from_contig(
                    chrom=window.chrom,
                    window_start=window.start,
                    ref_seq=ref_seq,
                    contig_seq=contig,
                )
                for variant in variants:
                    key = (variant.chrom, variant.pos, variant.ref, variant.alt)
                    if key in self.baseline_variants:
                        continue
                    variant.source = "assembly"
                    variant.assembly_support = support
                    variant.contig_id = f"{window.chrom}:{window.start}-{window.end}:contig{idx}"
                    variant.window = window
                    existing = self.assembly_variants.get(key)
                    if existing:
                        existing.assembly_support = max(
                            existing.assembly_support, support
                        )
                    else:
                        self.assembly_variants[key] = variant

    def _extract_read_sequences(self, window: WindowMetrics) -> List[str]:
        sequences: List[str] = []
        with pysam.AlignmentFile(
            self.args.bam, "rb", reference_filename=self.args.reference
        ) as bam:
            for read in bam.fetch(window.chrom, window.start, window.end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < self.args.min_mapq:
                    continue
                if read.query_sequence:
                    sequences.append(read.query_sequence.upper())
        return sequences

    def _derive_variants_from_contig(
        self,
        chrom: str,
        window_start: int,
        ref_seq: str,
        contig_seq: str,
    ) -> List[VariantCandidate]:
        matcher = SequenceMatcher(None, ref_seq, contig_seq)
        variants: List[VariantCandidate] = []
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "equal":
                continue
            global_pos = window_start + i1
            if tag == "replace":
                ref = ref_seq[i1:i2]
                alt = contig_seq[j1:j2]
                if not ref or not alt:
                    continue
                variant = VariantCandidate(
                    chrom=chrom,
                    pos=global_pos + 1,
                    ref=ref,
                    alt=alt,
                    source="assembly",
                )
                variants.append(variant)
            elif tag == "delete":
                deleted = ref_seq[i1:i2]
                if i1 == 0:
                    continue
                anchor_base = ref_seq[i1 - 1]
                ref = anchor_base + deleted
                alt = anchor_base
                pos = window_start + i1
                variant = VariantCandidate(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    source="assembly",
                )
                variants.append(variant)
            elif tag == "insert":
                inserted = contig_seq[j1:j2]
                if i1 == 0:
                    continue
                anchor_base = ref_seq[i1 - 1]
                ref = anchor_base
                alt = anchor_base + inserted
                pos = window_start + i1
                variant = VariantCandidate(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    source="assembly",
                )
                variants.append(variant)
        return variants

    # -- scoring ------------------------------------------------------------------------
    def _score_variants(self) -> None:
        all_variants = list(self.baseline_variants.values()) + list(
            self.assembly_variants.values()
        )
        with pysam.AlignmentFile(
            self.args.bam, "rb", reference_filename=self.args.reference
        ) as bam:
            for variant in all_variants:
                support = self._compute_read_support(bam, variant)
                variant.alt_depth = support["alt_depth"]
                variant.total_depth = support["total_depth"]
                variant.allele_balance = support["allele_balance"]
                variant.base_q = support["mean_baseq"]
                variant.map_q = support["mean_mapq"]
                variant.window = variant.window or self.scanner.get_window_for_position(
                    variant.chrom, variant.pos
                )
                entropy_component = variant.window.entropy_z if variant.window else 0.0
                entropy_component = clamp(
                    entropy_component / max(self.args.entropy_z, 1e-3), 0.0, 1.5
                )
                baseq_component = clamp(variant.base_q / 40.0, 0.0, 1.0)
                mapq_component = clamp(variant.map_q / 60.0, 0.0, 1.0)
                variant.score = (
                    self.args.weight_allele_balance * variant.allele_balance
                    + self.args.weight_baseq * baseq_component
                    + self.args.weight_mapq * mapq_component
                    + self.args.weight_entropy * entropy_component
                    + self.args.weight_assembly
                    * clamp(variant.assembly_support / 20.0, 0.0, 1.0)
                )
                if variant.source == "baseline":
                    variant.classification = "standard"
                elif variant.score >= self.args.novel_score_threshold:
                    variant.classification = "entropy_supported"
                elif variant.score >= self.args.low_conf_threshold:
                    variant.classification = "low_confidence"
                else:
                    variant.classification = "low_confidence"

    def _compute_read_support(
        self, bam: pysam.AlignmentFile, variant: VariantCandidate
    ) -> Dict[str, float]:
        start0 = variant.pos - 1
        span = max(len(variant.ref), len(variant.alt))
        region_start = max(0, start0 - 2)
        region_end = start0 + span + 2
        alt_reads = 0
        ref_reads = 0
        base_qualities: List[float] = []
        map_qualities: List[float] = []
        for read in bam.fetch(variant.chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < self.args.min_mapq:
                continue
            supports_alt, quals = self._read_supports_variant(read, variant)
            if supports_alt:
                alt_reads += 1
                if read.mapping_quality > 0:
                    map_qualities.append(read.mapping_quality)
                base_qualities.extend(quals)
            else:
                ref_reads += 1
        total_reads = alt_reads + ref_reads
        allele_balance = alt_reads / total_reads if total_reads else 0.0
        mean_baseq = mean(base_qualities) if base_qualities else 0.0
        mean_mapq = mean(map_qualities) if map_qualities else 0.0
        return {
            "alt_depth": alt_reads,
            "total_depth": total_reads,
            "allele_balance": allele_balance,
            "mean_baseq": mean_baseq,
            "mean_mapq": mean_mapq,
        }

    def _read_supports_variant(
        self, read: pysam.AlignedSegment, variant: VariantCandidate
    ) -> Tuple[bool, List[float]]:
        if len(variant.ref) == len(variant.alt):
            return self._read_supports_substitution(read, variant)
        if len(variant.ref) < len(variant.alt):
            return self._read_supports_insertion(read, variant)
        return self._read_supports_deletion(read, variant)

    def _read_supports_substitution(
        self, read: pysam.AlignedSegment, variant: VariantCandidate
    ) -> Tuple[bool, List[float]]:
        if read.query_sequence is None:
            return False, []
        seq = read.query_sequence
        quals = read.query_qualities
        pos0 = variant.pos - 1
        needed = len(variant.ref)
        bases = [""] * needed
        qvals = [0.0] * needed
        count = 0
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if ref_pos is None or query_pos is None:
                continue
            if pos0 <= ref_pos < pos0 + needed:
                idx = ref_pos - pos0
                bases[idx] = seq[query_pos]
                if quals:
                    qvals[idx] = quals[query_pos]
                count += 1
        if count != needed:
            return False, []
        if "".join(bases).upper() == variant.alt.upper():
            return True, qvals
        return False, []

    def _read_supports_insertion(
        self, read: pysam.AlignedSegment, variant: VariantCandidate
    ) -> Tuple[bool, List[float]]:
        seq = read.query_sequence
        if seq is None:
            return False, []
        quals = read.query_qualities
        anchor_pos = variant.pos - 1
        insertion = variant.alt[1:]
        found_anchor = False
        inserted: List[str] = []
        inserted_quals: List[float] = []
        for query_pos, ref_pos in read.get_aligned_pairs(with_seq=True):
            if ref_pos == anchor_pos and query_pos is not None:
                found_anchor = True
                continue
            if not found_anchor:
                continue
            if ref_pos is None and query_pos is not None:
                inserted.append(seq[query_pos])
                if quals:
                    inserted_quals.append(quals[query_pos])
            else:
                break
        if inserted and "".join(inserted).upper() == insertion.upper():
            return True, inserted_quals
        return False, []

    def _read_supports_deletion(
        self, read: pysam.AlignedSegment, variant: VariantCandidate
    ) -> Tuple[bool, List[float]]:
        seq = read.query_sequence
        if seq is None:
            return False, []
        quals = read.query_qualities
        anchor_pos = variant.pos - 1
        deleted_len = len(variant.ref) - 1
        found_anchor = False
        deleted = 0
        anchor_qual = 0.0
        for query_pos, ref_pos in read.get_aligned_pairs(with_seq=True):
            if ref_pos == anchor_pos and query_pos is not None:
                found_anchor = True
                if quals:
                    anchor_qual = quals[query_pos]
                continue
            if not found_anchor:
                continue
            if ref_pos is not None and query_pos is None:
                deleted += 1
            else:
                break
        if deleted == deleted_len and deleted_len > 0:
            return True, [anchor_qual]
        return False, []

    # -- outputs ------------------------------------------------------------------------
    def _write_outputs(self) -> None:
        prefix = self.args.output_prefix
        os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
        baseline_out = f"{prefix}.baseline.vcf.gz"
        shutil.copyfile(self.args.baseline_vcf, baseline_out)
        baseline_index = f"{baseline_out}.tbi"
        if os.path.exists(self.args.baseline_vcf_index):
            shutil.copyfile(self.args.baseline_vcf_index, baseline_index)

        header = pysam.VariantFile(self.args.baseline_vcf).header.copy()
        info_fields = {
            "COMPLEX_SCORE": "##INFO=<ID=COMPLEX_SCORE,Number=1,Type=Float,Description=\"Composite entropy score\">",
            "LOCAL_ENTROPY_Z": "##INFO=<ID=LOCAL_ENTROPY_Z,Number=1,Type=Float,Description=\"Window entropy z-score\">",
            "ALT_DEPTH": "##INFO=<ID=ALT_DEPTH,Number=1,Type=Integer,Description=\"Alt-supporting reads\">",
            "TOTAL_DEPTH": "##INFO=<ID=TOTAL_DEPTH,Number=1,Type=Integer,Description=\"Total usable reads\">",
            "ASSEMBLY_SUPPORT": "##INFO=<ID=ASSEMBLY_SUPPORT,Number=1,Type=Float,Description=\"Assembly contig support weight\">",
            "CLASSIFICATION": "##INFO=<ID=CLASSIFICATION,Number=1,Type=String,Description=\"Variant classification\">",
        }
        for info_id, line in info_fields.items():
            if info_id not in header.info:
                header.add_line(line)

        format_fields = {
            "AD": "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">",
            "DP": "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
            "GT": "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        }
        for fmt_id, line in format_fields.items():
            if fmt_id not in header.formats:
                header.add_line(line)

        augmented_path = f"{prefix}.augmented_variants.vcf.gz"
        novel_path = f"{prefix}.novel_entropy_variants.vcf.gz"
        low_conf_path = f"{prefix}.low_confidence_variants.vcf.gz"
        augmented_vcf = pysam.VariantFile(augmented_path, "wz", header=header)
        novel_vcf = pysam.VariantFile(novel_path, "wz", header=header)
        low_conf_vcf = pysam.VariantFile(low_conf_path, "wz", header=header)

        ordered_baseline = self._iter_baseline_records()
        novel_variants = [
            v for v in self.assembly_variants.values() if v.classification == "entropy_supported"
        ]
        low_conf_variants = [
            v for v in self.assembly_variants.values() if v.classification != "entropy_supported"
        ]

        for record in ordered_baseline:
            for alt in record.alts or []:
                key = (record.chrom, record.pos, record.ref, alt)
                candidate = self.baseline_variants.get(key)
                if not candidate:
                    continue
                self._annotate_record(record, candidate)
            augmented_vcf.write(record)

        for variant in novel_variants:
            record = self._build_record(header, variant)
            novel_vcf.write(record)
            augmented_vcf.write(record)

        for variant in low_conf_variants:
            record = self._build_record(header, variant)
            low_conf_vcf.write(record)

        augmented_vcf.close()
        novel_vcf.close()
        low_conf_vcf.close()
        pysam.tabix_index(augmented_path, preset="vcf", force=True)
        pysam.tabix_index(novel_path, preset="vcf", force=True)
        pysam.tabix_index(low_conf_path, preset="vcf", force=True)

        self.outputs = {
            "baseline_vcf": baseline_out,
            "augmented_vcf": augmented_path,
            "novel_vcf": novel_path,
            "low_conf_vcf": low_conf_path,
        }

    def _iter_baseline_records(self) -> Iterator[pysam.VariantRecord]:
        vcf = pysam.VariantFile(self.args.baseline_vcf)
        for record in vcf:
            yield record
        vcf.close()

    def _annotate_record(
        self, record: pysam.VariantRecord, variant: VariantCandidate
    ) -> None:
        record.info["COMPLEX_SCORE"] = round(variant.score, 3)
        record.info["LOCAL_ENTROPY_Z"] = (
            round(variant.window.entropy_z, 3) if variant.window else 0.0
        )
        record.info["ALT_DEPTH"] = variant.alt_depth
        record.info["TOTAL_DEPTH"] = variant.total_depth
        record.info["ASSEMBLY_SUPPORT"] = round(variant.assembly_support, 3)
        record.info["CLASSIFICATION"] = variant.classification

    def _build_record(
        self, header: pysam.VariantHeader, variant: VariantCandidate
    ) -> pysam.VariantRecord:
        record = header.new_record(
            contig=variant.chrom,
            start=variant.pos - 1,
            stop=variant.pos - 1 + len(variant.ref),
            alleles=(variant.ref, variant.alt),
        )
        record.info["COMPLEX_SCORE"] = round(variant.score, 3)
        record.info["LOCAL_ENTROPY_Z"] = (
            round(variant.window.entropy_z, 3) if variant.window else 0.0
        )
        record.info["ALT_DEPTH"] = variant.alt_depth
        record.info["TOTAL_DEPTH"] = variant.total_depth
        record.info["ASSEMBLY_SUPPORT"] = round(variant.assembly_support, 3)
        record.info["CLASSIFICATION"] = variant.classification
        genotype = (None, None)
        if variant.allele_balance >= 0.8:
            genotype = (1, 1)
        elif variant.allele_balance >= 0.3:
            genotype = (0, 1)
        else:
            genotype = (0, 0)
        ref_depth = max(variant.total_depth - variant.alt_depth, 0)
        for sample in self.sample_ids:
            record.samples[sample]["GT"] = genotype
            record.samples[sample]["DP"] = variant.total_depth
            record.samples[sample]["AD"] = (ref_depth, variant.alt_depth)
        return record

    # -- summary ------------------------------------------------------------------------
    def _write_summary(self) -> None:
        prefix = self.args.output_prefix
        summary_path = f"{prefix}.entropy_summary.json"
        baseline = list(self.baseline_variants.values())
        novel = [
            v for v in self.assembly_variants.values() if v.classification == "entropy_supported"
        ]
        low_conf = [
            v for v in self.assembly_variants.values() if v.classification != "entropy_supported"
        ]
        augmented = baseline + novel

        flagged_bp = sum(w.end - w.start for w in self.flagged_windows)
        non_flagged_bp = max(self.scanner.covered_bases - flagged_bp, 1)
        flagged_variants = [
            v
            for v in augmented
            if self.window_index.contains(v.chrom, v.pos - 1)
        ]
        normal_variants = [v for v in augmented if v not in flagged_variants]

        summary = {
            "counts": {
                "baseline": len(baseline),
                "entropy_supported": len(novel),
                "low_confidence": len(low_conf),
            },
            "ti_tv": {
                "baseline": ti_tv_ratio(baseline),
                "augmented": ti_tv_ratio(augmented),
                "novel": ti_tv_ratio(novel),
            },
            "indel_lengths": {
                "novel": indel_length_distribution(novel),
                "low_confidence": indel_length_distribution(low_conf),
            },
            "variant_density": {
                "flagged_regions_per_mb": len(flagged_variants)
                / max(flagged_bp / 1_000_000, 1e-6),
                "normal_regions_per_mb": len(normal_variants)
                / max(non_flagged_bp / 1_000_000, 1e-6),
            },
            "entropy_thresholds": self.stats,
            "novel_overlap_with_baseline": len(
                set(self.baseline_variants.keys())
                & set(self.assembly_variants.keys())
            ),
            "novel_fraction_in_flagged": (
                sum(
                    1
                    for v in novel
                    if self.window_index.contains(v.chrom, v.pos - 1)
                )
                / len(novel)
                if novel
                else 0.0
            ),
            "novel_read_support": {
                "mean_alt_depth": mean([v.alt_depth for v in novel]) if novel else 0.0,
                "mean_allele_balance": mean([v.allele_balance for v in novel])
                if novel
                else 0.0,
                "mean_mapq": mean([v.map_q for v in novel]) if novel else 0.0,
                "mean_baseq": mean([v.base_q for v in novel]) if novel else 0.0,
            },
        }
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)
        self.outputs["summary"] = summary_path


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Entropy-aware variant augmentation pipeline"
    )
    parser.add_argument("--bam", required=True, help="Input BAM/CRAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument(
        "--baseline-vcf", required=True, help="Baseline VCF from standard pipeline"
    )
    parser.add_argument(
        "--baseline-vcf-index",
        required=True,
        help="Index (tbi) for the baseline VCF",
    )
    parser.add_argument(
        "--output-prefix", required=True, help="Prefix for augmented outputs"
    )
    parser.add_argument("--window-size", type=int, default=120)
    parser.add_argument("--window-step", type=int, default=60)
    parser.add_argument("--min-window-reads", type=int, default=5)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--kmer", type=int, default=25)
    parser.add_argument("--max-contigs", type=int, default=30)
    parser.add_argument("--entropy-z", type=float, default=0.75)
    parser.add_argument("--error-z", type=float, default=1.0)
    parser.add_argument("--novel-score-threshold", type=float, default=0.7)
    parser.add_argument("--low-conf-threshold", type=float, default=0.4)
    parser.add_argument("--weight-allele-balance", type=float, default=0.4)
    parser.add_argument("--weight-baseq", type=float, default=0.15)
    parser.add_argument("--weight-mapq", type=float, default=0.15)
    parser.add_argument("--weight-entropy", type=float, default=0.2)
    parser.add_argument("--weight-assembly", type=float, default=0.1)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    pipeline = EntropyVariantPipeline(args)
    pipeline.run()


if __name__ == "__main__":
    main()
