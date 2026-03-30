"""
RBP Motif Enrichment Analysis following rMAPS method.
Reference: Park et al. (2016) NAR, Hwang et al. (2020) NAR
"""

import os
import re
from typing import Optional, List, Dict, Tuple
from collections import defaultdict
import numpy as np
from scipy import stats

from .motif_analysis import get_ensembl_db


SEQUENCE_FLANK_SIZE = 250
EXON_EDGE_SIZE = 50
SLIDING_WINDOW_SIZE = 50
SPLICE_SITE_3P_MASK = 20
SPLICE_SITE_5P_MASK = 6


RBP_MOTIFS = {
    "PTB": ["UCUCUCU", "UUCUUC"],
    "RBM25": ["ACUGU"],
    "QKI": ["ACUAUAC", "ACUGU"],
    "NOVA1": ["UCAU", "YCAY"],
    "MBNL1": ["GCUGC", "UGCU"],
    "CELF1": ["GUGU", "UGUU"],
    "CELF2": ["GUGU", "UGUU"],
    "TIA1": ["UUAUU", "UURU"],
    "TIAR": ["UUAUU", "UURU"],
    "HuR": ["UUU", "UUAU"],
    "TDP43": ["UGUGU", "GUUG"],
    "FUS": ["GGUGG"],
    "SF3B4": ["CCU", "UCC"],
    "U2AF65": ["AGG", "UAGG"],
    "HNRNPC": ["UAG", "CUNK"],
    "hnRNPA1": ["UAGGG", "AGG"],
    "RPS14": ["ACU", "GACU"],
    "ESRP1": ["GSCG", "GCCG"],
    "ESRP2": ["GSCG", "GCCG"],
    "RBFOX1": ["GCAUG"],
    "RBFOX2": ["GCAUG"],
    "RBFOX3": ["GCAUG"],
    "MBNL2": ["GCUGC", "UGCU"],
    "STAR": ["ACUAA"],
    "SRSF1": ["AGGAC"],
    "SRSF2": ["SCG"],
    "SRSF3": ["GGC"],
    "TRA2B": ["GAAGA"],
    "hnRNPF": ["GGG", "GGGG"],
    "hnRNPH": ["GGGR"],
    "U2AF1": ["YAG"],
    "KHSRP": ["YUGNNY"],
    "PABPC1": ["AAUAAA"],
    "ELAVL1": ["UUU", "UUUU"],
    "ELAVL2": ["UUU", "UUUU"],
    "ELAVL3": ["UUU", "UUUU"],
    "ELAVL4": ["UUU", "UUUU"],
    "IGF2BP1": ["CAUH"],
    "IGF2BP2": ["CAUH"],
    "IGF2BP3": ["CAUH"],
    "QKI_para": ["ACUAUAA"],
    "RBM47": ["UAG"],
    "ZC3H10": ["UUAUU"],
    "RBPMS": ["UUGURY"],
    "PCBP1": ["CCCC"],
    "PCBP2": ["CCCC"],
    "NCOA7": ["UAAU"],
    "LIN28A": ["CCNNGG"],
}


def calculate_motif_density(sequence: str, motif: str) -> float:
    """Calculate motif density (% nucleotides covered by motif occurrences)."""
    if not sequence or not motif:
        return 0.0

    sequence = sequence.upper()
    motif = motif.upper()

    count = 0
    for i in range(len(sequence) - len(motif) + 1):
        match = True
        for j, m in enumerate(motif):
            if m not in ["A", "C", "G", "U"]:
                continue
            if sequence[i + j] != m:
                match = False
                break
        if match:
            count += len(motif)

    return count / len(sequence) * 100


def sliding_window_density(
    sequence: str, motif: str, window_size: int = SLIDING_WINDOW_SIZE
) -> List[float]:
    """
    Calculate motif density in sliding windows.
    Returns list of densities for each window position.
    """
    if not sequence or len(sequence) < window_size:
        return [calculate_motif_density(sequence, motif)]

    densities = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        densities.append(calculate_motif_density(window, motif))

    return densities


def extract_genomic_region(
    chromosome: str,
    start: int,
    end: int,
    strand: str,
    flank_size: int = 0,
    ensembl_db=None,
) -> Optional[str]:
    """Extract genomic DNA sequence for a region using pyensembl."""
    if ensembl_db is None:
        ensembl_db = get_ensembl_db(109, "human")

    if ensembl_db is None:
        return None

    try:
        chrom = chromosome.replace("chr", "")
        if chrom in ["X", "Y", "M", "MT"]:
            chrom = chrom if chrom != "M" else "MT"
        else:
            try:
                chrom = str(int(chrom))
            except ValueError:
                pass

        region_start = max(0, start - flank_size)
        region_end = end + flank_size

        sequence = ensembl_db.region(chrom, region_start, region_end)

        if sequence:
            sequence = sequence.upper()
            if strand == "-":
                sequence = reverse_complement(sequence)

        return sequence

    except Exception:
        return None


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA/RNA sequence."""
    complement = {
        "A": "U",
        "T": "A",
        "G": "C",
        "C": "G",
        "U": "A",
        "a": "u",
        "t": "a",
        "g": "c",
        "c": "g",
        "u": "a",
    }
    return "".join(complement.get(base, base) for base in reversed(seq))


def mask_splice_sites(
    sequence: str, is_upstream: bool = True, strand: str = "+"
) -> str:
    """
    Mask splice site regions that are highly constrained.
    - 3' splice site: exclude 20bp
    - 5' splice site: exclude 6bp
    """
    if is_upstream:
        return sequence[SPLICE_SITE_3P_MASK:]
    else:
        return (
            sequence[:-SPLICE_SITE_5P_MASK]
            if len(sequence) > SPLICE_SITE_5P_MASK
            else sequence
        )


class SplicingEventRegions:
    """Container for extracted regions of a splicing event."""

    def __init__(self, event_id: str, event_type: str):
        self.event_id = event_id
        self.event_type = event_type
        self.upstream_exon = ""
        self.upstream_intron = ""
        self.target_exon = ""
        self.downstream_intron = ""
        self.downstream_exon = ""

    def all_sequences(self) -> List[str]:
        """Get all non-empty sequences."""
        return [
            s
            for s in [
                self.upstream_exon,
                self.upstream_intron,
                self.target_exon,
                self.downstream_intron,
                self.downstream_exon,
            ]
            if s
        ]


def extract_se_regions(event: dict, ensembl_db=None) -> Optional[SplicingEventRegions]:
    """
    Extract genomic regions for Skipped Exon (SE) event.

    Region structure:
    [upstream exon][upstream intron][target exon][downstream intron][downstream exon]
    """
    coords = event.get("coordinates", {})
    chromosome = event.get("chr", "")
    strand = event.get("strand", "+")

    upstream_start = coords.get("upstream_start", 0)
    upstream_end = coords.get("upstream_end", 0)
    exon_start = coords.get("start", 0)
    exon_end = coords.get("end", 0)
    downstream_start = coords.get("downstream_start", 0)
    downstream_end = coords.get("downstream_end", 0)

    if not all(
        [
            upstream_start,
            upstream_end,
            exon_start,
            exon_end,
            downstream_start,
            downstream_end,
        ]
    ):
        return None

    regions = SplicingEventRegions(event.get("id", ""), "SE")

    upstream_exon_seq = extract_genomic_region(
        chromosome, upstream_start, upstream_end, strand, ensembl_db=ensembl_db
    )
    upstream_intron_seq = extract_genomic_region(
        chromosome, upstream_end, exon_start, strand, ensembl_db=ensembl_db
    )
    target_exon_seq = extract_genomic_region(
        chromosome, exon_start, exon_end, strand, ensembl_db=ensembl_db
    )
    downstream_intron_seq = extract_genomic_region(
        chromosome, exon_end, downstream_start, strand, ensembl_db=ensembl_db
    )
    downstream_exon_seq = extract_genomic_region(
        chromosome, downstream_start, downstream_end, strand, ensembl_db=ensembl_db
    )

    if upstream_intron_seq and len(upstream_intron_seq) > SPLICE_SITE_3P_MASK:
        upstream_intron_seq = mask_splice_sites(
            upstream_intron_seq, is_upstream=True, strand=strand
        )

    if downstream_intron_seq and len(downstream_intron_seq) > SPLICE_SITE_5P_MASK:
        downstream_intron_seq = mask_splice_sites(
            downstream_intron_seq, is_upstream=False, strand=strand
        )

    regions.upstream_exon = upstream_exon_seq or ""
    regions.upstream_intron = upstream_intron_seq or ""
    regions.target_exon = target_exon_seq or ""
    regions.downstream_intron = downstream_intron_seq or ""
    regions.downstream_exon = downstream_exon_seq or ""

    return regions


def extract_a3ss_regions(
    event: dict, ensembl_db=None
) -> Optional[SplicingEventRegions]:
    """Extract regions for Alternative 3' Splice Site (A3SS) event."""
    coords = event.get("coordinates", {})
    chromosome = event.get("chr", "")
    strand = event.get("strand", "+")

    long_start = coords.get("long_start", 0)
    long_end = coords.get("long_end", 0)
    short_start = coords.get("short_start", 0)
    short_end = coords.get("short_end", 0)
    flanking_start = coords.get("flanking_start", 0)
    flanking_end = coords.get("flanking_end", 0)

    if not all(
        [long_start, long_end, short_start, short_end, flanking_start, flanking_end]
    ):
        return None

    regions = SplicingEventRegions(event.get("id", ""), "A3SS")

    flanking_seq = extract_genomic_region(
        chromosome, flanking_start, flanking_end, strand, ensembl_db=ensembl_db
    )
    long_exon_seq = extract_genomic_region(
        chromosome, long_start, long_end, strand, ensembl_db=ensembl_db
    )
    short_exon_seq = extract_genomic_region(
        chromosome, short_start, short_end, strand, ensembl_db=ensembl_db
    )

    regions.upstream_exon = flanking_seq or ""
    regions.target_exon = long_exon_seq or ""
    regions.downstream_exon = short_exon_seq or ""

    return regions


def extract_a5ss_regions(
    event: dict, ensembl_db=None
) -> Optional[SplicingEventRegions]:
    """Extract regions for Alternative 5' Splice Site (A5SS) event."""
    coords = event.get("coordinates", {})
    chromosome = event.get("chr", "")
    strand = event.get("strand", "+")

    long_start = coords.get("long_start", 0)
    long_end = coords.get("long_end", 0)
    short_start = coords.get("short_start", 0)
    short_end = coords.get("short_end", 0)
    flanking_start = coords.get("flanking_start", 0)
    flanking_end = coords.get("flanking_end", 0)

    if not all(
        [long_start, long_end, short_start, short_end, flanking_start, flanking_end]
    ):
        return None

    regions = SplicingEventRegions(event.get("id", ""), "A5SS")

    long_exon_seq = extract_genomic_region(
        chromosome, long_start, long_end, strand, ensembl_db=ensembl_db
    )
    short_exon_seq = extract_genomic_region(
        chromosome, short_start, short_end, strand, ensembl_db=ensembl_db
    )
    flanking_seq = extract_genomic_region(
        chromosome, flanking_start, flanking_end, strand, ensembl_db=ensembl_db
    )

    regions.upstream_exon = long_exon_seq or ""
    regions.target_exon = short_exon_seq or ""
    regions.downstream_exon = flanking_seq or ""

    return regions


def extract_mxe_regions(event: dict, ensembl_db=None) -> Optional[SplicingEventRegions]:
    """Extract regions for Mutually Exclusive Exons (MXE) event."""
    coords = event.get("coordinates", {})
    chromosome = event.get("chr", "")
    strand = event.get("strand", "+")

    first_start = coords.get("first_start", 0)
    first_end = coords.get("first_end", 0)
    second_start = coords.get("second_start", 0)
    second_end = coords.get("second_end", 0)
    upstream_start = coords.get("upstream_start", 0)
    upstream_end = coords.get("upstream_end", 0)
    downstream_start = coords.get("downstream_start", 0)
    downstream_end = coords.get("downstream_end", 0)

    if not all(
        [
            first_start,
            first_end,
            second_start,
            second_end,
            upstream_start,
            upstream_end,
            downstream_start,
            downstream_end,
        ]
    ):
        return None

    regions = SplicingEventRegions(event.get("id", ""), "MXE")

    upstream_exon_seq = extract_genomic_region(
        chromosome, upstream_start, upstream_end, strand, ensembl_db=ensembl_db
    )
    first_exon_seq = extract_genomic_region(
        chromosome, first_start, first_end, strand, ensembl_db=ensembl_db
    )
    second_exon_seq = extract_genomic_region(
        chromosome, second_start, second_end, strand, ensembl_db=ensembl_db
    )
    downstream_exon_seq = extract_genomic_region(
        chromosome, downstream_start, downstream_end, strand, ensembl_db=ensembl_db
    )

    regions.upstream_exon = upstream_exon_seq or ""
    regions.target_exon = first_exon_seq or ""
    regions.downstream_exon = second_exon_seq or ""
    regions.downstream_intron = downstream_exon_seq or ""

    return regions


def extract_ri_regions(event: dict, ensembl_db=None) -> Optional[SplicingEventRegions]:
    """Extract regions for Retained Intron (RI) event."""
    coords = event.get("coordinates", {})
    chromosome = event.get("chr", "")
    strand = event.get("strand", "+")

    ri_start = coords.get("ri_start", 0)
    ri_end = coords.get("ri_end", 0)
    upstream_start = coords.get("upstream_start", 0)
    upstream_end = coords.get("upstream_end", 0)
    downstream_start = coords.get("downstream_start", 0)
    downstream_end = coords.get("downstream_end", 0)

    if not all(
        [
            ri_start,
            ri_end,
            upstream_start,
            upstream_end,
            downstream_start,
            downstream_end,
        ]
    ):
        return None

    regions = SplicingEventRegions(event.get("id", ""), "RI")

    upstream_exon_seq = extract_genomic_region(
        chromosome, upstream_start, upstream_end, strand, ensembl_db=ensembl_db
    )
    intron_seq = extract_genomic_region(
        chromosome, ri_start, ri_end, strand, ensembl_db=ensembl_db
    )
    downstream_exon_seq = extract_genomic_region(
        chromosome, downstream_start, downstream_end, strand, ensembl_db=ensembl_db
    )

    regions.upstream_exon = upstream_exon_seq or ""
    regions.upstream_intron = intron_seq or ""
    regions.target_exon = downstream_exon_seq or ""

    return regions


def extract_event_regions(
    event: dict, ensembl_db=None
) -> Optional[SplicingEventRegions]:
    """Extract genomic regions for a splicing event based on event type."""
    event_type = event.get("event_type", "SE")

    extractors = {
        "SE": extract_se_regions,
        "A3SS": extract_a3ss_regions,
        "A5SS": extract_a5ss_regions,
        "MXE": extract_mxe_regions,
        "RI": extract_ri_regions,
    }

    extractor = extractors.get(event_type)
    if extractor:
        return extractor(event, ensembl_db)
    return None


def calculate_motif_scores_for_event(
    regions: SplicingEventRegions, rbp_motifs: Dict[str, List[str]] = None
) -> Dict[str, Dict[str, float]]:
    """
    Calculate motif density scores for each region and RBP.

    Returns:
        Dict[RBP_name, Dict[region_name, density_score]]
    """
    if rbp_motifs is None:
        rbp_motifs = RBP_MOTIFS

    scores = {}

    region_map = {
        "upstream_exon": regions.upstream_exon,
        "upstream_intron": regions.upstream_intron,
        "target_exon": regions.target_exon,
        "downstream_intron": regions.downstream_intron,
        "downstream_exon": regions.downstream_exon,
    }

    for rbp_name, motifs in rbp_motifs.items():
        rbp_scores = {}
        for region_name, sequence in region_map.items():
            if sequence:
                densities = [sliding_window_density(sequence, m) for m in motifs]
                max_density = max(max(d) if d else 0 for d in densities)
                rbp_scores[region_name] = max_density
            else:
                rbp_scores[region_name] = 0.0

        scores[rbp_name] = rbp_scores

    return scores


def wilcoxon_rank_sum_test(
    upregulated_scores: List[float], background_scores: List[float]
) -> Tuple[float, float]:
    """
    Perform Wilcoxon rank sum test (Mann-Whitney U test).

    Returns:
        (p_value_upregulated_vs_background, p_value_downregulated_vs_background)
    """
    if len(upregulated_scores) < 3 or len(background_scores) < 3:
        return 1.0, 1.0

    try:
        stat, p_value = stats.mannwhitneyu(
            upregulated_scores, background_scores, alternative="greater"
        )
        return p_value, 1.0
    except Exception:
        return 1.0, 1.0


def run_rmaps_style_analysis(
    upregulated_events: List[dict],
    downregulated_events: List[dict],
    background_events: List[dict],
    sample_size: int = 500,
    flank_size: int = SEQUENCE_FLANK_SIZE,
) -> dict:
    """
    Run RBP motif enrichment analysis following rMAPS method.

    Args:
        upregulated_events: List of upregulated splicing events
        downregulated_events: List of downregulated splicing events
        background_events: List of non-regulated background events
        sample_size: Maximum number of events to analyze per group
        flank_size: Size of flanking intronic regions

    Returns:
        Dictionary with RNA map data and enriched RBPs
    """
    if not upregulated_events and not downregulated_events:
        return {"error": "No events provided", "motifs": [], "rna_map": {}}

    ensembl_db = get_ensembl_db(109, "human")
    if ensembl_db is None:
        return {"error": "Ensembl database not available", "motifs": [], "rna_map": {}}

    up_events = upregulated_events[: min(sample_size, len(upregulated_events))]
    down_events = downregulated_events[: min(sample_size, len(downregulated_events))]
    bg_events = background_events[: min(sample_size * 2, len(background_events))]

    up_regions = []
    down_regions = []
    bg_regions = []

    for event in up_events:
        regions = extract_event_regions(event, ensembl_db)
        if regions and regions.all_sequences():
            up_regions.append(regions)

    for event in down_events:
        regions = extract_event_regions(event, ensembl_db)
        if regions and regions.all_sequences():
            down_regions.append(regions)

    for event in bg_events:
        regions = extract_event_regions(event, ensembl_db)
        if regions and regions.all_sequences():
            bg_regions.append(regions)

    if len(up_regions) < 3 and len(down_regions) < 3:
        return {
            "error": f"Not enough events with extractable sequences (up: {len(up_regions)}, down: {len(down_regions)})",
            "motifs": [],
            "rna_map": {},
        }

    rbp_motifs = RBP_MOTIFS
    region_names = [
        "upstream_exon",
        "upstream_intron",
        "target_exon",
        "downstream_intron",
        "downstream_exon",
    ]

    rna_map = {
        rbp: {r: {"up": [], "down": [], "bg": []} for r in region_names}
        for rbp in rbp_motifs
    }

    for regions in up_regions:
        scores = calculate_motif_scores_for_event(regions, rbp_motifs)
        for rbp_name, rbp_scores in scores.items():
            for region_name in region_names:
                if region_name in rbp_scores:
                    rna_map[rbp_name][region_name]["up"].append(rbp_scores[region_name])

    for regions in down_regions:
        scores = calculate_motif_scores_for_event(regions, rbp_motifs)
        for rbp_name, rbp_scores in scores.items():
            for region_name in region_names:
                if region_name in rbp_scores:
                    rna_map[rbp_name][region_name]["down"].append(
                        rbp_scores[region_name]
                    )

    for regions in bg_regions:
        scores = calculate_motif_scores_for_event(regions, rbp_motifs)
        for rbp_name, rbp_scores in scores.items():
            for region_name in region_names:
                if region_name in rbp_scores:
                    rna_map[rbp_name][region_name]["bg"].append(rbp_scores[region_name])

    enriched_motifs = []

    for rbp_name in rbp_motifs:
        for region_name in region_names:
            up_scores = rna_map[rbp_name][region_name]["up"]
            down_scores = rna_map[rbp_name][region_name]["down"]
            bg_scores = rna_map[rbp_name][region_name]["bg"]

            if len(up_scores) >= 3 and len(bg_scores) >= 3:
                p_val_up, _ = wilcoxon_rank_sum_test(up_scores, bg_scores)

                if p_val_up < 0.05:
                    avg_up = np.mean(up_scores) if up_scores else 0
                    avg_bg = np.mean(bg_scores) if bg_scores else 0
                    fold_change = avg_up / avg_bg if avg_bg > 0 else 0

                    enriched_motifs.append(
                        {
                            "rbp": rbp_name,
                            "region": region_name,
                            "direction": "upregulated",
                            "avg_density_up": round(avg_up, 2),
                            "avg_density_bg": round(avg_bg, 2),
                            "fold_change": round(fold_change, 2),
                            "pvalue": round(p_val_up, 6),
                            "n_up": len(up_scores),
                            "n_bg": len(bg_scores),
                        }
                    )

            if len(down_scores) >= 3 and len(bg_scores) >= 3:
                _, p_val_down = wilcoxon_rank_sum_test(down_scores, bg_scores)

                if p_val_down < 0.05:
                    avg_down = np.mean(down_scores) if down_scores else 0
                    avg_bg = np.mean(bg_scores) if bg_scores else 0
                    fold_change = avg_down / avg_bg if avg_bg > 0 else 0

                    enriched_motifs.append(
                        {
                            "rbp": rbp_name,
                            "region": region_name,
                            "direction": "downregulated",
                            "avg_density_down": round(avg_down, 2),
                            "avg_density_bg": round(avg_bg, 2),
                            "fold_change": round(fold_change, 2),
                            "pvalue": round(p_val_down, 6),
                            "n_down": len(down_scores),
                            "n_bg": len(bg_scores),
                        }
                    )

    enriched_motifs.sort(key=lambda x: x.get("pvalue", 1.0))

    if enriched_motifs:
        pvalues = [m.get("pvalue", 1.0) for m in enriched_motifs]
        adjusted = benjamini_hochberg(pvalues)
        for i, m in enumerate(enriched_motifs):
            m["adj_pvalue"] = round(adjusted[i], 6)

    return {
        "motifs": enriched_motifs[:50],
        "enriched_count": len(enriched_motifs),
        "rna_map": rna_map,
        "method": "rMAPS-style motif enrichment",
        "diagnostics": {
            "up_events_extracted": len(up_regions),
            "down_events_extracted": len(down_regions),
            "bg_events_extracted": len(bg_regions),
            "total_rbps_tested": len(rbp_motifs),
        },
    }


def benjamini_hochberg(p_values: list, alpha: float = 0.05) -> list:
    """Apply Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []

    indexed = [(i, p) for i, p in enumerate(p_values)]
    indexed.sort(key=lambda x: x[1])

    adjusted = [1.0] * n
    for rank, (orig_idx, p) in enumerate(indexed):
        bh_value = p * n / (rank + 1)
        adjusted[orig_idx] = min(1.0, bh_value)

    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    return adjusted


def run_simple_motif_scan(
    sequences: List[str],
    rbp_motifs: Dict[str, List[str]] = None,
    pvalue_threshold: float = 0.05,
) -> dict:
    """
    Simple motif scanning without event structure.
    Useful when only foreground sequences are available.

    Args:
        sequences: List of sequences to scan
        rbp_motifs: Dict of RBP name -> list of motif patterns
        pvalue_threshold: P-value threshold for significance

    Returns:
        Dictionary with enriched motifs
    """
    if rbp_motifs is None:
        rbp_motifs = RBP_MOTIFS

    if not sequences:
        return {"motifs": [], "enriched_count": 0, "method": "Simple motif scan"}

    motif_counts = defaultdict(lambda: {"foreground": 0, "total": 0})

    for sequence in sequences:
        if not sequence or len(sequence) < 5:
            continue
        sequence = sequence.upper()

        for rbp_name, motifs in rbp_motifs.items():
            has_motif = False
            for motif in motifs:
                motif_upper = motif.upper()
                if motif_upper in sequence:
                    has_motif = True
                    break

            motif_counts[rbp_name]["total"] += 1
            if has_motif:
                motif_counts[rbp_name]["foreground"] += 1

    total_seqs = len(sequences)
    background_freq = 0.05

    enrichment_results = []

    for rbp_name, counts in motif_counts.items():
        if counts["total"] < 3:
            continue

        fg_freq = counts["foreground"] / counts["total"]
        expected = background_freq * counts["total"]

        if fg_freq > background_freq * 2:
            contingency = [
                [counts["foreground"], counts["total"] - counts["foreground"]],
                [
                    int(expected),
                    int(
                        counts["total"] * (1 - background_freq)
                        - (
                            expected
                            if expected <= counts["total"] * (1 - background_freq)
                            else 0
                        )
                    ),
                ],
            ]

            try:
                odds_ratio, p_value = stats.fisher_exact(contingency)
                fold_enrichment = (
                    fg_freq / background_freq if background_freq > 0 else 0
                )

                if p_value < pvalue_threshold:
                    enrichment_results.append(
                        {
                            "rbp": rbp_name,
                            "count": counts["foreground"],
                            "total": counts["total"],
                            "percentage": round(fg_freq * 100, 1),
                            "pvalue": round(p_value, 6),
                            "fold_enrichment": round(fold_enrichment, 2),
                        }
                    )
            except Exception:
                pass

    enrichment_results.sort(key=lambda x: x.get("pvalue", 1.0))

    if enrichment_results:
        pvalues = [r.get("pvalue", 1.0) for r in enrichment_results]
        adjusted = benjamini_hochberg(pvalues)
        for i, r in enumerate(enrichment_results):
            r["adj_pvalue"] = round(adjusted[i], 6)

    return {
        "motifs": enrichment_results,
        "enriched_count": len(enrichment_results),
        "method": "Simple motif scan with Fisher's exact test",
    }
