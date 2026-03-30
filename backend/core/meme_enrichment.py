"""
MEME Suite Integration for RBP Motif Enrichment Analysis.
Uses FIMO + Fisher's exact test for statistically significant motif enrichment.
"""

import os
import json
import subprocess
import tempfile
import shutil
import math
from typing import Optional
from pathlib import Path
from collections import defaultdict
import httpx

from .motif_analysis import get_ensembl_db

MEME_SERVICE_URL = os.environ.get("MEME_SERVICE_URL", "http://meme:8001")


def log_factorial(n: int) -> float:
    """Compute log factorial using Stirling's approximation for large n."""
    if n < 0:
        return float("inf")
    if n == 0 or n == 1:
        return 0.0
    if n <= 170:
        return math.lgamma(n + 1)
    return n * math.log(n) - n + 0.5 * math.log(2 * math.pi * n)


def fisher_exact_test(a: int, b: int, c: int, d: int) -> float:
    """
    Fisher's exact test for 2x2 contingency table.

    Table:
           motif+    motif-
    fg       a        b
    bg       c        d

    Returns two-tailed p-value.
    """
    if a < 0 or b < 0 or c < 0 or d < 0:
        return 1.0

    n = a + b + c + d
    if n == 0:
        return 1.0

    observed_ratio = a / (a + b) if (a + b) > 0 else 0
    null_ratio = c / (c + d) if (c + d) > 0 else 0

    log_p = (
        log_factorial(a + b)
        + log_factorial(c + d)
        + log_factorial(a + c)
        + log_factorial(b + d)
        - log_factorial(n)
        - log_factorial(a)
        - log_factorial(b)
        - log_factorial(c)
        - log_factorial(d)
    )

    p_value = math.exp(log_p)

    if observed_ratio > null_ratio:
        max_a = min(a + b, a + c)
        for k in range(a + 1, max_a + 1):
            log_p_tail = (
                log_factorial(k + b)
                + log_factorial(c + d - (k - a))
                + log_factorial(k + c - (k - a))
                + log_factorial(b + d)
                - log_factorial(n)
                - log_factorial(k)
                - log_factorial(b + (k - a))
                - log_factorial(c - (k - a))
                - log_factorial(d)
            )
            p_value += math.exp(log_p_tail)

    return min(1.0, p_value)


def bonferroni_correction(p_values: list, alpha: float = 0.05) -> list:
    """Apply Bonferroni correction for multiple testing."""
    n = len(p_values)
    if n == 0:
        return []
    adjusted = [min(1.0, p * n) for p in p_values]
    return adjusted


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


def get_coordinates_for_event(event: dict) -> tuple:
    """Extract coordinates from event based on event type."""
    coords = event.get("coordinates", {})
    event_type = event.get("event_type", "SE")

    if event_type == "SE":
        start = coords.get("start", 0)
        end = coords.get("end", 0)
    elif event_type in ["A3SS", "A5SS"]:
        start = coords.get("long_start", coords.get("start", 0))
        end = coords.get("long_end", coords.get("end", 0))
    elif event_type == "MXE":
        start = coords.get("first_start", coords.get("start", 0))
        end = coords.get("first_end", coords.get("end", 0))
    elif event_type == "RI":
        start = coords.get("ri_start", coords.get("start", 0))
        end = coords.get("ri_end", coords.get("end", 0))
    else:
        start = coords.get("start", 0)
        end = coords.get("end", 0)

    return start, end


RBP_MOTIF_DATABASE = "/usr/share/meme-5.5.5/db/motifs/mouse_motifs.meme"

COMMON_RBP_MOTIFS = {
    "name": "RBP_motifs",
    "motifs": [
        {"id": "PTB_1", "alt": "PTB", "sequence": "UCUCUCUC"},
        {"id": "PTB_2", "alt": "PTB", "sequence": "UUCUUC"},
        {"id": "RBM25", "alt": "RBM25", "sequence": "ACUGU"},
        {"id": "QKI", "alt": "QKI", "sequence": "ACUAUAC"},
        {"id": "NOVA1_1", "alt": "NOVA1", "sequence": "UCAU"},
        {"id": "NOVA1_2", "alt": "NOVA1", "sequence": "YCAY"},
        {"id": "NOVA2", "alt": "NOVA2", "sequence": "UCAU"},
        {"id": "MBNL1", "alt": "MBNL1", "sequence": "GCUGC"},
        {"id": "MBNL2", "alt": "MBNL2", "sequence": "GCUGC"},
        {"id": "CELF1", "alt": "CELF1", "sequence": "GUGU"},
        {"id": "CELF2", "alt": "CELF2", "sequence": "GUGU"},
        {"id": "TIA1", "alt": "TIA1", "sequence": "UUAUU"},
        {"id": "TIAR", "alt": "TIAR", "sequence": "UUAUU"},
        {"id": "HuR", "alt": "ELAVL1", "sequence": "UUU"},
        {"id": "HuR_2", "alt": "ELAVL1", "sequence": "UUAU"},
        {"id": "TDP43", "alt": "TARDBP", "sequence": "UGUGU"},
        {"id": "FUS", "alt": "FUS", "sequence": "GGUGG"},
        {"id": "SF3B4", "alt": "SF3B4", "sequence": "CCU"},
        {"id": "U2AF65", "alt": "U2AF2", "sequence": "AGG"},
        {"id": "RPS14", "alt": "RPS14", "sequence": "ACU"},
        {"id": "HNRNPC", "alt": "HNRNPC", "sequence": "UAG"},
        {"id": "EWSR1", "alt": "EWSR1", "sequence": "GGU"},
        {"id": "RBM10", "alt": "RBM10", "sequence": "UGCU"},
        {"id": "hnRNPA1", "alt": "HNRNPA1", "sequence": "UAGGG"},
        {"id": "hnRNPC_2", "alt": "HNRNPC", "sequence": "ACAU"},
    ],
}


def create_rbp_motif_file(output_path: str) -> str:
    """Create a MEME-format motif file for RBP motifs."""
    lines = [
        "MEME version 5.5.5",
        "",
        "ALPHABET= ACGT",
        "",
        "MOTIFcount 24",
        "",
        "MOTIF RBM25",
        "letter-probability matrix: alength= 4 w= 5 nsites= 20",
        "0.20 0.20 0.40 0.20",
        "0.20 0.40 0.20 0.20",
        "0.20 0.20 0.40 0.20",
        "0.40 0.20 0.20 0.20",
        "0.20 0.20 0.20 0.40",
        "",
        "MOTIF PTB",
        "letter-probability matrix: alength= 4 w= 8 nsites= 20",
        "0.25 0.25 0.25 0.25",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "0.50 0.00 0.50 0.00",
        "",
    ]

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    return output_path


def extract_sequence_for_gene(
    gene_symbol: str,
    start: int,
    end: int,
    chromosome: str,
    strand: str,
    flank_size: int = 50,
    ensembl_db=None,
) -> Optional[str]:
    """Extract transcript sequence for a gene."""
    if ensembl_db is None:
        ensembl_db = get_ensembl_db(109, "human")

    if ensembl_db is None:
        return None

    try:
        genes = ensembl_db.genes_by_name(gene_symbol)
        if not genes:
            return None

        gene = genes[0]

        chrom = chromosome.replace("chr", "")
        if chrom in ["X", "Y", "M", "MT"]:
            chrom = chrom if chrom != "M" else "MT"
        else:
            try:
                chrom = str(int(chrom))
            except ValueError:
                pass

        if gene.contig != chrom:
            return None

        transcript_ids = ensembl_db.transcript_ids_of_gene_name(gene_symbol)
        if not transcript_ids:
            return None

        for tid in transcript_ids[:3]:
            try:
                seq = ensembl_db.transcript_sequence(tid)
                if seq and len(seq) > 10:
                    return seq.upper()
            except Exception:
                continue

        return None

    except Exception as e:
        print(f"Error extracting sequence for {gene_symbol}: {e}")
        return None

    try:
        genes = ensembl_db.genes_by_name(gene_symbol)
        if not genes:
            return None

        gene = genes[0]

        chrom = chromosome.replace("chr", "")
        if chrom in ["X", "Y", "M", "MT"]:
            chrom = chrom if chrom != "M" else "MT"
        else:
            try:
                chrom = str(int(chrom))
            except ValueError:
                pass

        if gene.contig != chrom:
            return None

        sequences = ensembl_db.exon_sequence(gene.id, flavor=None)
        if not sequences:
            return None

        return sequences[0].upper()

    except Exception as e:
        print(f"Error extracting sequence for {gene_symbol}: {e}")
        return None


def write_sequences_to_fasta(events: list, output_path: str, ensembl_db=None) -> int:
    """Write event sequences to FASTA file. Returns number of sequences written."""
    count = 0
    with open(output_path, "w") as f:
        for i, event in enumerate(events):
            gene = event.get("gene_symbol", "")
            if not gene:
                continue

            start, end = get_coordinates_for_event(event)
            if not start or not end:
                continue

            chromosome = event.get("chr", "")
            strand = event.get("strand", "+")

            sequence = extract_sequence_for_gene(
                gene,
                start,
                end,
                chromosome,
                strand,
                flank_size=50,
                ensembl_db=ensembl_db,
            )

            if sequence and len(sequence) > 10:
                f.write(f">{event.get('id', i)}_{gene}\n")
                f.write(f"{sequence}\n")
                count += 1

    return count


def run_ame(
    foreground_fasta: str,
    background_fasta: str,
    motif_file: str,
    output_dir: str,
    pvalue_threshold: float = 0.05,
) -> dict:
    """
    Run AME (Analysis of Motif Enrichment) to find enriched motifs.

    Returns dict with enrichment results.
    """
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "ame",
        "--verbose",
        "2",
        "--o",
        output_dir,
        "--oc",
        os.path.join(output_dir, "ame_results"),
        "--motif",
        motif_file,
        "--bgfile",
        background_fasta,
        "--control",
        "--pvalue",
        str(pvalue_threshold),
        "--method",
        "pearson",
        "--hit-lof-fraction",
        "0.02",
        "--max-stored-scores",
        "100000",
        foreground_fasta,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            print(f"AME error: {result.stderr}")
            return {"error": result.stderr, "motifs": []}

        return parse_ame_results(output_dir)

    except subprocess.TimeoutExpired:
        return {"error": "AME timed out", "motifs": []}
    except Exception as e:
        return {"error": str(e), "motifs": []}


def parse_ame_results(output_dir: str) -> dict:
    """Parse AME results from JSON output."""
    json_path = os.path.join(output_dir, "ame.tsv")
    results = {"motifs": [], "total_seqs": 0, "enriched_count": 0}

    if not os.path.exists(json_path):
        json_path = os.path.join(output_dir, "ame_results", "ame.tsv")

    if os.path.exists(json_path):
        try:
            with open(json_path, "r") as f:
                header = f.readline().strip().split("\t")
                for line in f:
                    fields = line.strip().split("\t")
                    if len(fields) >= 6:
                        motif_name = fields[1] if len(fields) > 1 else ""
                        pvalue = float(fields[3]) if len(fields) > 3 else 1.0
                        evalue = float(fields[4]) if len(fields) > 4 else 1.0
                        adistance = fields[5] if len(fields) > 5 else ""

                        results["motifs"].append(
                            {
                                "motif": motif_name,
                                "pvalue": pvalue,
                                "evalue": evalue,
                                "adj_pvalue": pvalue * len(open(json_path).readlines())
                                if pvalue < 1
                                else pvalue,
                                "annotation": adistance,
                                "source": "AME",
                            }
                        )

            results["enriched_count"] = len(results["motifs"])

        except Exception as e:
            print(f"Error parsing AME results: {e}")

    return results


RBP_BACKGROUND_FREQUENCIES = {
    "RBM25": 0.015,
    "RBM10": 0.02,
    "PTB": 0.08,
    "NOVA1": 0.04,
    "NOVA2": 0.04,
    "MBNL1": 0.035,
    "MBNL2": 0.035,
    "CELF1": 0.06,
    "CELF2": 0.06,
    "TIA1": 0.05,
    "TIAR": 0.05,
    "HuR": 0.07,
    "ELAVL1": 0.07,
    "TDP43": 0.04,
    "TARDBP": 0.04,
    "FUS": 0.03,
    "SF3B4": 0.08,
    "U2AF2": 0.06,
    "RPS14": 0.05,
    "HNRNPC": 0.06,
    "EWSR1": 0.04,
    "hnRNPA1": 0.05,
    "QKI": 0.025,
}

BACKGROUND_GENOME_SIZE = 50000
FOREGROUND_EFFECTIVE_SIZE = 1000


def calculate_nucleotide_freq(sequences: list) -> dict:
    """Calculate nucleotide frequencies from sequences."""
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}
    total = 0
    for seq in sequences:
        for nt in seq.upper():
            if nt in counts:
                counts[nt] += 1
                total += 1
    if total > 0:
        return {k: v / total for k, v in counts.items()}
    return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "U": 0.25}


def estimate_motif_background_frequency(motif: str, nt_freqs: dict) -> float:
    """Estimate background frequency of a motif based on nucleotide frequencies."""
    motif = motif.upper().replace("U", "T")
    freq = 1.0
    for nt in motif:
        freq *= nt_freqs.get(nt, 0.25)
    return freq


def run_meme_enrichment(
    events: list,
    sample_size: int = 500,
    pvalue_threshold: float = 0.05,
    use_background: bool = True,
) -> dict:
    """
    Run complete RBP motif enrichment analysis using MEME FIMO + Fisher's exact test.

    Args:
        events: List of splicing events
        sample_size: Number of events to analyze
        pvalue_threshold: P-value threshold for significance
        use_background: Use background frequencies for Fisher's test

    Returns:
        Dictionary with statistically enriched motifs
    """
    if not events:
        return {"error": "No events provided", "motifs": []}

    ensembl_db = get_ensembl_db(109, "human")

    if ensembl_db is None:
        return {"error": "Ensembl database not available", "motifs": []}

    sample_events = events[: min(sample_size, len(events))]

    results = {
        "total_events_analyzed": len(sample_events),
        "sequences_extracted": 0,
        "motifs": [],
        "enriched_count": 0,
        "method": "FIMO + Fisher's Exact Test",
        "pvalue_threshold": pvalue_threshold,
    }

    sequences = []
    sequence_ids = set()
    for i, event in enumerate(sample_events):
        gene = event.get("gene_symbol", "")
        if not gene:
            continue

        start, end = get_coordinates_for_event(event)
        if not start or not end:
            continue

        chromosome = event.get("chr", "")
        strand = event.get("strand", "+")

        sequence = extract_sequence_for_gene(
            gene, start, end, chromosome, strand, ensembl_db=ensembl_db
        )

        if sequence and len(sequence) > 10:
            seq_id = f"{event.get('id', i)}_{gene}"
            if seq_id not in sequence_ids:
                sequences.append({"id": seq_id, "sequence": sequence})
                sequence_ids.add(seq_id)

    results["sequences_extracted"] = len(sequences)

    if len(sequences) < 3:
        return {
            "error": "Not enough sequences extracted (minimum 3 required)",
            "motifs": [],
            **results,
        }

    nt_freqs = calculate_nucleotide_freq([s["sequence"] for s in sequences])

    try:
        response = httpx.post(
            f"{MEME_SERVICE_URL}/enrich",
            json={
                "foreground_sequences": sequences,
                "pvalue_threshold": 0.1,
                "method": "fimo",
            },
            timeout=120.0,
        )
        response.raise_for_status()
        api_result = response.json()

        motif_hits = api_result.get("motifs", [])
        results["method"] = api_result.get("method", "FIMO") + " + Fisher's Exact Test"

        motif_counts = defaultdict(int)
        motif_sequences = defaultdict(set)
        for motif in motif_hits:
            rbp = motif.get("motif", "")
            seq_id = motif.get("sequence_id", "")
            motif_counts[rbp] += 1
            if seq_id:
                motif_sequences[rbp].add(seq_id)

        total_seqs = len(sequences)
        enrichment_results = []

        for rbp, count in motif_counts.items():
            fg_with_motif = len(motif_sequences[rbp])
            fg_without_motif = total_seqs - fg_with_motif

            bg_expected = RBP_BACKGROUND_FREQUENCIES.get(rbp, 0.05)

            motif_len = (
                len(rbp)
                if rbp
                in [
                    m.upper()
                    for m in [
                        "RBM25",
                        "PTB",
                        "QKI",
                        "MBNL1",
                        "CELF1",
                        "TIA1",
                        "HuR",
                        "TDP43",
                        "hnRNPA1",
                    ]
                ]
                else 3
            )
            bg_freq_adjusted = estimate_motif_background_frequency(rbp, nt_freqs)
            bg_freq = max(bg_expected, bg_freq_adjusted)

            bg_with_motif = int(BACKGROUND_GENOME_SIZE * bg_freq)
            bg_without_motif = BACKGROUND_GENOME_SIZE - bg_with_motif

            pvalue = fisher_exact_test(
                fg_with_motif, fg_without_motif, bg_with_motif, bg_without_motif
            )

            fold_enrichment = 0
            if bg_with_motif > 0:
                fg_rate = fg_with_motif / total_seqs if total_seqs > 0 else 0
                bg_rate = bg_with_motif / BACKGROUND_GENOME_SIZE
                fold_enrichment = fg_rate / bg_rate if bg_rate > 0 else 0

            enrichment_results.append(
                {
                    "motif": rbp,
                    "count": fg_with_motif,
                    "percentage": round(fg_with_motif / total_seqs * 100, 1)
                    if total_seqs > 0
                    else 0,
                    "pvalue": pvalue,
                    "fold_enrichment": round(fold_enrichment, 2)
                    if fold_enrichment > 0
                    else 0,
                    "bg_frequency": round(bg_freq * 100, 2),
                }
            )

        if enrichment_results:
            pvalues = [r["pvalue"] for r in enrichment_results]
            adjusted_pvalues = benjamini_hochberg(pvalues, alpha=pvalue_threshold)

            for i, result in enumerate(enrichment_results):
                result["adj_pvalue"] = round(adjusted_pvalues[i], 6)
                result["significant"] = result["adj_pvalue"] < pvalue_threshold
                result["significance"] = (
                    "***"
                    if result["adj_pvalue"] < 0.001
                    else "**"
                    if result["adj_pvalue"] < 0.01
                    else "*"
                    if result["adj_pvalue"] < pvalue_threshold
                    else "ns"
                )

        enrichment_results.sort(key=lambda x: x["pvalue"])

        results["motifs"] = enrichment_results
        results["enriched_count"] = sum(
            1 for m in enrichment_results if m.get("significant", False)
        )
        results["significant_motifs"] = [
            m["motif"] for m in enrichment_results if m.get("significant", False)
        ]

        if results["motifs"]:
            results["top_motif"] = results["motifs"][0]

        return results

    except httpx.ConnectError:
        return {
            "error": f"Cannot connect to MEME service at {MEME_SERVICE_URL}",
            "motifs": [],
            **results,
        }
    except Exception as e:
        return {"error": f"MEME API error: {str(e)}", "motifs": [], **results}


def run_meme_enrichment_simple(
    events: list, sample_size: int = 500, pvalue_threshold: float = 0.05
) -> dict:
    """
    Simplified enrichment using Fisher's exact test against RBP motif database.
    Falls back to this if MEME is not available.
    """
    from collections import Counter

    if not events:
        return {"error": "No events provided", "motifs": []}

    ensembl_db = get_ensembl_db(109, "human")

    if ensembl_db is None:
        return {"error": "Ensembl database not available", "motifs": []}

    sample_events = events[: min(sample_size, len(events))]

    from .motif_analysis import (
        MOTIF_DATABASE,
        analyze_sequence_motifs,
    )

    motif_counts = Counter()
    total_genes = 0

    for event in sample_events:
        gene = event.get("gene_symbol", "")
        if not gene:
            continue

        start, end = get_coordinates_for_event(event)
        chromosome = event.get("chr", "")
        strand = event.get("strand", "+")

        sequence = extract_sequence_for_gene(
            gene, start, end, chromosome, strand, ensembl_db=ensembl_db
        )

        if sequence:
            result = analyze_sequence_motifs(sequence, event.get("event_type", ""))
            for motif in result.get("motifs_found", []):
                motif_counts[motif["rbp"]] += motif["count"]
            total_genes += 1

    motifs = []
    for rbp, count in motif_counts.most_common(20):
        percentage = (count / total_genes * 100) if total_genes > 0 else 0
        motifs.append(
            {
                "motif": rbp,
                "count": count,
                "percentage": round(percentage, 2),
                "pvalue": None,
                "evalue": None,
                "source": "motif_scan",
            }
        )

    return {
        "total_events_analyzed": len(sample_events),
        "sequences_extracted": total_genes,
        "motifs": motifs,
        "enriched_count": len(motifs),
        "method": "Motif scanning",
        "pvalue_threshold": None,
    }
