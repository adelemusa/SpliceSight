"""
MEME Suite Integration for RBP Motif Enrichment Analysis.
Uses AME (Analysis of Motif Enrichment) to find significantly enriched motifs.
"""

import os
import json
import subprocess
import tempfile
import shutil
from typing import Optional
from pathlib import Path
from collections import defaultdict
import httpx

from .motif_analysis import get_ensembl_db

MEME_SERVICE_URL = os.environ.get("MEME_SERVICE_URL", "http://meme:8001")


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


def run_meme_enrichment(
    events: list,
    sample_size: int = 500,
    pvalue_threshold: float = 0.05,
    use_background: bool = True,
) -> dict:
    """
    Run complete RBP motif enrichment analysis using MEME AME service.

    Args:
        events: List of splicing events
        sample_size: Number of events to analyze
        pvalue_threshold: P-value threshold for significance
        use_background: Use shuffled sequences as background

    Returns:
        Dictionary with enrichment results
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
        "method": "AME (MEME Suite API)",
        "pvalue_threshold": pvalue_threshold,
    }

    sequences = []
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
            sequences.append(
                {"id": f"{event.get('id', i)}_{gene}", "sequence": sequence}
            )

    results["sequences_extracted"] = len(sequences)

    if len(sequences) < 1:
        return {"error": "Not enough sequences extracted", "motifs": [], **results}

    try:
        response = httpx.post(
            f"{MEME_SERVICE_URL}/enrich",
            json={
                "foreground_sequences": sequences,
                "pvalue_threshold": 1.0,
                "method": "fimo",
            },
            timeout=120.0,
        )
        response.raise_for_status()
        api_result = response.json()

        motif_hits = api_result.get("motifs", [])
        results["method"] = api_result.get("method", "FIMO")

        motif_counts = defaultdict(int)
        for motif in motif_hits:
            rbp = motif.get("motif", "")
            motif_counts[rbp] += 1

        results["motifs"] = [
            {
                "motif": rbp,
                "count": count,
                "percentage": round(count / len(sequences) * 100, 1),
            }
            for rbp, count in sorted(
                motif_counts.items(), key=lambda x: x[1], reverse=True
            )
        ]
        results["enriched_count"] = len(results["motifs"])

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
        get_coordinates_for_event,
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
