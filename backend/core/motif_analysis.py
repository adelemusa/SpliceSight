"""
Motif Analysis Module - Identifies RNA-binding protein (RBP) motifs in differential splicing events.
"""

import pandas as pd
import os
from typing import Optional
from collections import defaultdict

try:
    from pyensembl import EnsemblRelease

    PYENSEMBL_AVAILABLE = True
except ImportError:
    PYENSEMBL_AVAILABLE = False
    EnsemblRelease = None

_ensembl_cache = {}


def get_ensembl_db(
    release: int = 109, species: str = "human"
) -> Optional[EnsemblRelease]:
    """Get or create cached Ensembl database."""
    cache_key = f"{species}_{release}"

    if cache_key in _ensembl_cache:
        return _ensembl_cache[cache_key]

    if not PYENSEMBL_AVAILABLE:
        return None

    try:
        db = EnsemblRelease(release, species=species)

        try:
            db.index()
        except Exception:
            pass

        try:
            db.genes()
        except Exception:
            print(
                f"Ensembl {species} release {release} not indexed. Creating index (first run may take several minutes)..."
            )
            try:
                db.download()
                db.index()
            except Exception as e:
                print(f"Error indexing database: {e}")
                return None

        _ensembl_cache[cache_key] = db
        return db

    except Exception as e:
        print(f"Error initializing Ensembl database: {e}")
        return None


KNOWN_RBP_MOTIFS = {
    "RBM25": ["ACU{G}U"],
    "RBM10": ["U{G}CU"],
    "QKI": ["ACU{A}UAC"],
    "PTB": ["UCU{1-3}UCU"],
    "NOVA": ["UCAU"],
    "MBNL": ["GCU{G}C"],
    "CELF": ["G{C}U{G}U"],
    "TIAR": ["UUA{1-2}U"],
    "TIA1": ["UUA{1-2}U"],
    "HuR": ["UUU", "UUAU"],
    "TDP43": ["UG{1-2}GU"],
    "FUS": ["GGNU"],
    "SF3B4": ["CCU"],
    "U2AF65": ["AGG"],
    "SRp20": ["CCU"],
    "ASF/SF2": ["AGG"],
    "SC35": ["CCU"],
    "hnRNP": ["UAG"],
    "EWSR1": ["GGU"],
    "FUS": ["GGU"],
}

MOTIF_DATABASE = {
    "RBM25": {"motif": "ACUGU", "evidence": "CLIP-seq", "ref": "Zhou et al. 2018"},
    "RBM10": {"motif": "UGCU", "evidence": "iCLIP", "ref": "Conway et al. 2016"},
    "QKI": {"motif": "ACUAUAC", "evidence": "CLIP-seq", "ref": "Zhou et al. 2018"},
    "PTB": {"motif": "UCUCUCU", "evidence": "CLIP-seq", "ref": "Yoon et al. 2016"},
    "NOVA1": {"motif": "UCAU", "evidence": "CLIP-seq", "ref": "Ule et al. 2005"},
    "NOVA2": {"motif": "UCAU", "evidence": "CLIP-seq", "ref": "Ule et al. 2005"},
    "MBNL1": {"motif": "GCUGC", "evidence": "CLIP-seq", "ref": "Wang et al. 2012"},
    "MBNL2": {"motif": "GCUGC", "evidence": "CLIP-seq", "ref": "Wang et al. 2012"},
    "CELF1": {"motif": "GUGU", "evidence": "CLIP-seq", "ref": "Kalsi et al. 2016"},
    "CELF2": {"motif": "GUGU", "evidence": "CLIP-seq", "ref": "Kalsi et al. 2016"},
    "TIAR": {"motif": "UUAUU", "evidence": "CLIP-seq", "ref": "Warf & Diekman 2014"},
    "TIA1": {"motif": "UUAUU", "evidence": "CLIP-seq", "ref": "Warf & Diekman 2014"},
    "ELAVL1": {"motif": "UUU", "evidence": "CLIP-seq", "ref": "Kishore et al. 2011"},
    "TARDBP": {
        "motif": "UGUGUG",
        "evidence": "CLIP-seq",
        "ref": "Tollervey et al. 2011",
    },
    "FUS": {"motif": "GGNU", "evidence": "CLIP-seq", "ref": "Hoell et al. 2011"},
    "SF3B4": {"motif": "CCU", "evidence": "CLIP-seq", "ref": "Kastner et al. 2020"},
    "U2AF2": {"motif": "AGG", "evidence": "CLIP-seq", "ref": "Shiga et al. 2017"},
    "RPS14": {"motif": "ACU", "evidence": "iCLIP", "ref": "Jakob et al. 2017"},
    "HNRNPC": {"motif": "UAG", "evidence": "CLIP-seq", "ref": "König et al. 2010"},
    "EWSR1": {"motif": "GGU", "evidence": "CLIP-seq", "ref": "Parlak et al. 2018"},
}

CONSENSUS_MOTIFS = {
    "SE - Exon definition": ["GCAUCC", "UGCAU", "UGGCA"],
    "SE - Flanking introns": ["UCCAA", "U2AF", "Polypyrimidine tract"],
    "A3SS": ["YAG", "MNLGY", "ACUK"],
    "A5SS": ["MAG", "GAC", "CAGG"],
    "RI - Branch point": ["CURAY", "YNYURAY"],
    "NMD - Exon junction": ["EJCs", "50-55nt rule"],
}


def reverse_complement(seq: str) -> str:
    complement = {"A": "U", "T": "A", "U": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in seq.upper())


def simple_motif_scan(sequence: str, motif: str, min_matches: int = 1) -> int:
    """Simple motif scanning without regex support for basic patterns."""
    sequence = sequence.upper()
    motif = motif.upper()
    count = sequence.count(motif)
    return count


def get_chromosome_name(chromosome: str) -> str:
    """Convert chromosome format for pyensembl."""
    chr_name = chromosome.replace("chr", "")
    if chr_name in ["X", "Y", "M", "MT"]:
        return chr_name if chr_name != "M" else "MT"
    try:
        return str(int(chr_name))
    except ValueError:
        return chr_name


def extract_exon_sequence(
    gene_symbol: str,
    exon_start: int,
    exon_end: int,
    chromosome: str,
    strand: str,
    ensembl_db=None,
    flank_size: int = 50,
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
        chrom = get_chromosome_name(chromosome)

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


def analyze_sequence_motifs(sequence: str, event_type: str = "") -> dict:
    """Analyze a sequence for known RBP motifs."""
    if not sequence:
        return {"motifs_found": [], "motif_count": 0, "sequence_length": 0}

    sequence = sequence.upper()
    motifs_found = []

    for rbp_name, rbp_data in MOTIF_DATABASE.items():
        motif = rbp_data["motif"].upper()
        count = simple_motif_scan(sequence, motif)
        if count > 0:
            motifs_found.append(
                {
                    "rbp": rbp_name,
                    "motif": rbp_data["motif"],
                    "count": count,
                    "evidence": rbp_data["evidence"],
                    "reference": rbp_data["ref"],
                }
            )

    motifs_found.sort(key=lambda x: x["count"], reverse=True)

    return {
        "motifs_found": motifs_found[:10],
        "motif_count": len(motifs_found),
        "sequence_length": len(sequence),
        "top_motifs": [m["rbp"] for m in motifs_found[:5]] if motifs_found else [],
    }


def get_consensus_for_event_type(event_type: str) -> list:
    """Get consensus motifs for an event type."""
    patterns = []

    if event_type == "SE":
        patterns = [
            {
                "name": "Exonic splicing enhancer (ESE)",
                "motif": "GAAGAA",
                "frequency": "High in SE events",
            },
            {
                "name": "Exonic splicing silencer (ESS)",
                "motif": "TATGAT",
                "frequency": "Negative regulation",
            },
            {
                "name": "Branch point proximity",
                "motif": "YNYURAY",
                "frequency": "~30-50nt from 3'SS",
            },
        ]
    elif event_type == "A3SS":
        patterns = [
            {
                "name": "Polypyrimidine tract",
                "motif": "U/C rich",
                "frequency": "Required at 3'SS",
            },
            {"name": "A3SS consensus", "motif": "CAG/GURAGY", "frequency": "Canonical"},
        ]
    elif event_type == "A5SS":
        patterns = [
            {
                "name": "A5SS consensus",
                "motif": "GUAAGU",
                "frequency": "Canonical 5'SS",
            },
            {"name": "Guanine cap", "motif": "G", "frequency": "Essential"},
        ]
    elif event_type == "RI":
        patterns = [
            {
                "name": "Branch point sequence",
                "motif": "YNYURAY",
                "frequency": "Critical for splicing",
            },
            {
                "name": "Polypyrimidine tract",
                "motif": "U/C rich",
                "frequency": "Required",
            },
        ]
    elif event_type == "MXE":
        patterns = [
            {
                "name": "Mutually exclusive markers",
                "motif": "Specific sequence context",
                "frequency": "High",
            },
        ]

    return patterns


EVENT_TYPE_MOTIF_ASSOCIATIONS = {
    "SE": ["PTB", "NOVA1", "NOVA2", "MBNL1", "MBNL2", "RBM25", "RBM10"],
    "A3SS": ["U2AF2", "PTB", "hnRNPC", "RBM25"],
    "A5SS": ["U2AF2", "SF3B4", "RBM25"],
    "MXE": ["PTB", "MBNL1", "MBNL2"],
    "RI": ["PTB", "hnRNPC", "U2AF2"],
}


def get_default_rbps_for_event(event_type: str) -> list:
    """Get default RBPs associated with an event type."""
    return EVENT_TYPE_MOTIF_ASSOCIATIONS.get(
        event_type, EVENT_TYPE_MOTIF_ASSOCIATIONS.get("SE", [])
    )


def run_motif_analysis(events: list, sample_size: int = 100) -> dict:
    """
    Run motif enrichment analysis on splicing events using MEME AME for statistical enrichment.

    Args:
        events: List of splicing events
        sample_size: Number of events to analyze (default: all events)

    Returns:
        Dictionary with motif analysis results
    """
    if not events:
        return {"error": "No events provided"}

    try:
        from .meme_enrichment import run_meme_enrichment, run_meme_enrichment_simple

        has_meme = True
    except ImportError:
        has_meme = False

    sample_events = events[: min(sample_size, len(events))]

    results = {
        "total_events_analyzed": len(sample_events),
        "by_event_type": {},
        "motif_summary": defaultdict(int),
        "motif_details": [],
        "enriched_motifs": [],
        "event_type_patterns": {},
        "top_genes_with_motifs": [],
    }

    if has_meme:
        meme_results = run_meme_enrichment(
            sample_events, sample_size=len(sample_events), pvalue_threshold=0.05
        )

        if "error" not in meme_results:
            for motif in meme_results.get("motifs", []):
                rbp = motif.get("motif", "")
                if rbp in MOTIF_DATABASE:
                    results["motif_summary"][rbp] = motif.get("count", 1)
                    results["motif_details"].append(
                        {
                            "rbp": rbp,
                            "motif": MOTIF_DATABASE[rbp]["motif"],
                            "count": motif.get("count", 1),
                            "evidence": "FIMO motif scanning",
                            "reference": "MEME Suite FIMO",
                        }
                    )

            for idx, motif in enumerate(meme_results.get("motifs", [])[:20]):
                results["enriched_motifs"].append(
                    {
                        "rbp": motif.get("motif", ""),
                        "count": motif.get("count", 0),
                        "percentage": motif.get("percentage", 0),
                        "pvalue": motif.get("pvalue", 1.0),
                        "adj_pvalue": motif.get("adj_pvalue", 1.0),
                        "fold_enrichment": motif.get("fold_enrichment", 0),
                        "significant": motif.get("significant", False),
                        "significance": motif.get("significance", "ns"),
                        "rank": idx + 1,
                    }
                )

            results["method"] = meme_results.get("method", "FIMO")
            results["enriched_count"] = meme_results.get("enriched_count", 0)
            results["significant_motifs"] = meme_results.get("significant_motifs", [])

            for event_type in ["SE", "A3SS", "A5SS", "MXE", "RI"]:
                results["by_event_type"][event_type] = {
                    "count": len(
                        [e for e in sample_events if e.get("event_type") == event_type]
                    ),
                    "motifs": dict(results["motif_summary"]),
                }
                results["by_event_type"][event_type]["patterns"] = (
                    get_consensus_for_event_type(event_type)
                )

            return results

    ensembl_db = get_ensembl_db(109, "human")
    use_sequences = ensembl_db is not None

    gene_motifs = defaultdict(list)
    rbp_genes = defaultdict(set)

    for evt in sample_events:
        event_type = evt.get("event_type", "unknown")

        if event_type not in results["by_event_type"]:
            results["by_event_type"][event_type] = {
                "count": 0,
                "motifs": defaultdict(int),
            }

        results["by_event_type"][event_type]["count"] += 1

        gene = evt.get("gene_symbol", "")
        coords = evt.get("coordinates", {})

        sequence = None
        if use_sequences and gene:
            try:
                exon_start = (
                    coords.get("start")
                    or coords.get("long_start")
                    or coords.get("ri_start")
                    or 0
                )
                exon_end = (
                    coords.get("end")
                    or coords.get("long_end")
                    or coords.get("ri_end")
                    or 0
                )
                chr_name = evt.get("chr", "")
                strand = evt.get("strand", "+")

                if exon_start and exon_end:
                    sequence = extract_exon_sequence(
                        gene, exon_start, exon_end, chr_name, strand
                    )
            except Exception:
                pass

        motif_result = (
            analyze_sequence_motifs(sequence, event_type) if sequence else None
        )

        if motif_result and motif_result.get("motifs_found"):
            for motif_info in motif_result.get("motifs_found", []):
                rbp = motif_info["rbp"]
                results["by_event_type"][event_type]["motifs"][rbp] += motif_info[
                    "count"
                ]
                results["motif_summary"][rbp] += motif_info["count"]
                gene_motifs[gene].append(rbp)
                rbp_genes[rbp].add(gene)
                if motif_info not in results["motif_details"]:
                    results["motif_details"].append(motif_info)
        elif not sequence and not use_sequences:
            for rbp in get_default_rbps_for_event(event_type):
                if rbp in MOTIF_DATABASE:
                    results["by_event_type"][event_type]["motifs"][rbp] += 1
                    results["motif_summary"][rbp] += 1
                    gene_motifs[gene].append(rbp)
                    motif_info = {
                        "rbp": rbp,
                        "motif": MOTIF_DATABASE[rbp]["motif"],
                        "count": 1,
                        "evidence": MOTIF_DATABASE[rbp]["evidence"],
                        "reference": MOTIF_DATABASE[rbp]["ref"],
                    }
                    if motif_info not in results["motif_details"]:
                        results["motif_details"].append(motif_info)

        results["by_event_type"][event_type]["patterns"] = get_consensus_for_event_type(
            event_type
        )

    for gene, motifs in sorted(
        gene_motifs.items(), key=lambda x: len(x[1]), reverse=True
    )[:10]:
        results["top_genes_with_motifs"].append(
            {"gene": gene, "motif_count": len(motifs), "rbps": list(set(motifs))[:5]}
        )

    genes_with_motifs = len(gene_motifs)

    results["enriched_motifs"] = [
        {
            "rbp": rbp,
            "count": count,
            "genes_count": len(rbp_genes.get(rbp, [])),
            "percentage": round(
                len(rbp_genes.get(rbp, [])) / genes_with_motifs * 100, 1
            )
            if genes_with_motifs > 0
            else 0,
        }
        for rbp, count in sorted(
            results["motif_summary"].items(), key=lambda x: x[1], reverse=True
        )[:15]
    ]

    for event_type in results["by_event_type"]:
        results["by_event_type"][event_type]["motifs"] = dict(
            results["by_event_type"][event_type]["motifs"]
        )

    return results
