import pandas as pd
import os
from dataclasses import dataclass, field
from typing import Optional, Union

try:
    from .genome_manager import get_ensembl
except ImportError:
    get_ensembl = None

try:
    from .enrichment import run_enrichment
except ImportError:

    def run_enrichment(gene_symbols, gene_sets=None, pvalue_cutoff=0.1):
        return {"pathways": [], "gene_count": 0}


try:
    from .motif_analysis import run_motif_analysis
except ImportError:

    def run_motif_analysis(events, sample_size=100):
        return {}


EVENT_TYPES = {
    "SE": {
        "pattern": [
            "SE.MATS.JC",
            "SE.MATS.JCEC",
            "fromGTF.SE",
            "fromGTF.novelSpliceSite.SE",
            "fromGTF.novelJunction.SE",
        ],
        "exon_cols": {
            "start": "exonStart_0base",
            "end": "exonEnd",
            "upstream_start": "upstreamES",
            "upstream_end": "upstreamEE",
            "downstream_start": "downstreamES",
            "downstream_end": "downstreamEE",
        },
        "description": "Skipped Exon",
    },
    "A3SS": {
        "pattern": [
            "A3SS.MATS.JC",
            "A3SS.MATS.JCEC",
            "fromGTF.A3SS",
            "fromGTF.novelSpliceSite.A3SS",
            "fromGTF.novelJunction.A3SS",
        ],
        "exon_cols": {
            "long_start": "longExonStart_0base",
            "long_end": "longExonEnd",
            "short_start": "shortES",
            "short_end": "shortEE",
            "flanking_start": "flankingES",
            "flanking_end": "flankingEE",
        },
        "description": "Alternative 3' Splice Site",
    },
    "A5SS": {
        "pattern": [
            "A5SS.MATS.JC",
            "A5SS.MATS.JCEC",
            "fromGTF.A5SS",
            "fromGTF.novelSpliceSite.A5SS",
            "fromGTF.novelJunction.A5SS",
        ],
        "exon_cols": {
            "long_start": "longExonStart_0base",
            "long_end": "longExonEnd",
            "short_start": "shortES",
            "short_end": "shortEE",
            "flanking_start": "flankingES",
            "flanking_end": "flankingEE",
        },
        "description": "Alternative 5' Splice Site",
    },
    "MXE": {
        "pattern": [
            "MXE.MATS.JC",
            "MXE.MATS.JCEC",
            "fromGTF.MXE",
            "fromGTF.novelSpliceSite.MXE",
            "fromGTF.novelJunction.MXE",
        ],
        "exon_cols": {
            "first_start": "1stExonStart_0base",
            "first_end": "1stExonEnd",
            "second_start": "2ndExonStart_0base",
            "second_end": "2ndExonEnd",
            "upstream_start": "upstreamES",
            "upstream_end": "upstreamEE",
            "downstream_start": "downstreamES",
            "downstream_end": "downstreamEE",
        },
        "description": "Mutually Exclusive Exons",
    },
    "RI": {
        "pattern": [
            "RI.MATS.JC",
            "RI.MATS.JCEC",
            "fromGTF.RI",
            "fromGTF.novelSpliceSite.RI",
            "fromGTF.novelJunction.RI",
        ],
        "exon_cols": {
            "ri_start": "riExonStart_0base",
            "ri_end": "riExonEnd",
            "upstream_start": "upstreamES",
            "upstream_end": "upstreamEE",
            "downstream_start": "downstreamES",
            "downstream_end": "downstreamEE",
        },
        "description": "Retained Intron",
    },
}


@dataclass
class RMatsConfig:
    fdr_threshold: float = 0.05
    dpsi_threshold: float = 0.1
    include_jc: bool = True
    include_jcec: bool = True
    species: str = "human"
    ensembl_release: int = 109


@dataclass
class SplicingEvent:
    event_id: str
    event_type: str
    gene_symbol: str
    gene_id: str
    chromosome: str
    strand: str
    fdr: float
    dpsi: float
    p_value: float
    inc_level_1: Optional[str] = None
    inc_level_2: Optional[str] = None
    coordinates: dict = field(default_factory=dict)
    proteome_impact: Optional[str] = None
    nmd_candidate: bool = False
    coordinates_1: Optional[str] = None
    coordinates_2: Optional[str] = None
    sample_info: dict = field(default_factory=dict)


def detect_event_type(filename: str) -> Optional[str]:
    for event_type, config in EVENT_TYPES.items():
        for pattern in config["pattern"]:
            if pattern in filename:
                return event_type
    return None


def get_file_type(filename: str) -> str:
    if "JCEC" in filename:
        return "JCEC"
    return "JC"


def parse_inc_level(value: str) -> list:
    if not value or value == "NA":
        return []
    try:
        return [float(x) if x != "NA" else None for x in value.split(",")]
    except (ValueError, AttributeError):
        return []


def check_proteome_impact(
    gene_symbol: str, exon_length: int, is_terminal: bool = False
) -> tuple:
    if exon_length % 3 != 0:
        impact = f"Frameshift ({exon_length}bp)"
        is_nmd = not is_terminal
    else:
        impact = f"In-frame indel ({exon_length}bp)"
        is_nmd = False
    return impact, is_nmd


def extract_coordinates(row: dict, event_type: str, cols_config: dict) -> dict:
    coords = {}
    for key, col in cols_config.items():
        if col in row:
            val = row[col]
            try:
                coords[key] = int(float(val)) if val is not None else 0
            except (ValueError, TypeError):
                pass
    return coords


def parse_single_file(file_path: str, config: RMatsConfig, ensembl_db=None) -> tuple:
    events = []
    significant_genes = set()

    try:
        df = pd.read_csv(file_path, sep="\t")

        required_cols = ["FDR", "IncLevelDifference"]
        if not all(col in df.columns for col in required_cols):
            return events, significant_genes

        df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
        df["IncLevelDifference"] = pd.to_numeric(
            df["IncLevelDifference"], errors="coerce"
        )

        filtered_df = df[
            (df["FDR"] < config.fdr_threshold)
            & (df["IncLevelDifference"].abs() > config.dpsi_threshold)
        ]

        event_type = detect_event_type(os.path.basename(file_path))
        file_type = get_file_type(file_path)

        if event_type is None:
            return events, significant_genes

        cols_config = EVENT_TYPES[event_type]["exon_cols"]

        for idx in filtered_df.index:
            row_dict = filtered_df.loc[idx].to_dict()
            gene_symbol = str(row_dict.get("geneSymbol", "Unknown") or "Unknown").strip(
                '"'
            )
            gene_id = str(row_dict.get("GeneID", "") or "").strip('"')
            chromosome = str(row_dict.get("chr", "") or "").replace("chr", "")
            strand = str(row_dict.get("strand", "+") or "+")

            significant_genes.add(gene_symbol)

            coords = extract_coordinates(row_dict, event_type, cols_config)

            exon_length = 0
            if event_type == "SE" and "start" in coords and "end" in coords:
                exon_length = coords["end"] - coords["start"]
            elif (
                event_type in ("A3SS", "A5SS")
                and "long_start" in coords
                and "long_end" in coords
            ):
                exon_length = coords["long_end"] - coords["long_start"]
            elif (
                event_type == "MXE"
                and "first_start" in coords
                and "first_end" in coords
            ):
                exon_length = coords["first_end"] - coords["first_start"]
            elif event_type == "RI" and "ri_start" in coords and "ri_end" in coords:
                exon_length = coords["ri_end"] - coords["ri_start"]

            impact, is_nmd = check_proteome_impact(gene_symbol, exon_length)

            sample_info = {}
            for col in df.columns:
                if col.startswith(("IJC_", "SJC_", "IncLevel")):
                    sample_info[col] = str(row_dict.get(col) or "")

            fdr_val = row_dict.get("FDR")
            dpsi_val = row_dict.get("IncLevelDifference")
            pval_val = row_dict.get("PValue", 1.0)
            inc1_val = row_dict.get("IncLevel1", "")
            inc2_val = row_dict.get("IncLevel2", "")

            event = SplicingEvent(
                event_id=str(row_dict.get("ID") or ""),
                event_type=event_type,
                gene_symbol=gene_symbol,
                gene_id=gene_id,
                chromosome=chromosome,
                strand=strand,
                fdr=float(fdr_val) if fdr_val is not None else 1.0,
                dpsi=float(dpsi_val) if dpsi_val is not None else 0.0,
                p_value=float(pval_val) if pval_val is not None else 1.0,
                inc_level_1=str(inc1_val) if inc1_val else "",
                inc_level_2=str(inc2_val) if inc2_val else "",
                coordinates=coords,
                proteome_impact=impact,
                nmd_candidate=is_nmd,
                sample_info=sample_info,
            )
            events.append(event)

    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")

    return events, significant_genes


def parse_rmats_files(
    input_path: str,
    fdr_threshold: float = 0.05,
    dpsi_threshold: float = 0.1,
    include_jc: bool = True,
    include_jcec: bool = True,
    config: Union[RMatsConfig, None] = None,
) -> dict:
    if config is None:
        config = RMatsConfig(
            fdr_threshold=fdr_threshold,
            dpsi_threshold=dpsi_threshold,
            include_jc=include_jc,
            include_jcec=include_jcec,
        )

    events = []
    significant_genes = set()

    if os.path.isfile(input_path):
        file_events, file_genes = parse_single_file(input_path, config)
        events.extend(file_events)
        significant_genes.update(file_genes)

    elif os.path.isdir(input_path):
        processed_types = set()

        file_names = sorted(
            os.listdir(input_path), key=lambda x: (0 if ".MATS." in x else 1, x)
        )

        for file_name in file_names:
            file_path = os.path.join(input_path, file_name)

            if not os.path.isfile(file_path):
                continue

            event_type = detect_event_type(file_name)
            file_type = get_file_type(file_name)

            if event_type is None:
                continue

            if file_type == "JCEC" and not config.include_jcec:
                continue
            if file_type == "JC" and not config.include_jc:
                continue

            type_key = (event_type, file_type)
            if type_key in processed_types:
                continue
            processed_types.add(type_key)

            file_events, file_genes = parse_single_file(file_path, config)
            events.extend(file_events)
            significant_genes.update(file_genes)

    enrichment_results = run_enrichment(list(significant_genes), pvalue_cutoff=0.1)

    events_data = [
        {
            "id": e.event_id,
            "event_type": e.event_type,
            "gene_symbol": e.gene_symbol,
            "gene_id": e.gene_id,
            "chr": e.chromosome,
            "strand": e.strand,
            "fdr": e.fdr,
            "dpsi": e.dpsi,
            "p_value": e.p_value,
            "coordinates": e.coordinates,
            "proteome_impact": e.proteome_impact,
            "nmd_candidate": e.nmd_candidate,
            "sample_info": e.sample_info,
        }
        for e in events
    ]

    motif_results = run_motif_analysis(events_data, sample_size=len(events_data))

    return {
        "events": events_data,
        "enrichment": enrichment_results,
        "motif_analysis": motif_results,
        "summary": {
            "total_events": len(events),
            "significant_genes": len(significant_genes),
            "by_event_type": _count_by_type(events),
        },
    }


def _count_by_type(events: list) -> dict:
    counts = {}
    for event in events:
        etype = event.event_type
        counts[etype] = counts.get(etype, 0) + 1
    return counts


def get_available_event_types() -> list:
    return [
        {
            "type": etype,
            "description": config["description"],
            "patterns": config["pattern"],
        }
        for etype, config in EVENT_TYPES.items()
    ]
