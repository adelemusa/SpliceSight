"""
MEME Service - REST API for MEME Suite tools.
Provides AME (Analysis of Motif Enrichment) for RBP motif analysis.
"""

import os
import subprocess
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import uvicorn

app = FastAPI(title="MEME Service")


class MotifEnrichmentRequest(BaseModel):
    foreground_sequences: List[dict]  # [{"id": str, "sequence": str}]
    background_sequences: Optional[List[dict]] = None
    pvalue_threshold: float = 0.05
    method: str = "ame"  # ame or fimo


class MotifEnrichmentResult(BaseModel):
    motifs: List[dict]
    enriched_count: int
    method: str
    error: Optional[str] = None


RBP_MOTIFS = """MEME version 5.5.5

ALPHABET= ACGT

MOTIF RBM25
letter-probability matrix: alength= 4 w= 5 nsites= 20
 0.15  0.35  0.35  0.15
 0.35  0.15  0.15  0.35
 0.15  0.15  0.35  0.35
 0.35  0.35  0.15  0.15
 0.15  0.15  0.15  0.55

MOTIF PTB
letter-probability matrix: alength= 4 w= 6 nsites= 20
 0.25  0.25  0.25  0.25
 0.50  0.00  0.50  0.00
 0.50  0.00  0.50  0.00
 0.50  0.00  0.50  0.00
 0.50  0.00  0.50  0.00
 0.50  0.00  0.50  0.00

MOTIF QKI
letter-probability matrix: alength= 4 w= 7 nsites= 20
 0.14  0.14  0.14  0.58
 0.14  0.14  0.58  0.14
 0.14  0.58  0.14  0.14
 0.58  0.14  0.14  0.14
 0.14  0.14  0.58  0.14
 0.14  0.58  0.14  0.14
 0.58  0.14  0.14  0.14

MOTIF NOVA1
letter-probability matrix: alength= 4 w= 4 nsites= 20
 0.00  0.25  0.50  0.25
 0.50  0.00  0.25  0.25
 0.25  0.50  0.00  0.25
 0.25  0.25  0.50  0.00

MOTIF MBNL1
letter-probability matrix: alength= 4 w= 5 nsites= 20
 0.25  0.00  0.50  0.25
 0.00  0.25  0.00  0.75
 0.25  0.00  0.50  0.25
 0.00  0.25  0.00  0.75
 0.25  0.00  0.50  0.25

MOTIF CELF1
letter-probability matrix: alength= 4 w= 4 nsites= 20
 0.00  0.00  0.75  0.25
 0.00  0.75  0.00  0.25
 0.75  0.00  0.00  0.25
 0.00  0.75  0.00  0.25

MOTIF TIA1
letter-probability matrix: alength= 4 w= 5 nsites= 20
 0.00  0.00  0.50  0.50
 0.00  0.00  0.50  0.50
 0.00  0.00  0.50  0.50
 0.50  0.50  0.00  0.00
 0.50  0.50  0.00  0.00

MOTIF HuR
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33

MOTIF TDP43
letter-probability matrix: alength= 4 w= 5 nsites= 20
 0.00  0.40  0.20  0.40
 0.40  0.20  0.00  0.40
 0.40  0.20  0.00  0.40
 0.00  0.40  0.20  0.40
 0.40  0.20  0.00  0.40

MOTIF U2AF2
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.00  0.67  0.33  0.00
 0.00  0.67  0.33  0.00
 0.33  0.33  0.00  0.33

MOTIF SF3B4
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33

MOTIF HNRNPC
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.00  0.67  0.00  0.33
 0.33  0.00  0.67  0.00
 0.00  0.67  0.00  0.33

MOTIF EWSR1
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33
 0.00  0.67  0.00  0.33

MOTIF RBM10
letter-probability matrix: alength= 4 w= 4 nsites= 20
 0.00  0.50  0.25  0.25
 0.00  0.25  0.50  0.25
 0.00  0.25  0.25  0.50
 0.00  0.25  0.25  0.50

MOTIF FUS
letter-probability matrix: alength= 4 w= 4 nsites= 20
 0.00  0.00  0.67  0.33
 0.00  0.00  0.67  0.33
 0.67  0.00  0.00  0.33
 0.00  0.00  0.67  0.33

MOTIF hnRNPA1
letter-probability matrix: alength= 4 w= 5 nsites= 20
 0.00  0.40  0.20  0.40
 0.00  0.40  0.20  0.40
 0.00  0.40  0.20  0.40
 0.00  0.40  0.20  0.40
 0.00  0.40  0.20  0.40

MOTIF RPS14
letter-probability matrix: alength= 4 w= 3 nsites= 20
 0.33  0.00  0.33  0.33
 0.00  0.67  0.00  0.33
 0.33  0.00  0.33  0.33
"""


def write_fasta(sequences: List[dict], output_path: str) -> int:
    """Write sequences to FASTA file. Returns number of sequences written."""
    count = 0
    with open(output_path, "w") as f:
        for seq in sequences:
            if "sequence" in seq and seq["sequence"]:
                seq_id = seq.get("id", f"seq_{count}")
                f.write(f">{seq_id}\n{seq['sequence']}\n")
                count += 1
    return count


def write_motifs(output_path: str) -> str:
    """Write RBP motif database to MEME format file."""
    with open(output_path, "w") as f:
        f.write(RBP_MOTIFS)
    return output_path


def run_fimo(
    sequences_fa: str,
    motif_file: str,
    output_dir: str,
    pvalue: float = 0.05,
) -> dict:
    """Run FIMO (Find Individual Motif Occurrences)."""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = [
        "fimo",
        "--o",
        output_dir,
        "--thresh",
        str(pvalue),
        motif_file,
        sequences_fa,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            return {"error": result.stderr, "motifs": []}

        return parse_fimo_output(output_dir)

    except subprocess.TimeoutExpired:
        return {"error": "FIMO timed out", "motifs": []}
    except Exception as e:
        return {"error": str(e), "motifs": []}


def parse_fimo_output(output_dir: str) -> dict:
    """Parse FIMO output TSV file."""
    tsv_path = os.path.join(output_dir, "fimo.tsv")

    if not os.path.exists(tsv_path):
        for root, dirs, files in os.walk(output_dir):
            if "fimo.tsv" in files:
                tsv_path = os.path.join(root, "fimo.tsv")
                break

    motif_counts = {}
    if os.path.exists(tsv_path):
        with open(tsv_path, "r") as f:
            header = f.readline()  # Skip header
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 2:
                    motif_name = fields[0] if fields[0] else ""
                    if motif_name:
                        motif_counts[motif_name] = motif_counts.get(motif_name, 0) + 1

    motifs = [
        {"motif": name, "count": count, "source": "FIMO"}
        for name, count in sorted(
            motif_counts.items(), key=lambda x: x[1], reverse=True
        )
    ]

    return {"motifs": motifs, "enriched_count": len(motifs)}


@app.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "service": "meme"}


@app.post("/enrich", response_model=MotifEnrichmentResult)
def enrich(request: MotifEnrichmentRequest):
    """
    Run motif scanning analysis using FIMO.

    Args:
        request: MotifEnrichmentRequest with sequences and parameters

    Returns:
        MotifEnrichmentResult with enriched motifs
    """
    if not request.foreground_sequences:
        raise HTTPException(status_code=400, detail="No foreground sequences provided")

    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = os.path.join(tmpdir, "sequences.fa")
        motif_path = os.path.join(tmpdir, "motifs.meme")
        out_dir = os.path.join(tmpdir, "fimo_output")

        fg_count = write_fasta(request.foreground_sequences, fg_path)
        if fg_count == 0:
            raise HTTPException(
                status_code=400, detail="No valid sequences in foreground"
            )

        write_motifs(motif_path)

        result = run_fimo(
            fg_path,
            motif_path,
            out_dir,
            request.pvalue_threshold,
        )

        return MotifEnrichmentResult(
            motifs=result.get("motifs", []),
            enriched_count=result.get("enriched_count", 0),
            method="FIMO",
            error=result.get("error"),
        )

        if request.background_sequences:
            bg_count = write_fasta(request.background_sequences, bg_path)
            if bg_count == 0:
                bg_path = fg_path
        else:
            bg_path = fg_path

        write_motifs(motif_path)

        result = run_ame(
            fg_path, bg_path, motif_path, out_dir, request.pvalue_threshold
        )

        return MotifEnrichmentResult(
            motifs=result.get("motifs", []),
            enriched_count=result.get("enriched_count", 0),
            method="AME",
            error=result.get("error"),
        )


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)
