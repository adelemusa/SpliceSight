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


def run_ame(
    foreground_fa: str,
    background_fa: str,
    motif_file: str,
    output_dir: str,
    pvalue: float = 0.05,
) -> dict:
    """Run AME (Analysis of Motif Enrichment)."""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    cmd = [
        "ame",
        "--o",
        output_dir,
        "--motif",
        motif_file,
        "--bgfile",
        background_fa,
        "--control",
        "--method",
        "fisher",
        "--pvalue",
        str(pvalue),
        "--hit-lof-fraction",
        "0.02",
        "--max-stored-scores",
        "100000",
        foreground_fa,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            return {"error": result.stderr, "motifs": [], "enriched_count": 0}

        return parse_ame_output(output_dir)

    except subprocess.TimeoutExpired:
        return {"error": "AME timed out", "motifs": [], "enriched_count": 0}
    except Exception as e:
        return {"error": str(e), "motifs": [], "enriched_count": 0}


def parse_ame_output(output_dir: str) -> dict:
    """Parse AME output (HTML or TSV)."""
    ame_tsv = os.path.join(output_dir, "ame.tsv")
    ame_html = os.path.join(output_dir, "ame.html")

    motifs = []

    if os.path.exists(ame_tsv):
        with open(ame_tsv, "r") as f:
            header = f.readline().strip().split("\t")
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 5:
                    motif_name = fields[1] if fields[1] else ""
                    pvalue = fields[3] if len(fields) > 3 else "1"
                    evalue = fields[4] if len(fields) > 4 else "1"

                    try:
                        motifs.append(
                            {
                                "motif": motif_name,
                                "pvalue": float(pvalue),
                                "evalue": float(evalue),
                                "source": "AME",
                            }
                        )
                    except ValueError:
                        continue

    motifs.sort(key=lambda x: x.get("pvalue", 1))

    return {"motifs": motifs, "enriched_count": len(motifs), "method": "AME"}


def create_background_from_foreground(sequences: List[dict], output_path: str) -> int:
    """Create a shuffled background from foreground sequences."""
    import random

    count = 0
    with open(output_path, "w") as f:
        for seq in sequences:
            if "sequence" in seq and seq["sequence"]:
                seq_id = seq.get("id", f"bg_{count}")
                sequence = seq["sequence"]
                shuffled = list(sequence)
                random.shuffle(shuffled)
                f.write(f">{seq_id}_shuffled\n{''.join(shuffled)}\n")
                count += 1
    return count


@app.get("/health")
def health():
    """Health check endpoint."""
    fimo_path = shutil.which("fimo")
    meme_installed = fimo_path is not None

    return {
        "status": "ok" if meme_installed else "degraded",
        "service": "meme",
        "fimo_available": meme_installed,
        "fimo_path": fimo_path,
        "version": "5.5.5" if meme_installed else None,
    }


@app.post("/enrich", response_model=MotifEnrichmentResult)
def enrich(request: MotifEnrichmentRequest):
    """
    Run motif enrichment analysis using AME (Analysis of Motif Enrichment).
    AME compares foreground sequences against background to find statistically enriched motifs.

    Args:
        request: MotifEnrichmentRequest with sequences and parameters

    Returns:
        MotifEnrichmentResult with enriched motifs and p-values
    """
    if not request.foreground_sequences:
        raise HTTPException(status_code=400, detail="No foreground sequences provided")

    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = os.path.join(tmpdir, "foreground.fa")
        bg_path = os.path.join(tmpdir, "background.fa")
        motif_path = os.path.join(tmpdir, "motifs.meme")
        out_dir = os.path.join(tmpdir, "ame_output")

        fg_count = write_fasta(request.foreground_sequences, fg_path)
        if fg_count == 0:
            raise HTTPException(
                status_code=400, detail="No valid sequences in foreground"
            )

        write_motifs(motif_path)

        if request.background_sequences:
            bg_count = write_fasta(request.background_sequences, bg_path)
            if bg_count == 0:
                create_background_from_foreground(request.foreground_sequences, bg_path)
        else:
            create_background_from_foreground(request.foreground_sequences, bg_path)

        result = run_ame(
            fg_path,
            bg_path,
            motif_path,
            out_dir,
            request.pvalue_threshold,
        )

        return MotifEnrichmentResult(
            motifs=result.get("motifs", []),
            enriched_count=result.get("enriched_count", 0),
            method=result.get("method", "AME"),
            error=result.get("error"),
        )


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)
