import gseapy as gp
import time


def run_enrichment(
    gene_symbols,
    gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2016"],
    pvalue_cutoff=0.1,
    max_retries=3,
):
    if not gene_symbols or len(gene_symbols) < 3:
        return {
            "error": "Not enough genes for enrichment analysis (minimum 3 required)",
            "pathways": [],
            "gene_count": len(gene_symbols) if gene_symbols else 0,
        }

    last_error = None

    for attempt in range(max_retries):
        try:
            enr = gp.enrichr(
                gene_list=gene_symbols,
                gene_sets=gene_sets,
                organism="human",
                outdir=None,
                cutoff=pvalue_cutoff,
            )
            if enr.results is not None and not enr.results.empty:
                results = enr.results[
                    [
                        "Term",
                        "Adjusted P-value",
                        "P-value",
                        "Gene_set",
                        "Overlap",
                        "Genes",
                    ]
                ].copy()
                results = results.sort_values("Adjusted P-value")
                return {
                    "pathways": results.head(30).to_dict("records"),
                    "total_results": len(results),
                    "gene_count": len(gene_symbols),
                }
            return {"pathways": [], "gene_count": len(gene_symbols)}
        except Exception as e:
            last_error = str(e)
            error_str = str(e).lower()

            if (
                "503" in error_str
                or "unavailable" in error_str
                or "timeout" in error_str
            ):
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 2
                    time.sleep(wait_time)
                    continue
                last_error = f"Enrichr service temporarily unavailable (503). Please try again in a few minutes. Error: {last_error}"
            else:
                break

    return {
        "error": f"Enrichr API error: {last_error}",
        "pathways": [],
        "gene_count": len(gene_symbols),
    }
