import gseapy as gp


def run_enrichment(
    gene_symbols,
    gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"],
    pvalue_cutoff=0.1,
):
    if not gene_symbols or len(gene_symbols) < 3:
        return {
            "error": "Not enough genes for enrichment analysis (minimum 3 required)",
            "pathways": [],
            "gene_count": len(gene_symbols) if gene_symbols else 0,
        }

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
                ["Term", "Adjusted P-value", "P-value", "Gene_set", "Overlap", "Genes"]
            ].copy()
            results = results.sort_values("Adjusted P-value")
            return {
                "pathways": results.head(20).to_dict("records"),
                "total_results": len(results),
                "gene_count": len(gene_symbols),
            }
        return {"pathways": [], "gene_count": len(gene_symbols)}
    except Exception as e:
        return {
            "error": f"Enrichr API error: {str(e)}",
            "pathways": [],
            "gene_count": len(gene_symbols),
        }
