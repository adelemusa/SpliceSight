import gseapy as gp

def run_enrichment(gene_symbols, gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"]):
    if not gene_symbols:
        return []
        
    try:
        enr = gp.enrichr(gene_list=gene_symbols, 
                         gene_sets=gene_sets, 
                         organism='human',
                         outdir=None) 
        if enr.results is not None and not enr.results.empty:
            sig_results = enr.results[enr.results["Adjusted P-value"] < 0.05]
            # Convert directly for frontend mapping
            if not sig_results.empty:
                return sig_results[["Term", "Adjusted P-value", "Gene_set"]].head(10).to_dict('records')
    except Exception as e:
        print(f"Error finding GSEA pathways: {str(e)}")
        
    return []
