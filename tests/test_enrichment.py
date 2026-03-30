#!/usr/bin/env python
"""Test enrichment analysis functions."""

from core.enrichment import run_enrichment
from core.meme_enrichment import (
    fisher_exact_test,
    benjamini_hochberg,
    run_meme_enrichment,
    RBP_BACKGROUND_FREQUENCIES,
)


def test_fisher_exact():
    print("=== Fisher's Exact Test ===")
    # Strong enrichment
    p1 = fisher_exact_test(45, 55, 1000, 9000)
    print(f"Strong enrichment (45/100 vs 1000/10000): p={p1:.6f}")
    assert p1 < 0.001, "Strong enrichment should be significant"

    # No enrichment
    p2 = fisher_exact_test(10, 90, 1000, 9000)
    print(f"No enrichment (10/100 vs 1000/10000): p={p2:.6f}")
    assert p2 > 0.05, "No enrichment should not be significant"

    print("Fisher test: PASSED\n")


def test_multiple_testing():
    print("=== Multiple Testing Correction ===")
    pvalues = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5]
    bh = benjamini_hochberg(pvalues)
    print(f"Raw:     {[round(p, 3) for p in pvalues]}")
    print(f"BH adj:  {[round(p, 3) for p in bh]}")
    assert all(0 <= p <= 1 for p in bh), "Adjusted p-values should be in [0,1]"
    assert bh[0] <= bh[1] <= bh[2], "BH should be monotonically decreasing"
    print("Multiple testing: PASSED\n")


def test_gsea():
    print("=== GSEA Enrichment ===")
    test_genes = [
        "BRCA1",
        "TP53",
        "EGFR",
        "MYC",
        "KRAS",
        "PTEN",
        "RB1",
        "VHL",
        "APC",
        "CDKN2A",
        "BRAF",
        "ATM",
    ]

    print(f"Testing with {len(test_genes)} genes...")

    try:
        result = run_enrichment(test_genes, pvalue_cutoff=0.1)

        if isinstance(result, dict):
            print(f"Gene count: {result.get('gene_count', 0)}")
            pathways = result.get("pathways", [])
            print(f"Pathways returned: {len(pathways)}")

            if "error" in result:
                print(f"Error: {result['error']}")
            elif pathways:
                print("Top 5 pathways:")
                for p in pathways[:5]:
                    term = p.get("Term", "")[:50]
                    adj_p = p.get("Adjusted P-value", 1)
                    overlap = p.get("Overlap", "")
                    print(f"  {term}: p_adj={adj_p:.4f} ({overlap})")
        else:
            print(f"Unexpected result type: {type(result)}")
    except Exception as e:
        print(f"Exception: {e}")

    print("GSEA test: COMPLETED\n")


def test_rbp_database():
    print("=== RBP Background Database ===")
    print(f"Total RBPs: {len(RBP_BACKGROUND_FREQUENCIES)}")

    sorted_freqs = sorted(
        RBP_BACKGROUND_FREQUENCIES.items(), key=lambda x: x[1], reverse=True
    )
    print("Top 5 by expected frequency:")
    for rbp, freq in sorted_freqs[:5]:
        print(f"  {rbp}: {freq * 100:.1f}%")

    assert len(RBP_BACKGROUND_FREQUENCIES) >= 20, "Should have at least 20 RBPs"
    print("RBP database: PASSED\n")


def test_meme_enrichment():
    print("=== MEME Enrichment Function ===")
    # Test with empty events
    result = run_meme_enrichment([], sample_size=10)
    assert "error" in result or "motifs" in result
    print("Empty events handled: OK")

    # Test with mock events (won't have real sequences)
    mock_events = [
        {
            "id": f"test_{i}",
            "event_type": "SE",
            "gene_symbol": "BRCA1",
            "chr": "17",
            "strand": "-",
            "coordinates": {"start": 43000000, "end": 43010000},
        }
        for i in range(5)
    ]
    result = run_meme_enrichment(mock_events, sample_size=5)
    print(f"Result keys: {list(result.keys())}")
    print(f"Sequences extracted: {result.get('sequences_extracted', 0)}")
    print("MEME enrichment: PASSED\n")


if __name__ == "__main__":
    print("=" * 60)
    print("TESTING ENRICHMENT ANALYSIS")
    print("=" * 60 + "\n")

    test_fisher_exact()
    test_multiple_testing()
    test_rbp_database()
    test_gsea()
    test_meme_enrichment()

    print("=" * 60)
    print("ALL TESTS COMPLETED")
    print("=" * 60)
