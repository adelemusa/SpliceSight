import React, { useState } from 'react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell } from 'recharts';
import { Info, TrendingUp, BookOpen, ChevronDown, ChevronUp, AlertCircle } from 'lucide-react';

const GSEAPanel = ({ enrichments = {} }) => {
  const [expanded, setExpanded] = useState(false);
  const [selectedGeneSet, setSelectedGeneSet] = useState('all');

  const pathways = Array.isArray(enrichments) ? enrichments : (enrichments.pathways || []);
  const geneCount = enrichments.gene_count || pathways.length;
  const error = enrichments.error;
  const totalResults = enrichments.total_results || pathways.length;

  const geneSets = ['all', ...new Set(pathways.map(e => e.Gene_set).filter(Boolean))];
  
  const filteredData = selectedGeneSet === 'all' 
    ? pathways 
    : pathways.filter(e => e.Gene_set === selectedGeneSet);

  const chartData = filteredData.slice(0, 20).map(enr => ({
    name: enr.Term?.split('(')[0]?.substring(0, 30) || 'Unknown',
    fullName: enr.Term || '',
    pValue: -Math.log10(enr["Adjusted P-value"] || enr["P-value"] || 1),
    adjPValue: enr["Adjusted P-value"] || 1,
    geneSet: enr.Gene_set || '',
    genes: enr.Genes || '',
    overlap: enr.Overlap || ''
  }));

  const getBarColor = (entry) => {
    const pVal = entry.adjPValue || 1;
    if (pVal < 0.001) return '#00cec9';
    if (pVal < 0.01) return '#6c5ce7';
    if (pVal < 0.05) return '#a29bfe';
    return '#74b9ff';
  };

  return (
    <div className="glass-panel" style={{ padding: '1.5rem', marginTop: '1rem' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
        <div>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <TrendingUp size={20} color="var(--accent)" />
            Functional Enrichment (GSEA)
          </h3>
          <p style={{ color: 'var(--text-muted)', fontSize: '0.85rem', marginTop: '0.5rem' }}>
            Gene Ontology and KEGG pathway enrichment analysis
            {geneCount > 0 && ` - ${geneCount} genes analyzed`}
          </p>
        </div>
        <button 
          onClick={() => setExpanded(!expanded)}
          style={{
            background: expanded ? "var(--accent)" : "rgba(108, 92, 231, 0.2)",
            border: "none",
            color: "white",
            padding: "0.5rem 1rem",
            borderRadius: "6px",
            cursor: "pointer",
            display: "flex",
            alignItems: "center",
            gap: "0.5rem",
            fontSize: "0.85rem"
          }}
        >
          {expanded ? <><ChevronUp size={16} /> Hide</> : <><ChevronDown size={16} /> Details</>}
        </button>
      </div>

      {/* Error Message */}
      {error && (
        <div style={{
          background: 'rgba(255, 107, 107, 0.1)',
          border: '1px solid rgba(255, 107, 107, 0.3)',
          borderRadius: '8px',
          padding: '1rem',
          marginBottom: '1rem',
          display: 'flex',
          alignItems: 'flex-start',
          gap: '0.75rem'
        }}>
          <AlertCircle size={20} color="#ff6b6b" style={{ flexShrink: 0, marginTop: '2px' }} />
          <div>
            <h4 style={{ margin: '0 0 0.25rem 0', color: '#ff6b6b', fontSize: '0.9rem' }}>
              Enrichment Analysis Error
            </h4>
            <p style={{ margin: 0, fontSize: '0.85rem', color: 'var(--text-muted)' }}>
              {error}
            </p>
          </div>
        </div>
      )}

      {/* Explanation */}
      <div style={{
        background: 'rgba(0, 206, 201, 0.1)',
        border: '1px solid rgba(0, 206, 201, 0.3)',
        borderRadius: '8px',
        padding: '1rem',
        marginBottom: '1rem'
      }}>
        <div style={{ display: 'flex', alignItems: 'flex-start', gap: '0.75rem' }}>
          <Info size={20} color="#00cec9" style={{ flexShrink: 0, marginTop: '2px' }} />
          <div>
            <h4 style={{ margin: '0 0 0.5rem 0', color: '#00cec9', fontSize: '0.95rem' }}>
              What is Gene Set Enrichment Analysis?
            </h4>
            <p style={{ margin: 0, fontSize: '0.85rem', color: 'var(--text-muted)', lineHeight: 1.5 }}>
              <strong>GSEA</strong> tests whether a predefined set of genes (e.g., genes in a KEGG pathway or GO biological process) 
              is statistically enriched in your list of differentially spliced genes. A pathway is "enriched" when it contains 
              more genes from your list than expected by chance. The X-axis shows <strong>-log10(FDR)</strong>: 
              higher values indicate greater statistical significance.
            </p>
          </div>
        </div>
      </div>

      {pathways.length > 0 ? (
        <>
          {/* Summary Stats */}
          <div style={{ 
            display: 'grid', 
            gridTemplateColumns: 'repeat(3, 1fr)', 
            gap: '1rem', 
            marginBottom: '1rem',
            fontSize: '0.85rem'
          }}>
            <div style={{ 
              background: 'rgba(108, 92, 231, 0.15)', 
              borderRadius: '8px', 
              padding: '0.75rem', 
              textAlign: 'center',
              border: '1px solid rgba(108, 92, 231, 0.3)'
            }}>
              <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#6c5ce7' }}>{geneCount}</div>
              <div style={{ color: 'var(--text-muted)' }}>Genes analyzed</div>
            </div>
            <div style={{ 
              background: 'rgba(0, 206, 201, 0.15)', 
              borderRadius: '8px', 
              padding: '0.75rem', 
              textAlign: 'center',
              border: '1px solid rgba(0, 206, 201, 0.3)'
            }}>
              <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#00cec9' }}>{totalResults}</div>
              <div style={{ color: 'var(--text-muted)' }}>Pathways found</div>
            </div>
            <div style={{ 
              background: 'rgba(253, 203, 110, 0.15)', 
              borderRadius: '8px', 
              padding: '0.75rem', 
              textAlign: 'center',
              border: '1px solid rgba(253, 203, 110, 0.3)'
            }}>
              <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#fdcb6e' }}>
                {pathways.filter(p => (p["Adjusted P-value"] || 1) < 0.05).length}
              </div>
              <div style={{ color: 'var(--text-muted)' }}>Significant (FDR &lt; 0.05)</div>
            </div>
          </div>

          {/* Gene Set Filter */}
          <div style={{ marginBottom: '1rem' }}>
            <label style={{ fontSize: '0.85rem', color: 'var(--text-muted)', marginRight: '0.5rem' }}>
              Filter by database:
            </label>
            <select
              value={selectedGeneSet}
              onChange={(e) => setSelectedGeneSet(e.target.value)}
              style={{
                background: 'rgba(0,0,0,0.3)',
                border: '1px solid var(--border-color)',
                color: 'white',
                padding: '0.4rem 0.8rem',
                borderRadius: '6px',
                fontSize: '0.85rem'
              }}
            >
              {geneSets.map(set => (
                <option key={set} value={set}>
                  {set === 'all' ? 'All databases' : set.replace(/_/g, ' ')}
                </option>
              ))}
            </select>
          </div>

          {/* Bar Chart */}
          <div style={{ height: chartData.length > 0 ? Math.min(450, chartData.length * 30 + 50) : 200 }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart 
                data={chartData} 
                layout="vertical"
                margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
              >
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  type="number" 
                  stroke="#888"
                  tickFormatter={(v) => v.toFixed(0)}
                  label={{ value: '-log10(FDR)', position: 'bottom', fill: '#888', fontSize: 12 }}
                />
                <YAxis 
                  type="category" 
                  dataKey="name" 
                  stroke="#888"
                  width={150}
                  tick={{ fontSize: 10 }}
                />
                <Tooltip 
                  contentStyle={{
                    background: 'rgba(30, 30, 50, 0.95)',
                    border: '1px solid rgba(108, 92, 231, 0.5)',
                    borderRadius: '8px',
                    fontSize: '0.8rem'
                  }}
                  formatter={(value, name, props) => {
                    const d = props.payload;
                    return [
                      <div key="tooltip" style={{ minWidth: 220 }}>
                        <div style={{ fontWeight: 'bold', marginBottom: 4, borderBottom: '1px solid rgba(255,255,255,0.1)', paddingBottom: 4 }}>
                          {d.fullName}
                        </div>
                        <div style={{ color: '#00cec9' }}>-log10(FDR): {value.toFixed(2)}</div>
                        <div style={{ color: '#a29bfe' }}>FDR: {d.adjPValue.toExponential(2)}</div>
                        <div style={{ color: '#fdcb6e' }}>Database: {d.geneSet}</div>
                        {d.overlap && <div style={{ color: '#74b9ff' }}>Overlap: {d.overlap}</div>}
                      </div>,
                      'Enrichment'
                    ];
                  }}
                />
                <Bar dataKey="pValue" radius={[0, 4, 4, 0]}>
                  {chartData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={getBarColor(entry)} />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>

          {/* Legend */}
          <div style={{ 
            display: 'flex', 
            gap: '1rem', 
            marginTop: '1rem', 
            flexWrap: 'wrap',
            fontSize: '0.8rem',
            color: 'var(--text-muted)'
          }}>
            <span style={{ display: 'flex', alignItems: 'center', gap: '0.3rem' }}>
              <span style={{ width: 12, height: 12, background: '#00cec9', borderRadius: 2 }}></span>
              FDR &lt; 0.001
            </span>
            <span style={{ display: 'flex', alignItems: 'center', gap: '0.3rem' }}>
              <span style={{ width: 12, height: 12, background: '#6c5ce7', borderRadius: 2 }}></span>
              FDR &lt; 0.01
            </span>
            <span style={{ display: 'flex', alignItems: 'center', gap: '0.3rem' }}>
              <span style={{ width: 12, height: 12, background: '#a29bfe', borderRadius: 2 }}></span>
              FDR &lt; 0.05
            </span>
            <span style={{ display: 'flex', alignItems: 'center', gap: '0.3rem' }}>
              <span style={{ width: 12, height: 12, background: '#74b9ff', borderRadius: 2 }}></span>
              FDR &gt; 0.05
            </span>
          </div>

          {/* Expanded List */}
          {expanded && (
            <div style={{ marginTop: '1.5rem', borderTop: '1px solid var(--border-color)', paddingTop: '1rem' }}>
              <h4 style={{ margin: '0 0 1rem 0', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <BookOpen size={16} color="var(--accent)" />
                All pathways ({filteredData.length})
              </h4>
              <div style={{ maxHeight: '400px', overflowY: 'auto' }}>
                {filteredData.map((enr, i) => (
                  <div 
                    key={i}
                    style={{
                      padding: '0.75rem',
                      marginBottom: '0.5rem',
                      background: 'rgba(108, 92, 231, 0.1)',
                      borderRadius: '6px',
                      borderLeft: `3px solid ${getBarColor({ adjPValue: enr["Adjusted P-value"] || 1 })}`
                    }}
                  >
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                      <div style={{ flex: 1 }}>
                        <div style={{ fontWeight: 'bold', fontSize: '0.9rem' }}>{enr.Term}</div>
                        <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)', marginTop: '0.25rem' }}>
                          {enr.Gene_set}
                        </div>
                        {enr.Genes && (
                          <div style={{ fontSize: '0.75rem', color: '#74b9ff', marginTop: '0.25rem', fontFamily: 'monospace' }}>
                            {String(enr.Genes).substring(0, 60)}{String(enr.Genes).length > 60 ? '...' : ''}
                          </div>
                        )}
                      </div>
                      <div style={{ textAlign: 'right', marginLeft: '1rem' }}>
                        <div style={{ 
                          fontWeight: 'bold', 
                          color: (enr["Adjusted P-value"] || 1) < 0.01 ? '#00cec9' : 'inherit'
                        }}>
                          FDR: {(enr["Adjusted P-value"] || 1).toExponential(2)}
                        </div>
                        {enr.Overlap && (
                          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)', marginTop: '0.25rem' }}>
                            {enr.Overlap}
                          </div>
                        )}
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </>
      ) : !error ? (
        <div style={{ 
          textAlign: 'center', 
          padding: '2rem',
          color: 'var(--text-muted)',
          background: 'rgba(0,0,0,0.2)',
          borderRadius: '8px'
        }}>
          <TrendingUp size={40} style={{ opacity: 0.3, marginBottom: '0.5rem' }} />
          <p style={{ margin: 0 }}>No pathways detected</p>
          <p style={{ fontSize: '0.85rem', marginTop: '0.5rem' }}>
            This may mean that the differentially spliced genes are not associated with known biological processes,
            or that the sample size is too small to detect significant enrichment.
          </p>
        </div>
      ) : null}
    </div>
  );
};

export default GSEAPanel;
