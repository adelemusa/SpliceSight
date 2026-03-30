import React, { useState } from 'react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell } from 'recharts';
import { Info, Dna, Search, Zap, ChevronDown, ChevronUp, Star } from 'lucide-react';

const MotifPanel = ({ motifData }) => {
  const [expanded, setExpanded] = useState(false);
  const [selectedEventType, setSelectedEventType] = useState('all');
  const [showSignificantOnly, setShowSignificantOnly] = useState(false);

  if (!motifData) {
    return (
      <div className="widget glass-panel">
        <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
          <Zap size={18} color="var(--accent)" />
          RBP Motif Analysis
        </h3>
        <p style={{ color: 'var(--text-muted)', marginTop: '1rem', textAlign: 'center' }}>
          Data not available. Run the analysis to view RBP motifs.
        </p>
      </div>
    );
  }

  const enrichedMotifs = motifData.enriched_motifs || [];
  const byEventType = motifData.by_event_type || {};
  const topGenes = motifData.top_genes_with_motifs || [];
  const totalAnalyzed = motifData.total_events_analyzed || 0;
  const significantMotifs = motifData.significant_motifs || [];
  const method = motifData.method || 'Motif scanning';

  const filteredMotifs = showSignificantOnly 
    ? enrichedMotifs.filter(m => m.significant)
    : enrichedMotifs;

  const chartData = filteredMotifs.slice(0, 15).map(m => ({
    name: m.rbp,
    count: m.count,
    percentage: m.percentage,
    pvalue: m.pvalue || 1,
    adjPvalue: m.adj_pvalue || 1,
    foldEnrichment: m.fold_enrichment || 0,
    significant: m.significant || false,
    significance: m.significance || 'ns'
  }));

  const getBarColor = (motif) => {
    if (motif.significant) return '#00cec9';
    if (motif.adjPvalue < 0.1) return '#6c5ce7';
    if (motif.adjPvalue < 0.25) return '#a29bfe';
    return '#74b9ff';
  };

  const getSignificanceColor = (sig) => {
    if (sig === '***') return '#00cec9';
    if (sig === '**') return '#6c5ce7';
    if (sig === '*') return '#fdcb6e';
    return '#888';
  };

  const eventTypes = Object.keys(byEventType);

  const renderEventTypePatterns = (eventType) => {
    const patterns = byEventType[eventType]?.patterns || [];
    return patterns.map((p, i) => (
      <div key={i} style={{ 
        background: 'rgba(0,0,0,0.2)', 
        padding: '0.5rem', 
        borderRadius: '6px', 
        marginBottom: '0.5rem',
        fontSize: '0.85rem'
      }}>
        <div style={{ fontWeight: 'bold', color: '#00cec9', marginBottom: '0.25rem' }}>{p.name}</div>
        <div style={{ fontFamily: 'monospace', color: '#fdcb6e', marginBottom: '0.25rem' }}>{p.motif}</div>
        <div style={{ color: 'var(--text-muted)', fontSize: '0.8rem' }}>{p.frequency}</div>
      </div>
    ));
  };

  return (
    <div className="widget glass-panel">
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
        <div>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Zap size={18} color="var(--accent)" />
            RBP Motif Analysis
          </h3>
          <p style={{ color: 'var(--text-muted)', fontSize: '0.8rem', marginTop: '0.5rem' }}>
            {method} - Statistical enrichment analysis
          </p>
        </div>
        <button
          onClick={() => setExpanded(!expanded)}
          style={{
            background: expanded ? "var(--accent)" : "rgba(108, 92, 231, 0.2)",
            border: "none",
            color: "white",
            padding: "0.4rem 0.8rem",
            borderRadius: "6px",
            cursor: "pointer",
            fontSize: "0.8rem",
            display: "flex",
            alignItems: "center",
            gap: "0.3rem"
          }}
        >
          {expanded ? <><ChevronUp size={14} /> Hide</> : <><ChevronDown size={14} /> Details</>}
        </button>
      </div>

      {/* Explanation */}
      {expanded && (
        <div style={{
          background: 'rgba(108, 92, 231, 0.1)',
          border: '1px solid rgba(108, 92, 231, 0.3)',
          borderRadius: '8px',
          padding: '1rem',
          marginBottom: '1rem'
        }}>
          <h4 style={{ margin: '0 0 0.5rem 0', color: '#6c5ce7', fontSize: '0.9rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Dna size={16} /> What are RBP motifs?
          </h4>
          <div style={{ fontSize: '0.85rem', color: 'var(--text-muted)', lineHeight: 1.6 }}>
            <p style={{ margin: '0 0 0.5rem 0' }}>
              <strong>RNA-Binding Proteins (RBPs)</strong> are proteins that recognize specific sequences 
              (motifs) in RNA. These motifs regulate splicing, mRNA stability, localization, and translation.
            </p>
            <p style={{ margin: 0 }}>
              This analysis uses <strong>Fisher's exact test</strong> to identify RBPs whose motifs are 
              statistically enriched in your splicing events compared to genome background frequencies.
            </p>
          </div>
        </div>
      )}

      {/* Statistics */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '1rem', marginBottom: '1.5rem' }}>
        <div style={{ background: 'rgba(108, 92, 231, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(108, 92, 231, 0.3)' }}>
          <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#6c5ce7' }}>{totalAnalyzed}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Events analyzed</div>
        </div>
        <div style={{ background: 'rgba(0, 206, 201, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(0, 206, 201, 0.3)' }}>
          <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#00cec9' }}>{significantMotifs.length}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Significant RBPs</div>
        </div>
        <div style={{ background: 'rgba(253, 203, 110, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(253, 203, 110, 0.3)' }}>
          <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#fdcb6e' }}>{enrichedMotifs.length}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Total motifs</div>
        </div>
        <div style={{ background: 'rgba(116, 185, 255, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(116, 185, 255, 0.3)' }}>
          <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#74b9ff' }}>{motifData.sequences_extracted || 0}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Sequences</div>
        </div>
      </div>

      {/* Chart */}
      {chartData.length > 0 ? (
        <div style={{ marginBottom: '1.5rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.75rem' }}>
            <h4 style={{ margin: 0, fontSize: '0.9rem', color: 'var(--text-muted)' }}>
              <Search size={14} style={{ marginRight: '0.5rem', verticalAlign: 'middle' }} />
              Enriched RBP motifs ({filteredMotifs.length} {showSignificantOnly ? 'significant' : 'total'})
            </h4>
            <label style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', fontSize: '0.8rem', color: 'var(--text-muted)', cursor: 'pointer' }}>
              <input 
                type="checkbox" 
                checked={showSignificantOnly}
                onChange={(e) => setShowSignificantOnly(e.target.checked)}
                style={{ accentColor: '#00cec9' }}
              />
              Significant only
            </label>
          </div>
          <div style={{ height: chartData.length > 10 ? 350 : 200 }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={chartData} layout="vertical" margin={{ top: 5, right: 30, left: 100, bottom: 5 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis type="number" stroke="#888" tickFormatter={(v) => `${v}%`} />
                <YAxis type="category" dataKey="name" stroke="#888" width={95} tick={{ fontSize: 11 }} />
                <Tooltip 
                  contentStyle={{ background: 'rgba(30, 30, 50, 0.95)', border: '1px solid rgba(108, 92, 231, 0.5)', borderRadius: '8px', fontSize: '0.8rem' }}
                  formatter={(value, name, props) => [
                    `${value}%`,
                    'Frequency'
                  ]}
                  labelFormatter={(label, payload) => {
                    if (payload && payload[0]) {
                      const d = payload[0].payload;
                      return (
                        <div style={{ minWidth: 200 }}>
                          <div style={{ fontWeight: 'bold', marginBottom: 4 }}>{d.name}</div>
                          <div>Count: {d.count} sequences</div>
                          <div>Frequency: {d.percentage}%</div>
                          <div>Fold Enrichment: {d.foldEnrichment}x</div>
                          <div style={{ color: getSignificanceColor(d.significance) }}>
                            {d.significance !== 'ns' ? `Significance: ${d.significance} (p_adj=${d.adjPvalue.toExponential(2)})` : `p_adj=${d.adjPvalue.toExponential(2)}`}
                          </div>
                        </div>
                      );
                    }
                    return label;
                  }}
                />
                <Bar dataKey="percentage" radius={[0, 4, 4, 0]}>
                  {chartData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={getBarColor(entry)} />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>
        </div>
      ) : (
        <div style={{ textAlign: 'center', padding: '2rem', color: 'var(--text-muted)', background: 'rgba(0,0,0,0.2)', borderRadius: '8px' }}>
          <Zap size={32} style={{ opacity: 0.3, marginBottom: '0.5rem' }} />
          <p style={{ margin: 0 }}>No{showSignificantOnly ? ' significant' : ''} motifs found</p>
        </div>
      )}

      {/* Significance Legend */}
      {expanded && (
        <div style={{ 
          display: 'flex', 
          gap: '1rem', 
          marginBottom: '1rem',
          flexWrap: 'wrap',
          fontSize: '0.8rem',
          color: 'var(--text-muted)',
          padding: '0.75rem',
          background: 'rgba(0,0,0,0.2)',
          borderRadius: '6px'
        }}>
          <span style={{ fontWeight: 'bold', marginRight: '0.5rem' }}>Significance:</span>
          <span style={{ color: '#00cec9' }}>*** p_adj &lt; 0.001</span>
          <span style={{ color: '#6c5ce7' }}>** p_adj &lt; 0.01</span>
          <span style={{ color: '#fdcb6e' }}>* p_adj &lt; 0.05</span>
          <span style={{ color: '#888' }}>ns = not significant</span>
        </div>
      )}

      {/* Details by event type */}
      {expanded && eventTypes.length > 0 && (
        <div style={{ marginTop: '1rem', borderTop: '1px solid var(--border-color)', paddingTop: '1rem' }}>
          <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: 'var(--text-muted)' }}>
            Motifs by event type
          </h4>
          <select
            value={selectedEventType}
            onChange={(e) => setSelectedEventType(e.target.value)}
            style={{
              width: '100%',
              background: 'rgba(0,0,0,0.3)',
              border: '1px solid var(--border-color)',
              color: 'white',
              padding: '0.5rem',
              borderRadius: '6px',
              fontSize: '0.85rem',
              marginBottom: '1rem'
            }}
          >
            <option value="all">All events</option>
            {eventTypes.map(type => (
              <option key={type} value={type}>{type} ({byEventType[type]?.count || 0} events)</option>
            ))}
          </select>

          {(selectedEventType === 'all' ? eventTypes : [selectedEventType]).map(type => {
            const typeData = byEventType[type];
            if (!typeData) return null;
            
            return (
              <div key={type} style={{ marginBottom: '1rem' }}>
                <h5 style={{ margin: '0 0 0.5rem 0', fontSize: '0.9rem', color: '#6c5ce7' }}>
                  {type} - Splicing patterns
                </h5>
                {renderEventTypePatterns(type)}
                
                {typeData.motifs && Object.keys(typeData.motifs).length > 0 && (
                  <div style={{ marginTop: '0.5rem' }}>
                    <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)', marginBottom: '0.25rem' }}>Associated RBP motifs:</div>
                    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem' }}>
                      {Object.entries(typeData.motifs).sort((a, b) => b[1] - a[1]).slice(0, 8).map(([rbp, count]) => (
                        <span key={rbp} style={{
                          background: 'rgba(0, 206, 201, 0.2)',
                          padding: '2px 8px',
                          borderRadius: '12px',
                          fontSize: '0.75rem',
                          color: '#00cec9'
                        }}>
                          {rbp} ({count})
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            );
          })}
        </div>
      )}

      {/* Top genes */}
      {expanded && topGenes.length > 0 && (
        <div style={{ marginTop: '1rem', borderTop: '1px solid var(--border-color)', paddingTop: '1rem' }}>
          <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: 'var(--text-muted)' }}>
            Top 10 genes with most motifs
          </h4>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(200px, 1fr))', gap: '0.5rem' }}>
            {topGenes.map((gene, i) => (
              <div key={i} style={{ 
                background: 'rgba(0,0,0,0.2)', 
                padding: '0.75rem', 
                borderRadius: '6px',
                borderLeft: `3px solid ${i < 3 ? '#fdcb6e' : 'rgba(108, 92, 231, 0.3)'}`
              }}>
                <div style={{ fontWeight: 'bold', marginBottom: '0.25rem' }}>{gene.gene}</div>
                <div style={{ fontSize: '0.75rem', color: 'var(--text-muted)' }}>
                  {gene.motif_count} motifs: {gene.rbps?.join(', ')}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Legend */}
      {expanded && enrichedMotifs.length > 0 && (
        <div style={{ marginTop: '1rem', borderTop: '1px solid var(--border-color)', paddingTop: '1rem' }}>
          <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: 'var(--text-muted)' }}>
            RBP reference database
          </h4>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)', lineHeight: 1.8 }}>
            <p style={{ margin: '0 0 0.5rem 0' }}>
              Motifs were identified using CLIP-seq and iCLIP data from published studies:
            </p>
            <ul style={{ margin: 0, paddingLeft: '1.5rem' }}>
              <li>Zhou et al. 2018 - RBM25, QKI binding sites</li>
              <li>Wang et al. 2012 - MBNL1/2 expanD motifs</li>
              <li>Yoon et al. 2016 - PTB/hnRNPs binding</li>
              <li>Ule et al. 2005 - NOVA1/2 neuronal splicing</li>
              <li>Kishore et al. 2011 - HuR/ELAVL1 targets</li>
            </ul>
          </div>
        </div>
      )}
    </div>
  );
};

export default MotifPanel;
