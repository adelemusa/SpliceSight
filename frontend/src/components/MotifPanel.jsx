import React, { useState } from 'react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell, LineChart, Line, Legend } from 'recharts';
import { Info, Dna, Search, Zap, ChevronDown, ChevronUp, Star, Activity, TrendingUp } from 'lucide-react';

const MotifPanel = ({ motifData, rmapsData }) => {
  const [expanded, setExpanded] = useState(false);
  const [selectedEventType, setSelectedEventType] = useState('all');
  const [showSignificantOnly, setShowSignificantOnly] = useState(false);
  const [selectedRBP, setSelectedRBP] = useState(null);
  const [activeTab, setActiveTab] = useState('summary');

  if (!motifData && !rmapsData) {
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

  const enrichedMotifs = motifData?.enriched_motifs || [];
  const byEventType = motifData?.by_event_type || {};
  const topGenes = motifData?.top_genes_with_motifs || [];
  const totalAnalyzed = motifData?.total_events_analyzed || 0;
  const significantMotifs = motifData?.significant_motifs || [];
  const method = motifData?.method || 'Motif scanning';

  const rmapsMotifs = rmapsData?.motifs || [];
  const rmapsMethod = rmapsData?.method || 'rMAPS-style analysis';
  const eventCounts = rmapsData?.event_counts || {};
  const diagnostics = rmapsData?.diagnostics || {};
  const rnaMap = rmapsData?.rna_map || {};

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

  const rmapsChartData = rmapsMotifs.slice(0, 20).map((m, i) => ({
    rank: i + 1,
    rbp: m.rbp,
    region: m.region || '',
    direction: m.direction || 'upregulated',
    density: m.avg_density_up || m.avg_density_down || 0,
    bgDensity: m.avg_density_bg || 0,
    foldChange: m.fold_change || 0,
    pvalue: m.pvalue || 1,
    adjPvalue: m.adj_pvalue || 1,
  }));

  const getRNAMapData = (rbpName) => {
    if (!rnaMap[rbpName]) return null;
    
    const regions = ['upstream_exon', 'upstream_intron', 'target_exon', 'downstream_intron', 'downstream_exon'];
    const regionLabels = ['Upstream Exon', 'Upstream Intron', 'Target Exon', 'Downstream Intron', 'Downstream Exon'];
    
    return regions.map((region, i) => {
      const regionData = rnaMap[rbpName]?.[region] || {};
      const upScores = regionData.up || [];
      const downScores = regionData.down || [];
      const bgScores = regionData.bg || [];
      
      return {
        name: regionLabels[i],
        shortName: region.split('_').pop(),
        upregulated: upScores.length > 0 ? upScores.reduce((a, b) => a + b, 0) / upScores.length : 0,
        downregulated: downScores.length > 0 ? downScores.reduce((a, b) => a + b, 0) / downScores.length : 0,
        background: bgScores.length > 0 ? bgScores.reduce((a, b) => a + b, 0) / bgScores.length : 0,
        nUp: upScores.length,
        nDown: downScores.length,
        nBg: bgScores.length,
      };
    });
  };

  const rnaMapChartData = selectedRBP ? getRNAMapData(selectedRBP) : null;

  return (
    <div className="widget glass-panel">
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
        <div>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Zap size={18} color="var(--accent)" />
            RBP Motif Analysis
          </h3>
          <p style={{ color: 'var(--text-muted)', fontSize: '0.8rem', marginTop: '0.5rem' }}>
            {rmapsData ? rmapsMethod : method} - Statistical enrichment analysis
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

      {/* Tab Navigation */}
      {rmapsData && (
        <div style={{ display: 'flex', gap: '0.5rem', marginBottom: '1rem', borderBottom: '1px solid rgba(255,255,255,0.1)', paddingBottom: '0.5rem' }}>
          <button
            onClick={() => setActiveTab('summary')}
            style={{
              background: activeTab === 'summary' ? 'rgba(0, 206, 201, 0.2)' : 'transparent',
              border: 'none',
              color: activeTab === 'summary' ? '#00cec9' : 'var(--text-muted)',
              padding: '0.5rem 1rem',
              borderRadius: '6px',
              cursor: 'pointer',
              fontSize: '0.85rem',
              fontWeight: activeTab === 'summary' ? 'bold' : 'normal',
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem'
            }}
          >
            <Activity size={14} /> rMAPS Summary
          </button>
          <button
            onClick={() => setActiveTab('rnamap')}
            style={{
              background: activeTab === 'rnamap' ? 'rgba(0, 206, 201, 0.2)' : 'transparent',
              border: 'none',
              color: activeTab === 'rnamap' ? '#00cec9' : 'var(--text-muted)',
              padding: '0.5rem 1rem',
              borderRadius: '6px',
              cursor: 'pointer',
              fontSize: '0.85rem',
              fontWeight: activeTab === 'rnamap' ? 'bold' : 'normal',
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem'
            }}
          >
            <TrendingUp size={14} /> RNA Map
          </button>
        </div>
      )}

      {/* rMAPS Summary Tab */}
      {activeTab === 'summary' && rmapsData && (
        <>
          {/* Event Counts */}
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '1rem', marginBottom: '1.5rem' }}>
            <div style={{ background: 'rgba(255, 107, 107, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(255, 107, 107, 0.3)' }}>
              <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#ff6b6b' }}>{eventCounts.upregulated || 0}</div>
              <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Upregulated</div>
            </div>
            <div style={{ background: 'rgba(116, 185, 255, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(116, 185, 255, 0.3)' }}>
              <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#74b9ff' }}>{eventCounts.downregulated || 0}</div>
              <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Downregulated</div>
            </div>
            <div style={{ background: 'rgba(46, 213, 115, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(46, 213, 115, 0.3)' }}>
              <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#2ed573' }}>{eventCounts.background || 0}</div>
              <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Background</div>
            </div>
            <div style={{ background: 'rgba(253, 203, 110, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(253, 203, 110, 0.3)' }}>
              <div style={{ fontSize: '1.8rem', fontWeight: 'bold', color: '#fdcb6e' }}>{rmapsMotifs.length}</div>
              <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Enriched RBPs</div>
            </div>
          </div>

          {/* rMAPS Motif Enrichment Table */}
          {rmapsChartData.length > 0 ? (
            <div style={{ marginBottom: '1.5rem' }}>
              <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: 'var(--text-muted)', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <Search size={14} /> Enriched RBP Motifs (rMAPS)
              </h4>
              <div style={{ maxHeight: 300, overflowY: 'auto', borderRadius: '8px', border: '1px solid rgba(255,255,255,0.1)' }}>
                <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.8rem' }}>
                  <thead style={{ position: 'sticky', top: 0, background: 'rgba(30, 30, 50, 0.95)' }}>
                    <tr style={{ borderBottom: '1px solid rgba(255,255,255,0.1)' }}>
                      <th style={{ padding: '0.75rem', textAlign: 'left', color: '#888' }}>#</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', color: '#888' }}>RBP</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', color: '#888' }}>Region</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', color: '#888' }}>Direction</th>
                      <th style={{ padding: '0.75rem', textAlign: 'right', color: '#888' }}>Density %</th>
                      <th style={{ padding: '0.75rem', textAlign: 'right', color: '#888' }}>Fold</th>
                      <th style={{ padding: '0.75rem', textAlign: 'right', color: '#888' }}>p-value</th>
                      <th style={{ padding: '0.75rem', textAlign: 'right', color: '#888' }}>p-adj</th>
                    </tr>
                  </thead>
                  <tbody>
                    {rmapsChartData.map((row, i) => (
                      <tr 
                        key={i} 
                        style={{ 
                          borderBottom: '1px solid rgba(255,255,255,0.05)',
                          cursor: 'pointer',
                          background: selectedRBP === row.rbp ? 'rgba(0, 206, 201, 0.1)' : 'transparent'
                        }}
                        onClick={() => {
                          setSelectedRBP(row.rbp);
                          setActiveTab('rnamap');
                        }}
                      >
                        <td style={{ padding: '0.5rem 0.75rem', color: '#888' }}>{row.rank}</td>
                        <td style={{ padding: '0.5rem 0.75rem', fontWeight: 'bold', color: '#00cec9' }}>{row.rbp}</td>
                        <td style={{ padding: '0.5rem 0.75rem', color: '#a29bfe' }}>{row.region.replace('_', ' ')}</td>
                        <td style={{ padding: '0.5rem 0.75rem' }}>
                          <span style={{
                            color: row.direction === 'upregulated' ? '#ff6b6b' : '#74b9ff',
                            fontWeight: 'bold'
                          }}>
                            {row.direction === 'upregulated' ? '↑ UP' : '↓ DOWN'}
                          </span>
                        </td>
                        <td style={{ padding: '0.5rem 0.75rem', textAlign: 'right', color: '#fff' }}>{row.density.toFixed(2)}</td>
                        <td style={{ padding: '0.5rem 0.75rem', textAlign: 'right', color: row.foldChange > 1 ? '#2ed573' : '#ff6b6b' }}>
                          {row.foldChange.toFixed(2)}x
                        </td>
                        <td style={{ padding: '0.5rem 0.75rem', textAlign: 'right', color: row.pvalue < 0.05 ? '#fdcb6e' : '#888' }}>
                          {row.pvalue < 0.001 ? row.pvalue.toExponential(2) : row.pvalue.toFixed(4)}
                        </td>
                        <td style={{ padding: '0.5rem 0.75rem', textAlign: 'right', color: row.adjPvalue < 0.05 ? '#00cec9' : '#888' }}>
                          {row.adjPvalue < 0.001 ? row.adjPvalue.toExponential(2) : row.adjPvalue.toFixed(4)}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          ) : (
            <div style={{ textAlign: 'center', padding: '2rem', color: 'var(--text-muted)', background: 'rgba(0,0,0,0.2)', borderRadius: '8px', marginBottom: '1.5rem' }}>
              <Activity size={32} style={{ opacity: 0.3, marginBottom: '0.5rem' }} />
              <p style={{ margin: 0 }}>No significant RBP motif enrichment detected</p>
            </div>
          )}

          {/* Region Legend */}
          <div style={{ 
            display: 'flex', 
            gap: '1rem', 
            flexWrap: 'wrap',
            fontSize: '0.8rem',
            color: 'var(--text-muted)',
            padding: '0.75rem',
            background: 'rgba(0,0,0,0.2)',
            borderRadius: '6px',
            marginBottom: '1rem'
          }}>
            <span style={{ fontWeight: 'bold', marginRight: '0.5rem' }}>Regions:</span>
            <span style={{ color: '#ff6b6b' }}>upstream_exon</span>
            <span style={{ color: '#a29bfe' }}>upstream_intron</span>
            <span style={{ color: '#00cec9' }}>target_exon</span>
            <span style={{ color: '#fdcb6e' }}>downstream_intron</span>
            <span style={{ color: '#74b9ff' }}>downstream_exon</span>
          </div>
        </>
      )}

      {/* RNA Map Tab */}
      {activeTab === 'rnamap' && rmapsData && (
        <>
          {/* RBP Selector */}
          <div style={{ marginBottom: '1rem' }}>
            <label style={{ fontSize: '0.85rem', color: 'var(--text-muted)', marginRight: '0.5rem' }}>Select RBP:</label>
            <select
              value={selectedRBP || ''}
              onChange={(e) => setSelectedRBP(e.target.value || null)}
              style={{
                background: 'rgba(0,0,0,0.3)',
                border: '1px solid var(--border-color)',
                color: 'white',
                padding: '0.5rem',
                borderRadius: '6px',
                fontSize: '0.85rem',
                minWidth: '200px'
              }}
            >
              <option value="">Select an RBP...</option>
              {rmapsMotifs.map((m, i) => (
                <option key={i} value={m.rbp}>{m.rbp} ({m.region})</option>
              ))}
            </select>
          </div>

          {/* RNA Map Visualization */}
          {selectedRBP && rnaMapChartData && (
            <div style={{ background: 'rgba(0,0,0,0.3)', borderRadius: '8px', padding: '1rem', marginBottom: '1rem' }}>
              <h4 style={{ margin: '0 0 1rem 0', fontSize: '0.9rem', color: '#00cec9' }}>
                RNA Map: {selectedRBP} Motif Density Distribution
              </h4>
              <div style={{ height: 250 }}>
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={rnaMapChartData} layout="vertical" margin={{ top: 5, right: 30, left: 100, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                    <XAxis type="number" stroke="#888" tickFormatter={(v) => `${v}%`} />
                    <YAxis type="category" dataKey="shortName" stroke="#888" width={95} tick={{ fontSize: 11 }} />
                    <Tooltip 
                      contentStyle={{ background: 'rgba(30, 30, 50, 0.95)', border: '1px solid rgba(108, 92, 231, 0.5)', borderRadius: '8px', fontSize: '0.8rem' }}
                      formatter={(value, name) => [`${value.toFixed(2)}%`, name]}
                      labelFormatter={(label, payload) => {
                        if (payload && payload[0]) {
                          const d = payload[0].payload;
                          return (
                            <div style={{ minWidth: 180 }}>
                              <div style={{ fontWeight: 'bold', marginBottom: 4, color: '#00cec9' }}>{d.name}</div>
                              <div>Upregulated: {d.upregulated.toFixed(2)}% (n={d.nUp})</div>
                              <div>Downregulated: {d.downregulated.toFixed(2)}% (n={d.nDown})</div>
                              <div>Background: {d.background.toFixed(2)}% (n={d.nBg})</div>
                            </div>
                          );
                        }
                        return label;
                      }}
                    />
                    <Legend />
                    <Bar dataKey="upregulated" fill="#ff6b6b" name="Upregulated" />
                    <Bar dataKey="downregulated" fill="#74b9ff" name="Downregulated" />
                    <Bar dataKey="background" fill="#888" name="Background" />
                  </BarChart>
                </ResponsiveContainer>
              </div>
              
              {/* Region Explanation */}
              <div style={{ marginTop: '1rem', fontSize: '0.8rem', color: 'var(--text-muted)', lineHeight: 1.6 }}>
                <p style={{ margin: '0 0 0.5rem 0' }}>
                  <strong style={{ color: '#fff' }}>RNA Map Interpretation:</strong>
                </p>
                <ul style={{ margin: 0, paddingLeft: '1.5rem' }}>
                  <li><span style={{ color: '#ff6b6b' }}>Red bars</span>: Motif density in {selectedRBP} upregulated events</li>
                  <li><span style={{ color: '#74b9ff' }}>Blue bars</span>: Motif density in {selectedRBP} downregulated events</li>
                  <li><span style={{ color: '#888' }}>Gray bars</span>: Motif density in background (non-regulated) events</li>
                </ul>
                <p style={{ margin: '0.5rem 0 0 0' }}>
                  Higher bars indicate higher motif density, suggesting {selectedRBP} may bind more frequently in that region.
                </p>
              </div>
            </div>
          )}

          {!selectedRBP && (
            <div style={{ textAlign: 'center', padding: '2rem', color: 'var(--text-muted)', background: 'rgba(0,0,0,0.2)', borderRadius: '8px' }}>
              <TrendingUp size={32} style={{ opacity: 0.3, marginBottom: '0.5rem' }} />
              <p style={{ margin: 0 }}>Select an RBP from the dropdown to view its RNA map</p>
            </div>
          )}
        </>
      )}

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
              This analysis uses <strong>rMAPS-style enrichment</strong> to identify RBPs whose motifs are 
              statistically enriched in splicing events, comparing upregulated/downregulated events against 
              background using Wilcoxon rank sum test.
            </p>
          </div>
        </div>
      )}

      {/* Old motif analysis stats (shown when no rMAPS data) */}
      {!rmapsData && motifData && (
        <>
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
            <div style={{ textAlign: 'center', padding: '2rem', color: 'var(--text-muted)', background: 'rgba(0,0,0,0.2)', borderRadius: '8px', marginBottom: '1.5rem' }}>
              <Zap size={32} style={{ opacity: 0.3, marginBottom: '0.5rem' }} />
              <p style={{ margin: 0 }}>No{showSignificantOnly ? ' significant' : ''} motifs found</p>
            </div>
          )}
        </>
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
