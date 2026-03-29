import React, { useState, useMemo } from 'react';
import { ArrowLeft, Activity, Target, Search, BarChart2, Filter, ArrowUpDown, ArrowUp, ArrowDown, Info, Dna, ChevronLeft, ChevronRight, ChevronFirst, ChevronLast, X } from 'lucide-react';
import SashimiPlot from './SashimiPlot';
import GSEAPanel from './GSEAPanel';
import MotifPanel from './MotifPanel';

const EVENT_TYPES = ['SE', 'A3SS', 'A5SS', 'MXE', 'RI'];
const PAGE_SIZE_OPTIONS = [25, 50, 100, 250, 500];
const DEFAULT_PAGE_SIZE = 50;

const ProteomePanel = ({ events }) => {
  const [showDetails, setShowDetails] = useState(false);

  const stats = useMemo(() => {
    const frameshift = events.filter(e => e.proteome_impact?.includes("Frameshift"));
    const inFrame = events.filter(e => e.proteome_impact?.includes("In-frame"));
    const nmd = events.filter(e => e.nmd_candidate);
    const byType = {};
    
    EVENT_TYPES.forEach(type => {
      const typeEvents = events.filter(e => e.event_type === type);
      byType[type] = {
        total: typeEvents.length,
        frameshift: typeEvents.filter(e => e.proteome_impact?.includes("Frameshift")).length,
        nmd: typeEvents.filter(e => e.nmd_candidate).length
      };
    });
    
    return { frameshift, inFrame, nmd, byType, total: events.length };
  }, [events]);

  return (
    <div className="widget glass-panel">
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <div>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Target size={18} color="var(--accent)" />
            Proteome Impact
          </h3>
          <p style={{ color: 'var(--text-muted)', fontSize: '0.8rem', marginTop: '0.5rem' }}>
            In silico analysis of splicing events impact on protein translation
          </p>
        </div>
        <button
          onClick={() => setShowDetails(!showDetails)}
          style={{
            background: showDetails ? "var(--accent)" : "rgba(108, 92, 231, 0.2)",
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
          <Info size={14} /> {showDetails ? "Hide" : "Details"}
        </button>
      </div>

      {showDetails && (
        <div style={{
          background: 'rgba(0, 206, 201, 0.1)',
          border: '1px solid rgba(0, 206, 201, 0.3)',
          borderRadius: '8px',
          padding: '1rem',
          marginTop: '1rem'
        }}>
          <h4 style={{ margin: '0 0 0.5rem 0', color: '#00cec9', fontSize: '0.9rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Dna size={16} /> How does proteome impact analysis work?
          </h4>
          <div style={{ fontSize: '0.85rem', color: 'var(--text-muted)', lineHeight: 1.6 }}>
            <p style={{ margin: '0 0 0.75rem 0' }}>
              <strong>Frameshift:</strong> When the exon length is not a multiple of 3, translation is altered. 
              The resulting protein has a completely different amino acid sequence downstream of the splicing site.
            </p>
            <p style={{ margin: '0 0 0.75rem 0' }}>
              <strong>In-frame:</strong> Insertions/deletions are multiples of 3, so the sequence remains in-frame. 
              Impact depends on position within functional domains.
            </p>
            <p style={{ margin: '0 0 0 0' }}>
              <strong>NMD:</strong> If a premature stop codon (PTC) is introduced more than 50-55 nt upstream 
              of the last exon-exon junction, the transcript is degraded before translation.
            </p>
          </div>
        </div>
      )}

      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '1rem', marginTop: '1.5rem' }}>
        <div style={{ background: 'rgba(255, 118, 117, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(255, 118, 117, 0.3)' }}>
          <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#ff7675' }}>{stats.frameshift.length}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Frameshift</div>
          <div style={{ fontSize: '0.75rem', color: '#ff7675' }}>({(stats.frameshift.length / stats.total * 100).toFixed(1)}%)</div>
        </div>
        <div style={{ background: 'rgba(85, 239, 196, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(85, 239, 196, 0.3)' }}>
          <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#55efc4' }}>{stats.inFrame.length}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>In-frame</div>
          <div style={{ fontSize: '0.75rem', color: '#55efc4' }}>({(stats.inFrame.length / stats.total * 100).toFixed(1)}%)</div>
        </div>
        <div style={{ background: 'rgba(232, 67, 147, 0.15)', borderRadius: '8px', padding: '1rem', textAlign: 'center', border: '1px solid rgba(232, 67, 147, 0.3)' }}>
          <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#e84393' }}>{stats.nmd.length}</div>
          <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>NMD candidates</div>
          <div style={{ fontSize: '0.75rem', color: '#e84393' }}>({(stats.nmd.length / stats.total * 100).toFixed(1)}%)</div>
        </div>
      </div>

      <div style={{ marginTop: '1.5rem' }}>
        <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: 'var(--text-muted)' }}>Distribution by event type</h4>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(5, 1fr)', gap: '0.5rem' }}>
          {EVENT_TYPES.map(type => {
            const data = stats.byType[type];
            const fsPercent = data.total > 0 ? (data.frameshift / data.total * 100) : 0;
            return (
              <div key={type} style={{ background: 'rgba(0,0,0,0.2)', borderRadius: '6px', padding: '0.75rem', textAlign: 'center' }}>
                <div style={{ fontWeight: 'bold', fontSize: '1.1rem' }}>{type}</div>
                <div style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>{data.total} events</div>
                <div style={{ marginTop: '0.5rem', display: 'flex', gap: '2px', height: '8px' }}>
                  <div style={{ flex: fsPercent, background: '#ff7675', borderRadius: '4px 0 0 4px' }} />
                  <div style={{ flex: 100 - fsPercent, background: '#55efc4', borderRadius: '0 4px 4px 0' }} />
                </div>
                <div style={{ fontSize: '0.7rem', marginTop: '0.25rem', color: '#ff7675' }}>{data.frameshift} FS | {data.nmd} NMD</div>
              </div>
            );
          })}
        </div>
      </div>

      <div style={{ marginTop: '1.5rem' }}>
        <h4 style={{ margin: '0 0 0.75rem 0', fontSize: '0.9rem', color: '#ff7675' }}>Frameshift events ({stats.frameshift.length})</h4>
        <div style={{ maxHeight: '180px', overflowY: 'auto', borderRadius: '6px', border: '1px solid var(--border-color)' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.8rem' }}>
            <thead style={{ background: 'rgba(255, 118, 117, 0.2)', position: 'sticky', top: 0 }}>
              <tr>
                <th style={{ padding: '6px 8px', textAlign: 'left', fontWeight: 'normal', color: 'var(--text-muted)' }}>Gene</th>
                <th style={{ padding: '6px 8px', textAlign: 'left', fontWeight: 'normal', color: 'var(--text-muted)' }}>Type</th>
                <th style={{ padding: '6px 8px', textAlign: 'right', fontWeight: 'normal', color: 'var(--text-muted)' }}>Impact</th>
              </tr>
            </thead>
            <tbody>
              {stats.frameshift.slice(0, 30).map((evt, idx) => (
                <tr key={idx} style={{ borderBottom: '1px solid rgba(255,255,255,0.05)' }}>
                  <td style={{ padding: '6px 8px', fontWeight: 'bold' }}>{evt.gene_symbol}</td>
                  <td style={{ padding: '6px 8px' }}>
                    <span style={{ background: 'rgba(255,255,255,0.1)', padding: '2px 4px', borderRadius: '4px', fontSize: '0.7rem' }}>{evt.event_type}</span>
                  </td>
                  <td style={{ padding: '6px 8px', textAlign: 'right', color: '#ff7675', fontSize: '0.75rem' }}>{evt.proteome_impact}</td>
                </tr>
              ))}
            </tbody>
          </table>
          {stats.frameshift.length > 30 && (
            <div style={{ padding: '0.5rem', textAlign: 'center', color: 'var(--text-muted)', fontSize: '0.75rem', background: 'rgba(0,0,0,0.2)' }}>
              ... and {stats.frameshift.length - 30} more events
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

const Dashboard = ({ data, reset }) => {
  const [selectedEvent, setSelectedEvent] = useState(null);
  const [searchTerm, setSearchTerm] = useState("");
  const [sortColumns, setSortColumns] = useState([{ field: "fdr", direction: "asc" }]);
  const [filterType, setFilterType] = useState("all");
  const [filterNMD, setFilterNMD] = useState(false);
  const [filterFDR, setFilterFDR] = useState(0.05);
  const [showFilters, setShowFilters] = useState(false);
  const [currentPage, setCurrentPage] = useState(1);
  const [pageSize, setPageSize] = useState(DEFAULT_PAGE_SIZE);
  const [showSortModal, setShowSortModal] = useState(false);
  
  const payload = data?.data || {};
  const events = payload.events || [];
  const enrichments = payload.enrichment || [];
  const motifData = payload.motif_analysis || null;
  
  const filteredAndSortedEvents = useMemo(() => {
    let result = [...events];
    
    if (searchTerm) {
      const term = searchTerm.toLowerCase();
      result = result.filter(e => 
        e.gene_symbol?.toLowerCase().includes(term) || 
        e.event_type?.toLowerCase().includes(term)
      );
    }
    
    if (filterType !== "all") result = result.filter(e => e.event_type === filterType);
    if (filterNMD) result = result.filter(e => e.nmd_candidate);
    result = result.filter(e => e.fdr < filterFDR);
    
    result.sort((a, b) => {
      for (const sort of sortColumns) {
        let aVal = a[sort.field];
        let bVal = b[sort.field];
        if (sort.field === "abs_dpsi") {
          aVal = Math.abs(a.dpsi);
          bVal = Math.abs(b.dpsi);
        }
        if (aVal == null) aVal = sort.direction === "asc" ? Infinity : -Infinity;
        if (bVal == null) bVal = sort.direction === "asc" ? Infinity : -Infinity;
        if (typeof aVal === "string") { aVal = aVal.toLowerCase(); bVal = bVal.toLowerCase(); }
        if (aVal < bVal) return sort.direction === "asc" ? -1 : 1;
        if (aVal > bVal) return sort.direction === "asc" ? 1 : -1;
      }
      return 0;
    });
    
    return result;
  }, [events, searchTerm, filterType, filterNMD, filterFDR, sortColumns]);

  const totalPages = Math.ceil(filteredAndSortedEvents.length / pageSize);
  const paginatedEvents = filteredAndSortedEvents.slice((currentPage - 1) * pageSize, currentPage * pageSize);

  useMemo(() => {
    if (currentPage > totalPages) setCurrentPage(1);
  }, [totalPages, currentPage]);

  const handleSort = (field) => {
    const existing = sortColumns.find(s => s.field === field);
    if (existing) {
      if (existing.direction === "asc") {
        if (sortColumns.length === 1) {
          setSortColumns([{ field, direction: "desc" }]);
        } else {
          setSortColumns(sortColumns.filter(s => s.field !== field).map((s, i) => i === 0 ? { ...s, direction: "desc" } : s));
        }
      } else {
        setSortColumns(sortColumns.filter(s => s.field !== field));
      }
    } else {
      setSortColumns([...sortColumns, { field, direction: "asc" }]);
    }
    setCurrentPage(1);
  };

  const SortIcon = ({ field }) => {
    const idx = sortColumns.findIndex(s => s.field === field);
    if (idx === -1) return <ArrowUpDown size={12} style={{ opacity: 0.3 }} />;
    return (
      <span style={{ display: 'flex', alignItems: 'center', gap: '2px' }}>
        {sortColumns.length > 1 && <span style={{ fontSize: '10px', opacity: 0.7 }}>{idx + 1}</span>}
        {sortColumns[idx].direction === "asc" ? <ArrowUp size={12} /> : <ArrowDown size={12} />}
      </span>
    );
  };

  const addSortColumn = (field) => {
    if (!sortColumns.find(s => s.field === field)) {
      setSortColumns([...sortColumns, { field, direction: "asc" }]);
    }
    setShowSortModal(false);
  };

  const removeSortColumn = (field, e) => {
    e.stopPropagation();
    setSortColumns(sortColumns.filter(s => s.field !== field));
  };

  const getEventTypeColor = (type) => ({ SE: '#6c5ce7', A3SS: '#00cec9', A5SS: '#fdcb6e', MXE: '#e84393', RI: '#74b9ff' }[type] || '#888');

  const sortableColumns = [
    { field: 'gene_symbol', label: 'Gene' },
    { field: 'event_type', label: 'Type' },
    { field: 'fdr', label: 'FDR' },
    { field: 'dpsi', label: 'dPSI' },
    { field: 'abs_dpsi', label: '|dPSI|' },
    { field: 'p_value', label: 'P-value' },
  ];

  const stats = useMemo(() => ({
    total: events.length,
    byType: EVENT_TYPES.reduce((acc, type) => { acc[type] = events.filter(e => e.event_type === type).length; return acc; }, {}),
    nmdCount: events.filter(e => e.nmd_candidate).length
  }), [events]);

  return (
    <div className="dashboard-view animation-fade-in" style={{ maxWidth: '1400px', margin: '0 auto', padding: '1rem' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.5rem' }}>
        <button className="btn-upload" style={{ padding: '0.5rem 1rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }} onClick={reset}>
          <ArrowLeft size={16} /> New Analysis
        </button>
        <h2 style={{ margin: 0, color: 'var(--accent)', fontSize: '1.2rem' }}>
          {events.length} significant events | {stats.nmdCount} NMD candidates
        </h2>
      </div>

      <ProteomePanel events={events} />

      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '1rem', marginTop: '1.5rem' }}>
        <div className="widget glass-panel" style={{ textAlign: 'center' }}>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', justifyContent: 'center', gap: '0.5rem' }}>
            <Activity size={18} color="var(--accent)" /> NMD Prediction
          </h3>
          <div style={{ marginTop: '1rem' }}>
            <div style={{ fontSize: '2.5rem', fontWeight: 'bold', color: '#e84393' }}>{stats.nmdCount}</div>
            <div style={{ color: 'var(--text-muted)', fontSize: '0.85rem' }}>{(stats.nmdCount / stats.total * 100).toFixed(1)}% of total</div>
          </div>
        </div>

        <div className="widget glass-panel">
          <h3 style={{ margin: '0 0 1rem 0', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <span style={{ color: 'var(--accent)' }}>Event Distribution</span>
          </h3>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
            {EVENT_TYPES.map(type => (
              <div key={type} style={{ display: 'flex', alignItems: 'center', gap: '0.75rem' }}>
                <span style={{ width: '45px', fontSize: '0.85rem', fontWeight: 'bold', color: getEventTypeColor(type) }}>{type}</span>
                <div style={{ flex: 1, background: 'rgba(255,255,255,0.1)', borderRadius: '4px', height: '18px', overflow: 'hidden' }}>
                  <div style={{ width: `${(stats.byType[type] / stats.total * 100)}%`, height: '100%', background: getEventTypeColor(type) }} />
                </div>
                <span style={{ width: '50px', textAlign: 'right', fontSize: '0.85rem' }}>{stats.byType[type]}</span>
              </div>
            ))}
          </div>
        </div>
      </div>

      <GSEAPanel enrichments={enrichments} />

      <MotifPanel motifData={motifData} />

      <div className="widget glass-panel" style={{ marginTop: '1.5rem' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem', flexWrap: 'wrap', gap: '10px' }}>
          <h3 style={{ margin: 0, display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <Search size={18} color="var(--accent)" />
            Explore Events
            <span style={{ fontSize: '0.85rem', color: 'var(--text-muted)', fontWeight: 'normal' }}>
              ({filteredAndSortedEvents.length} events)
            </span>
          </h3>
          <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
            <input 
              type="text" 
              placeholder="Search gene..." 
              value={searchTerm}
              onChange={(e) => { setSearchTerm(e.target.value); setCurrentPage(1); }}
              style={{ background: 'rgba(0,0,0,0.2)', border: '1px solid var(--border-color)', color: 'white', padding: '0.5rem 0.75rem', borderRadius: '6px', fontSize: '0.85rem', width: '180px', outline: 'none' }}
            />
            <select value={filterType} onChange={(e) => { setFilterType(e.target.value); setCurrentPage(1); }}
              style={{ background: 'rgba(0,0,0,0.2)', border: '1px solid var(--border-color)', color: 'white', padding: '0.5rem', borderRadius: '6px', fontSize: '0.85rem', outline: 'none' }}>
              <option value="all">All types</option>
              {EVENT_TYPES.map(type => <option key={type} value={type}>{type} ({stats.byType[type]})</option>)}
            </select>
            <label style={{ display: 'flex', alignItems: 'center', gap: '0.3rem', fontSize: '0.85rem', color: 'var(--text-muted)', cursor: 'pointer' }}>
              <input type="checkbox" checked={filterNMD} onChange={(e) => { setFilterNMD(e.target.checked); setCurrentPage(1); }} />
              NMD
            </label>
            <button onClick={() => setShowFilters(!showFilters)}
              style={{ background: showFilters ? "var(--accent)" : "rgba(108, 92, 231, 0.2)", border: "none", color: "white", padding: "0.5rem 1rem", borderRadius: "6px", cursor: "pointer", display: "flex", alignItems: "center", gap: "0.5rem", fontSize: "0.85rem" }}>
              <Filter size={14} /> Filters
            </button>
          </div>
        </div>

        {showFilters && (
          <div style={{ background: "rgba(0,0,0,0.2)", padding: "1rem", borderRadius: "8px", marginBottom: "1rem", display: "flex", alignItems: "center", gap: "2rem", flexWrap: "wrap" }}>
            <div>
              <label style={{ fontSize: "0.85rem", color: "var(--text-muted)", display: "block", marginBottom: "0.25rem" }}>Max FDR: {filterFDR.toFixed(3)}</label>
              <input type="range" min="0.001" max="0.05" step="0.001" value={filterFDR} onChange={(e) => { setFilterFDR(parseFloat(e.target.value)); setCurrentPage(1); }} style={{ width: '150px' }} />
            </div>
            <button onClick={() => { setFilterType("all"); setFilterNMD(false); setFilterFDR(0.05); setSearchTerm(""); setCurrentPage(1); }}
              style={{ background: "rgba(255,255,255,0.1)", border: "none", color: "white", padding: "0.5rem 1rem", borderRadius: "6px", cursor: "pointer", fontSize: "0.85rem" }}>
              Reset
            </button>
          </div>
        )}

        <div style={{ overflowX: 'auto', borderRadius: '8px', border: '1px solid var(--border-color)' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse', textAlign: 'left', fontSize: '0.85rem' }}>
            <thead style={{ background: 'rgba(108, 92, 231, 0.2)', position: 'sticky', top: 0 }}>
              <tr>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)', cursor: 'pointer', whiteSpace: 'nowrap' }} onClick={() => handleSort("gene_symbol")}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}>Gene <SortIcon field="gene_symbol" /></span>
                </th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)', cursor: 'pointer', whiteSpace: 'nowrap' }} onClick={() => handleSort("event_type")}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}>Type <SortIcon field="event_type" /></span>
                </th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)', cursor: 'pointer', whiteSpace: 'nowrap' }} onClick={() => handleSort("fdr")}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}>FDR <SortIcon field="fdr" /></span>
                </th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)', cursor: 'pointer', whiteSpace: 'nowrap' }} onClick={() => handleSort("dpsi")}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}>dPSI <SortIcon field="dpsi" /></span>
                </th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)', cursor: 'pointer', whiteSpace: 'nowrap' }} onClick={() => handleSort("abs_dpsi")}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}>|dPSI| <SortIcon field="abs_dpsi" /></span>
                </th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)' }}>Impact</th>
                <th style={{ padding: '10px 12px', borderBottom: '1px solid var(--border-color)' }}>Actions</th>
              </tr>
            </thead>
            <tbody>
              {paginatedEvents.map((evt, idx) => (
                <tr key={idx} style={{ background: selectedEvent === evt ? "rgba(108, 92, 231, 0.15)" : "transparent", borderBottom: '1px solid rgba(255,255,255,0.05)' }}>
                  <td style={{ padding: '8px 12px', fontWeight: 'bold' }}>{evt.gene_symbol}</td>
                  <td style={{ padding: '8px 12px' }}>
                    <span style={{ background: getEventTypeColor(evt.event_type), padding: '2px 6px', borderRadius: '4px', fontSize: '0.75rem', color: 'white', fontWeight: 'bold' }}>{evt.event_type}</span>
                  </td>
                  <td style={{ padding: '8px 12px', color: evt.fdr < 0.01 ? "#00cec9" : "inherit", fontFamily: 'monospace', fontSize: '0.8rem' }}>{Number(evt.fdr).toExponential(2)}</td>
                  <td style={{ padding: '8px 12px', color: evt.dpsi > 0 ? "#55efc4" : "#ff7675", fontFamily: 'monospace', fontSize: '0.8rem' }}>{Number(evt.dpsi) >= 0 ? '+' : ''}{Number(evt.dpsi).toFixed(3)}</td>
                  <td style={{ padding: '8px 12px', fontFamily: 'monospace', fontSize: '0.8rem', fontWeight: 'bold' }}>{Math.abs(Number(evt.dpsi)).toFixed(3)}</td>
                  <td style={{ padding: '8px 12px', fontSize: '0.8rem' }}>
                    <span style={{ color: evt.proteome_impact?.includes("Frameshift") ? "#ff7675" : "#55efc4" }}>{evt.proteome_impact?.split(' ')[0]}</span>
                    {evt.nmd_candidate && <span style={{ color: "#e84393", fontWeight: "bold", marginLeft: "4px", fontSize: "0.7rem" }}>NMD</span>}
                  </td>
                  <td style={{ padding: '8px 12px' }}>
                    <button onClick={() => setSelectedEvent(evt)} style={{ background: "var(--accent)", color: "white", border: "none", padding: "4px 10px", borderRadius: "4px", cursor: "pointer", fontSize: "0.8rem", display: 'flex', alignItems: 'center', gap: '4px' }}>
                      <BarChart2 size={12} /> Plot
                    </button>
                  </td>
                </tr>
              ))}
              {paginatedEvents.length === 0 && (
                <tr><td colSpan="6" style={{ padding: '30px', textAlign: 'center', color: 'var(--text-muted)' }}>No events found</td></tr>
              )}
            </tbody>
          </table>
        </div>

        {totalPages > 0 && (
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginTop: '1rem', paddingTop: '1rem', borderTop: '1px solid var(--border-color)', flexWrap: 'wrap', gap: '1rem' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '1rem', flexWrap: 'wrap' }}>
              <span style={{ fontSize: '0.85rem', color: 'var(--text-muted)' }}>
                Show
                <select 
                  value={pageSize} 
                  onChange={(e) => { setPageSize(Number(e.target.value)); setCurrentPage(1); }}
                  style={{ background: 'rgba(0,0,0,0.3)', border: '1px solid var(--border-color)', color: 'white', padding: '0.25rem 0.5rem', borderRadius: '4px', margin: '0 0.5rem', fontSize: '0.85rem' }}
                >
                  {PAGE_SIZE_OPTIONS.map(size => <option key={size} value={size}>{size}</option>)}
                </select>
                of {filteredAndSortedEvents.length} events
              </span>
              <span style={{ fontSize: '0.85rem', color: 'var(--text-muted)' }}>
                Page {currentPage} of {totalPages}
              </span>
            </div>
            
            <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
              {sortColumns.length > 0 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', flexWrap: 'wrap' }}>
                  <span style={{ fontSize: '0.8rem', color: 'var(--text-muted)' }}>Sort:</span>
                  {sortColumns.map((s, idx) => (
                    <span key={s.field} style={{ background: 'rgba(108, 92, 231, 0.3)', padding: '2px 8px', borderRadius: '4px', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '4px' }}>
                      {idx + 1}. {sortableColumns.find(c => c.field === s.field)?.label} {s.direction === 'asc' ? '↑' : '↓'}
                      <X size={12} style={{ cursor: 'pointer', opacity: 0.7 }} onClick={(e) => removeSortColumn(s.field, e)} />
                    </span>
                  ))}
                </div>
              )}
              <button 
                onClick={() => setShowSortModal(!showSortModal)}
                style={{ background: 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.75rem', borderRadius: '4px', cursor: 'pointer', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '0.3rem' }}
              >
                <ArrowUpDown size={12} /> Sort
              </button>
              {showSortModal && (
                <div style={{ position: 'absolute', background: 'rgba(30, 30, 50, 0.98)', border: '1px solid var(--border-color)', borderRadius: '8px', padding: '0.5rem', zIndex: 100, boxShadow: '0 4px 12px rgba(0,0,0,0.5)' }}>
                  {sortableColumns.map(col => {
                    const isActive = sortColumns.find(s => s.field === col.field);
                    return (
                      <div 
                        key={col.field}
                        onClick={() => addSortColumn(col.field)}
                        style={{ padding: '0.5rem 1rem', cursor: 'pointer', borderRadius: '4px', background: isActive ? 'rgba(108, 92, 231, 0.3)' : 'transparent', color: 'white', fontSize: '0.85rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}
                      >
                        <span>{col.label}</span>
                        {isActive && <span style={{ fontSize: '0.75rem', color: 'var(--accent)' }}>Click to remove</span>}
                      </div>
                    );
                  })}
                </div>
              )}
            </div>
            
            <div style={{ display: 'flex', gap: '0.25rem' }}>
              <button onClick={() => setCurrentPage(1)} disabled={currentPage === 1} style={{ background: 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.6rem', borderRadius: '4px', cursor: currentPage === 1 ? 'not-allowed' : 'pointer', opacity: currentPage === 1 ? 0.5 : 1 }}>
                <ChevronFirst size={16} />
              </button>
              <button onClick={() => setCurrentPage(p => Math.max(1, p - 1))} disabled={currentPage === 1} style={{ background: 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.6rem', borderRadius: '4px', cursor: currentPage === 1 ? 'not-allowed' : 'pointer', opacity: currentPage === 1 ? 0.5 : 1 }}>
                <ChevronLeft size={16} />
              </button>
              {[...Array(Math.min(5, totalPages))].map((_, i) => {
                let page = i + 1;
                if (totalPages > 5) {
                  if (currentPage <= 3) page = i + 1;
                  else if (currentPage >= totalPages - 2) page = totalPages - 4 + i;
                  else page = currentPage - 2 + i;
                }
                return (
                  <button key={page} onClick={() => setCurrentPage(page)} style={{ background: currentPage === page ? 'var(--accent)' : 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.75rem', borderRadius: '4px', cursor: 'pointer', fontSize: '0.85rem' }}>
                    {page}
                  </button>
                );
              })}
              <button onClick={() => setCurrentPage(p => Math.min(totalPages, p + 1))} disabled={currentPage === totalPages} style={{ background: 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.6rem', borderRadius: '4px', cursor: currentPage === totalPages ? 'not-allowed' : 'pointer', opacity: currentPage === totalPages ? 0.5 : 1 }}>
                <ChevronRight size={16} />
              </button>
              <button onClick={() => setCurrentPage(totalPages)} disabled={currentPage === totalPages} style={{ background: 'rgba(108, 92, 231, 0.2)', border: 'none', color: 'white', padding: '0.4rem 0.6rem', borderRadius: '4px', cursor: currentPage === totalPages ? 'not-allowed' : 'pointer', opacity: currentPage === totalPages ? 0.5 : 1 }}>
                <ChevronLast size={16} />
              </button>
            </div>
          </div>
        )}
      </div>
        
      {selectedEvent && (
        <div className="widget glass-panel" style={{ marginTop: '1.5rem' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
            <h3 style={{ margin: 0 }}>
              Sashimi Plot: <span style={{ color: "#00cec9" }}>{selectedEvent.gene_symbol}</span> ({selectedEvent.event_type})
            </h3>
            <button onClick={() => setSelectedEvent(null)} style={{ background: 'rgba(255,255,255,0.1)', border: 'none', color: 'white', padding: '6px 12px', borderRadius: '4px', cursor: 'pointer' }}>Close</button>
          </div>
          <SashimiPlot event={selectedEvent} />
        </div>
      )}
    </div>
  );
};

export default Dashboard;
