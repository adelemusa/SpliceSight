import React from 'react';
import { AreaChart, Area, XAxis, YAxis, Tooltip, ResponsiveContainer } from 'recharts';

// Questo è un mockup del sashimi plot usando un grafico ad area. 
// Un VERO sashimi plot richiederebbe curve D3 custom, ma per una MVP di Recharts funziona.
const SashimiPlot = ({ event }) => {
  // Generiamo dati mockup a forma di read coverage
  const mockCoverage = Array.from({length: 100}).map((_, i) => {
    let base = 10;
    // Esone 1
    if (i > 10 && i < 30) base = 50 + Math.random() * 20;
    // Cassette Exon (Alternativo)
    if (i > 45 && i < 55) base = event.dpsi > 0 ? 80 : 20; 
    // Esone 2
    if (i > 70 && i < 90) base = 50 + Math.random() * 20;
    
    return {
      position: i,
      condA: base,
      condB: base * (1 - (event.dpsi / 2))
    }
  });

  return (
    <div style={{ height: "300px", width: "100%", marginTop: "1rem" }}>
      <div style={{display: "flex", justifyContent: "space-between", marginBottom: "1rem"}}>
        <span style={{color: "var(--text-muted)"}}>ΔPSI: <strong>{Number(event.dpsi).toFixed(3)}</strong></span>
        <span style={{color: "var(--text-muted)"}}>FDR: <strong>{Number(event.fdr).toExponential(2)}</strong></span>
      </div>
      <ResponsiveContainer width="100%" height="100%">
        <AreaChart data={mockCoverage} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
          <defs>
            <linearGradient id="colorA" x1="0" y1="0" x2="0" y2="1">
              <stop offset="5%" stopColor="#6c5ce7" stopOpacity={0.8}/>
              <stop offset="95%" stopColor="#6c5ce7" stopOpacity={0}/>
            </linearGradient>
            <linearGradient id="colorB" x1="0" y1="0" x2="0" y2="1">
              <stop offset="5%" stopColor="#00cec9" stopOpacity={0.8}/>
              <stop offset="95%" stopColor="#00cec9" stopOpacity={0}/>
            </linearGradient>
          </defs>
          <XAxis dataKey="position" hide />
          <YAxis hide />
          <Tooltip 
             contentStyle={{ backgroundColor: "rgba(15, 17, 26, 0.9)", border: "1px solid var(--border-color)", borderRadius: "8px" }}
             itemStyle={{ color: "#fff" }}
          />
          <Area type="monotone" dataKey="condA" stroke="#6c5ce7" fillOpacity={1} fill="url(#colorA)" />
          <Area type="monotone" dataKey="condB" stroke="#00cec9" fillOpacity={1} fill="url(#colorB)" />
        </AreaChart>
      </ResponsiveContainer>
    </div>
  )
}

export default SashimiPlot;
