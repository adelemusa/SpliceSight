import React, { useState } from 'react';
import './index.css';
import Dashboard from './components/Dashboard';
import { UploadCloud, Folder, FileText, Info } from 'lucide-react';
import axios from 'axios';

function App() {
  const [analyzing, setAnalyzing] = useState(false);
  const [taskId, setTaskId] = useState(null);
  const [data, setData] = useState(null);
  const [inputPath, setInputPath] = useState("GSCs_vs_T98G");
  const [countType, setCountType] = useState("jcec");

  const API_BASE = import.meta.env.VITE_API_URL || "http://backend:8000";

  const handleStartAnalysis = async () => {
    setAnalyzing(true);
    setData(null);
    try {
      const resp = await axios.post(`${API_BASE}/api/analyze`, {
        path: inputPath,
        fdr: 0.05,
        dpsi: 0.1,
        include_jc: countType === "jc",
        include_jcec: countType === "jcec"
      });
      setTaskId(resp.data.task_id);
      
      const poll = setInterval(async () => {
        const statusResp = await axios.get(`${API_BASE}/api/task/${resp.data.task_id}`);
        if(statusResp.data.status === "COMPLETED") {
          clearInterval(poll);
          setData(statusResp.data.result);
          setAnalyzing(false);
        } else if (statusResp.data.status === "FAILED") {
          clearInterval(poll);
          alert("Analysis error: " + statusResp.data.error);
          setAnalyzing(false);
        }
      }, 2000);
      
    } catch (e) {
      console.error(e);
      alert("Connection error to server");
      setAnalyzing(false);
    }
  };

  return (
    <div className="app-container">
      <header>
        <h1>SpliceSight</h1>
        <h2>rMATS downstream analysis & visualization</h2>
      </header>
      
      {!data && (
        <section className="glass-panel upload-section">
          <div className="flex-center" style={{ gap: "1rem", marginBottom: "1rem" }}>
            <Folder size={48} color="var(--accent)" />
            <FileText size={48} color="var(--accent)" />
          </div>
          <p style={{ color: "var(--text-muted)" }}>Enter the path to the rMATS folder or file</p>
          
          <div style={{ width: "100%", maxWidth: "500px", display: "flex", gap: "10px" }}>
            <input 
              type="text" 
              value={inputPath}
              onChange={(e) => setInputPath(e.target.value)}
              placeholder="e.g. GSCs_vs_T98G or GSCs_vs_T98G/SE.MATS.JC.txt"
              className="path-input"
              style={{
                flex: 1,
                background: "rgba(0,0,0,0.3)",
                border: "1px solid var(--border-color)",
                borderRadius: "8px",
                padding: "0.8rem",
                color: "white",
                outline: "none"
              }}
            />
            <button 
              className="btn-upload" 
              onClick={handleStartAnalysis}
              disabled={analyzing}
            >
              {analyzing ? "Analyzing..." : "Start"}
            </button>
          </div>

          <div style={{ width: "100%", maxWidth: "500px", marginTop: "1.5rem" }}>
            <div style={{ display: "flex", alignItems: "center", gap: "0.5rem", marginBottom: "0.75rem", color: "var(--text-muted)", fontSize: "0.9rem" }}>
              <Info size={16} />
              <strong>Read counting method</strong>
            </div>
            <div style={{ display: "flex", gap: "1rem", flexWrap: "wrap" }}>
              <label style={{ display: "flex", alignItems: "center", gap: "0.5rem", cursor: "pointer", padding: "0.5rem 1rem", background: countType === "jc" ? "rgba(108, 92, 231, 0.3)" : "rgba(0,0,0,0.2)", borderRadius: "8px", border: countType === "jc" ? "1px solid rgba(108, 92, 231, 0.5)" : "1px solid transparent" }}>
                <input type="radio" name="countType" value="jc" checked={countType === "jc"} onChange={(e) => setCountType(e.target.value)} />
                <span>JC (Junction Count)</span>
              </label>
              <label style={{ display: "flex", alignItems: "center", gap: "0.5rem", cursor: "pointer", padding: "0.5rem 1rem", background: countType === "jcec" ? "rgba(108, 92, 231, 0.3)" : "rgba(0,0,0,0.2)", borderRadius: "8px", border: countType === "jcec" ? "1px solid rgba(108, 92, 231, 0.5)" : "1px solid transparent" }}>
                <input type="radio" name="countType" value="jcec" checked={countType === "jcec"} onChange={(e) => setCountType(e.target.value)} />
                <span>JCEC (JC + Exon Count)</span>
              </label>
            </div>
            <div style={{ marginTop: "0.75rem", fontSize: "0.8rem", color: "var(--text-muted)", lineHeight: 1.5 }}>
              <strong>JC:</strong> Counts only reads spanning splicing junctions. More conservative, recommended for clean datasets.<br />
              <strong>JCEC:</strong> Also includes reads fully mapped to exons. More sensitive, recommended for complex samples.
            </div>
          </div>
          
          {analyzing && (
            <p style={{ marginTop: "1rem", color: "var(--accent)" }}>
              Analysis in progress... This may take a few minutes for genome indexing.
            </p>
          )}
        </section>
      )}

      {data && (
        <Dashboard data={data} reset={() => setData(null)} />
      )}
    </div>
  );
}

export default App;
