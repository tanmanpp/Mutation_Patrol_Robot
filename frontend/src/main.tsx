import React, { useEffect, useMemo, useState } from "react";
import { createRoot } from "react-dom/client";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  ComposedChart,
  LabelList,
  Legend,
  Line,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import "./styles.css";

const API_BASE = window.location.port === "5173" ? "http://127.0.0.1:8000" : "";

type Database = {
  name: string;
  gene_count: number;
  warning_count: number;
  genes: Gene[];
};

type Gene = {
  symbol: string;
  description: string;
  chrom?: string;
  gene_start?: string | number;
  gene_end?: string | number;
  sequence_status?: string;
};

type ResultPayload = {
  run_id: string;
  sample_summary: Record<string, string>[];
  mutation_candidates: Record<string, string>[];
  site_query: Record<string, string>[];
};

function Field(props: { label: string; children: React.ReactNode }) {
  return (
    <label className="field">
      <span>{props.label}</span>
      {props.children}
    </label>
  );
}

function FolderInput(props: { name: string }) {
  return React.createElement("input", {
    name: props.name,
    type: "file",
    multiple: true,
    webkitdirectory: "",
    directory: "",
  } as React.InputHTMLAttributes<HTMLInputElement> & { webkitdirectory: string; directory: string });
}

function DataTable({ rows, empty }: { rows: Record<string, string>[]; empty: string }) {
  const columns = useMemo(() => (rows[0] ? Object.keys(rows[0]) : []), [rows]);
  if (!rows.length) return <div className="empty">{empty}</div>;
  return (
    <div className="tableWrap">
      <table>
        <thead>
          <tr>{columns.map((column) => <th key={column}>{column}</th>)}</tr>
        </thead>
        <tbody>
          {rows.map((row, rowIndex) => (
            <tr key={rowIndex}>
              {columns.map((column) => <td key={column}>{row[column]}</td>)}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

const PAGE_SIZE = 20;

function compareCellValues(left: string | undefined, right: string | undefined) {
  const leftText = left || "";
  const rightText = right || "";
  const leftNumber = Number(leftText);
  const rightNumber = Number(rightText);
  if (leftText !== "" && rightText !== "" && Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
    return leftNumber - rightNumber;
  }
  return leftText.localeCompare(rightText, undefined, { numeric: true, sensitivity: "base" });
}

function PaginatedSortableTable({ rows, empty }: { rows: Record<string, string>[]; empty: string }) {
  const columns = useMemo(() => (rows[0] ? Object.keys(rows[0]) : []), [rows]);
  const [sortColumn, setSortColumn] = useState("");
  const [sortDirection, setSortDirection] = useState<"asc" | "desc">("asc");
  const [page, setPage] = useState(1);

  const sortedRows = useMemo(() => {
    if (!sortColumn) return rows;
    const direction = sortDirection === "asc" ? 1 : -1;
    return [...rows].sort((left, right) => (
      compareCellValues(left[sortColumn], right[sortColumn]) * direction
    ));
  }, [rows, sortColumn, sortDirection]);

  const totalPages = Math.max(1, Math.ceil(sortedRows.length / PAGE_SIZE));
  const currentPage = Math.min(page, totalPages);
  const pageRows = sortedRows.slice((currentPage - 1) * PAGE_SIZE, currentPage * PAGE_SIZE);

  useEffect(() => {
    setPage(1);
  }, [rows, sortColumn, sortDirection]);

  function changeSort(column: string) {
    if (column === sortColumn) {
      setSortDirection((previous) => previous === "asc" ? "desc" : "asc");
    } else {
      setSortColumn(column);
      setSortDirection("asc");
    }
  }

  if (!rows.length) return <div className="empty">{empty}</div>;

  return (
    <div className="pagedTable">
      <div className="tableToolbar">
        <span>
          Showing {(currentPage - 1) * PAGE_SIZE + 1}-{Math.min(currentPage * PAGE_SIZE, sortedRows.length)} of {sortedRows.length}
        </span>
        <div className="pager">
          <button type="button" onClick={() => setPage(1)} disabled={currentPage === 1}>First</button>
          <button type="button" onClick={() => setPage((value) => Math.max(1, value - 1))} disabled={currentPage === 1}>Previous</button>
          <span>Page {currentPage} / {totalPages}</span>
          <button type="button" onClick={() => setPage((value) => Math.min(totalPages, value + 1))} disabled={currentPage === totalPages}>Next</button>
          <button type="button" onClick={() => setPage(totalPages)} disabled={currentPage === totalPages}>Last</button>
        </div>
      </div>
      <div className="tableWrap">
        <table>
          <thead>
            <tr>
              {columns.map((column) => (
                <th key={column}>
                  <button className="sortHeader" type="button" onClick={() => changeSort(column)}>
                    <span>{column}</span>
                    <span>{sortColumn === column ? (sortDirection === "asc" ? "Asc" : "Desc") : ""}</span>
                  </button>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {pageRows.map((row, rowIndex) => (
              <tr key={rowIndex}>
                {columns.map((column) => <td key={column}>{row[column]}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

type SiteChartDatum = {
  label: string;
  gene: string;
  geneDescription: string;
  position: string;
  allele: string;
  depth: number;
  alleleFreq: number;
  filter: string;
  variantType: string;
  aaChange: string;
};

type SiteCombinedDatum = {
  label: string;
  geneDescription: string;
  depth: number;
  mutationFreq: number;
  filter: string;
};

function parseNumber(value: string | undefined) {
  const parsed = Number(value || 0);
  return Number.isFinite(parsed) ? parsed : 0;
}

function siteQueryLabel(row: Record<string, string>) {
  const gene = row["Gene symbol"] || row.gene || "";
  const aaPos = row.aa_pos || "";
  const aaChange = row.aa_change && row.aa_change !== "NA" ? row.aa_change : "";
  const position = aaPos || row.pos || "";
  const allele = row.alt_codon || row.alt || "";
  if (aaChange) return `${gene} ${aaChange}`;
  if (position && allele) return `${gene} ${position} ${allele}`;
  return `${gene} ${position}`.trim();
}

function siteQueryChartData(rows: Record<string, string>[]): SiteChartDatum[] {
  return rows.map((row, index) => ({
    label: siteQueryLabel(row) || `site ${index + 1}`,
    gene: row["Gene symbol"] || row.gene || "",
    geneDescription: row["Gene name"] || "",
    position: row.aa_pos || row.pos || "",
    allele: row.alt_codon || row.alt || "",
    depth: parseNumber(row.depth),
    alleleFreq: parseNumber(row.allele_freq),
    filter: row.filter || "",
    variantType: row.variant_type || "",
    aaChange: row.aa_change || "",
  }));
}

function sitePositionKey(item: SiteChartDatum) {
  return `${item.gene}|${item.position}`;
}

function sitePositionLabel(item: SiteChartDatum) {
  return item.gene && item.position ? `${item.gene} ${item.position}` : item.label;
}

function siteQueryCombinedData(data: SiteChartDatum[]): SiteCombinedDatum[] {
  const grouped = new Map<string, SiteCombinedDatum>();
  for (const item of data) {
    const key = sitePositionKey(item);
    const current = grouped.get(key) || {
      label: sitePositionLabel(item),
      geneDescription: item.geneDescription,
      depth: item.depth,
      mutationFreq: 0,
      filter: item.filter,
    };
    current.depth = Math.max(current.depth, item.depth);
    current.geneDescription = current.geneDescription || item.geneDescription;
    if (current.filter !== "NO_MAPPING" && item.filter === "NO_MAPPING") {
      current.filter = item.filter;
    }
    if (item.variantType !== "REF") {
      current.mutationFreq += item.alleleFreq;
    }
    grouped.set(key, current);
  }
  return Array.from(grouped.values()).map((item) => ({
    ...item,
    mutationFreq: Math.min(item.mutationFreq, 1),
  }));
}

function chartColor(item: SiteChartDatum) {
  if (item.filter === "NO_MAPPING") return "#a3adb0";
  if (item.variantType === "REF") return "#79a99b";
  if (item.variantType === "CODON") return "#b85c38";
  return "#2f6f9f";
}

function SiteTooltip({ active, payload, label }: { active?: boolean; payload?: any[]; label?: string }) {
  if (!active || !payload?.length) return null;
  const item = payload[0]?.payload || {};
  return (
    <div className="chartTooltip">
      <strong>{label}</strong>
      {item.geneDescription && <span>{item.geneDescription}</span>}
      {item.gene && <span>Gene: {item.gene}</span>}
      {item.position && <span>Position: {item.position}</span>}
      {item.allele && <span>Allele: {item.allele}</span>}
      {item.filter && <span>Status: {item.filter}</span>}
      {payload.map((entry) => (
        <span key={entry.name}>
          {entry.name}: {entry.name?.includes("frequency")
            ? `${(Number(entry.value) * 100).toFixed(2)}%`
            : entry.value}
        </span>
      ))}
    </div>
  );
}

function SiteQueryCharts({ rows }: { rows: Record<string, string>[] }) {
  const data = useMemo(() => siteQueryChartData(rows), [rows]);
  const combinedData = useMemo(() => siteQueryCombinedData(data), [data]);
  if (!data.length) return <div className="empty">No site query chart data loaded.</div>;

  return (
    <div className="chartGrid">
      <section className="chartPanel">
        <h3>Depth by Queried Allele</h3>
        <ResponsiveContainer width="100%" height={260}>
          <BarChart data={data} margin={{ top: 12, right: 18, bottom: 54, left: 0 }}>
            <CartesianGrid strokeDasharray="3 3" vertical={false} />
            <XAxis dataKey="label" angle={-35} textAnchor="end" height={70} interval={0} tick={{ fontSize: 11 }} />
            <YAxis allowDecimals={false} />
            <Tooltip content={<SiteTooltip />} />
            <Bar dataKey="depth" name="Depth">
              {data.map((item, index) => <Cell key={index} fill={chartColor(item)} />)}
              <LabelList dataKey="depth" position="top" fontSize={11} />
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </section>

      <section className="chartPanel">
        <h3>Allele Frequency</h3>
        <ResponsiveContainer width="100%" height={260}>
          <BarChart data={data} margin={{ top: 12, right: 18, bottom: 54, left: 0 }}>
            <CartesianGrid strokeDasharray="3 3" vertical={false} />
            <XAxis dataKey="label" angle={-35} textAnchor="end" height={70} interval={0} tick={{ fontSize: 11 }} />
            <YAxis domain={[0, 1]} tickFormatter={(value) => `${Math.round(Number(value) * 100)}%`} />
            <Tooltip content={<SiteTooltip />} />
            <Bar dataKey="alleleFreq" name="Allele frequency">
              {data.map((item, index) => <Cell key={index} fill={chartColor(item)} />)}
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </section>

      <section className="chartPanel chartPanelWide">
        <h3>Depth and Mutation Frequency</h3>
        <ResponsiveContainer width="100%" height={300}>
          <ComposedChart data={combinedData} margin={{ top: 12, right: 24, bottom: 58, left: 0 }}>
            <CartesianGrid strokeDasharray="3 3" vertical={false} />
            <XAxis dataKey="label" angle={-30} textAnchor="end" height={76} interval={0} tick={{ fontSize: 11 }} />
            <YAxis yAxisId="depth" allowDecimals={false} />
            <YAxis
              yAxisId="freq"
              orientation="right"
              domain={[0, 1]}
              tickFormatter={(value) => `${Math.round(Number(value) * 100)}%`}
            />
            <Tooltip content={<SiteTooltip />} />
            <Legend />
            <Bar yAxisId="depth" dataKey="depth" name="Depth" fill="#79a99b" />
            <Line yAxisId="freq" dataKey="mutationFreq" name="Mutation frequency" stroke="#b85c38" strokeWidth={2} />
          </ComposedChart>
        </ResponsiveContainer>
      </section>
    </div>
  );
}

function App() {
  const [databases, setDatabases] = useState<Database[]>([]);
  const [selectedDb, setSelectedDb] = useState("");
  const [selectedGene, setSelectedGene] = useState("");
  const [active, setActive] = useState("database");
  const [status, setStatus] = useState("Ready");
  const [statusLog, setStatusLog] = useState<string[]>(["Ready"]);
  const [result, setResult] = useState<ResultPayload | null>(null);
  const [results, setResults] = useState<ResultPayload[]>([]);

  function pushStatus(message: string) {
    setStatus(message);
    setStatusLog((previous) => [message, ...previous].slice(0, 8));
  }

  async function refreshDatabases() {
    pushStatus("Refreshing database list...");
    const response = await fetch(`${API_BASE}/api/databases`);
    const payload = await response.json();
    const nextDatabases: Database[] = payload.databases || [];
    setDatabases(nextDatabases);
    if ((!selectedDb || !nextDatabases.some((db) => db.name === selectedDb)) && nextDatabases.length) {
      setSelectedDb(nextDatabases[0].name);
      setSelectedGene(nextDatabases[0].genes?.[0]?.symbol || "");
    }
    pushStatus(`Database list updated: ${nextDatabases.length} available.`);
  }

  useEffect(() => {
    refreshDatabases().catch(() => pushStatus("Backend is not reachable."));
    refreshResults().catch(() => pushStatus("Result list is not reachable."));
  }, []);

  const currentDb = databases.find((db) => db.name === selectedDb);
  const genes = currentDb?.genes || [];

  async function submitDatabase(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    pushStatus("Uploading reference and annotation files...");
    const form = new FormData(event.currentTarget);
    pushStatus("Building gene database...");
    const response = await fetch(`${API_BASE}/api/databases`, { method: "POST", body: form });
    const payload = await response.json();
    if (!response.ok) {
      pushStatus(payload.detail || "Database build failed.");
      return;
    }
    const builtDb: Database = {
      name: payload.name,
      gene_count: payload.gene_count || payload.genes?.length || 0,
      warning_count: payload.warning_count || 0,
      genes: payload.genes || [],
    };
    setDatabases((previous) => [builtDb, ...previous.filter((db) => db.name !== builtDb.name)]);
    setSelectedDb(payload.name);
    setSelectedGene(payload.genes?.[0]?.symbol || "");
    pushStatus(`Database built: ${payload.name} (${builtDb.gene_count} genes).`);
    await refreshDatabases();
  }

  async function submitAnalysis(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const formElement = event.currentTarget;
    const fastqInput = formElement.elements.namedItem("fastq") as HTMLInputElement | null;
    const bamInput = formElement.elements.namedItem("bam") as HTMLInputElement | null;
    const hasFastq = Boolean(fastqInput?.files?.length);
    const hasBam = Boolean(bamInput?.files?.length);
    if (hasFastq === hasBam) {
      pushStatus("Select exactly one input: FASTQ or BAM.");
      return;
    }
    pushStatus("Uploading sample file...");
    const form = new FormData(formElement);
    pushStatus("Running sample analysis...");
    const response = await fetch(`${API_BASE}/api/analyze`, { method: "POST", body: form });
    const payload = await response.json();
    if (!response.ok) {
      pushStatus(payload.detail || "Sample analysis failed.");
      return;
    }
    setResult(payload);
    pushStatus(`Analysis complete: ${payload.run_id}`);
    setActive("results");
  }

  async function submitSiteQuery(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = new FormData(event.currentTarget);
    pushStatus("Querying selected sites from analysis BAM...");
    const response = await fetch(`${API_BASE}/api/site-query`, { method: "POST", body: form });
    const payload = await response.json();
    if (!response.ok) {
      pushStatus(payload.detail || "Site query failed.");
      return;
    }
    setResult((previous) => ({
      run_id: payload.run_id,
      sample_summary: previous?.sample_summary || [],
      mutation_candidates: previous?.mutation_candidates || [],
      site_query: payload.rows || [],
    }));
    pushStatus(`Site query complete: ${payload.run_id} / ${payload.gene || payload.query_id}`);
    setActive("results");
  }

  async function refreshResults() {
    pushStatus("Refreshing result runs...");
    const response = await fetch(`${API_BASE}/api/results`);
    const payload = await response.json();
    const nextResults = payload.results || [];
    setResults(nextResults);
    if (nextResults.length && !result) {
      setResult(nextResults[nextResults.length - 1]);
    }
    pushStatus(`Result list updated: ${nextResults.length} available.`);
  }

  return (
    <main className="shell">
      <aside className="sidebar">
        <div>
          <h1>Mutation Patrol Robot</h1>
          <p>AMR gene database, sample analysis, and site inspection.</p>
        </div>
        <nav>
          {[
            ["database", "Gene Database"],
            ["analysis", "Sample Analysis"],
            ["query", "Site Query"],
            ["results", "Results"],
          ].map(([id, label]) => (
            <button className={active === id ? "active" : ""} key={id} onClick={() => setActive(id)}>
              {label}
            </button>
          ))}
        </nav>
        <div className="status">
          <strong>{status}</strong>
          <ul>
            {statusLog.map((item, index) => <li key={`${item}-${index}`}>{item}</li>)}
          </ul>
        </div>
      </aside>

      <section className="content">
        {active === "database" && (
          <section className="panel">
            <div className="panelHeader">
              <h2>Gene Database</h2>
              <button type="button" onClick={() => refreshDatabases()}>Refresh</button>
            </div>
            <form onSubmit={submitDatabase} className="gridForm">
              <Field label="Database name">
                <input name="name" required placeholder="pf_amr_genes" />
              </Field>
              <Field label="Reference FASTA">
                <input name="ref_fasta" type="file" required />
              </Field>
              <Field label="Gene JSONL files">
                <input name="annotations" type="file" multiple />
              </Field>
              <Field label="Gene JSONL folder">
                <FolderInput name="annotations" />
              </Field>
              <Field label="Gene symbols filter">
                <input name="genes" placeholder="optional, comma separated" />
              </Field>
              <input type="hidden" name="force" value="true" />
              <button className="primary" type="submit">Build Database</button>
            </form>

            <h3>Available Databases</h3>
            <div className="databaseList">
              {databases.map((db) => (
                <button key={db.name} className={selectedDb === db.name ? "selectedCard" : "card"} onClick={() => {
                  setSelectedDb(db.name);
                  setSelectedGene(db.genes?.[0]?.symbol || "");
                }}>
                  <strong>{db.name}</strong>
                  <span>{db.gene_count} genes · {db.warning_count} warnings</span>
                </button>
              ))}
            </div>
            {!databases.length && <div className="empty">No databases found. Build one above, or place a database under app_data/databases or databases.</div>}
            <DataTable rows={genes as unknown as Record<string, string>[]} empty="No database selected." />
          </section>
        )}

        {active === "analysis" && (
          <section className="panel">
            <h2>Sample Analysis</h2>
            <form onSubmit={submitAnalysis} className="gridForm">
              <Field label="Run name">
                <input name="run_name" required placeholder="sample_01" />
              </Field>
              <Field label="Gene database">
                <select name="db_name" value={selectedDb} onChange={(event) => setSelectedDb(event.target.value)} required>
                  <option value="">Select database</option>
                  {databases.map((db) => <option key={db.name}>{db.name}</option>)}
                </select>
              </Field>
              <Field label="FASTQ">
                <input name="fastq" type="file" />
              </Field>
              <Field label="BAM">
                <input name="bam" type="file" />
              </Field>
              <Field label="Candidate source">
                <select name="candidate_source" defaultValue="bam">
                  <option value="bam">BAM pileup</option>
                  <option value="vcf">VCF caller</option>
                </select>
              </Field>
              <Field label="Minimum depth">
                <input name="min_depth" type="number" min="0" defaultValue="10" />
              </Field>
              <Field label="Threads">
                <input name="threads" type="number" min="1" defaultValue="8" />
              </Field>
              <Field label="Minimum alt count">
                <input name="min_alt_count" type="number" min="0" defaultValue="1" />
              </Field>
              <Field label="Minimum alt frequency">
                <input name="min_alt_freq" type="number" min="0" max="1" step="0.01" defaultValue="0.05" />
              </Field>
              <input type="hidden" name="force" value="true" />
              <button className="primary" type="submit">Run Analysis</button>
            </form>
          </section>
        )}

        {active === "query" && (
          <section className="panel">
            <h2>Site Query</h2>
            <form onSubmit={submitSiteQuery} className="gridForm">
              <Field label="Analysis run">
                <select name="analysis_run" required>
                  <option value="">Select completed analysis</option>
                  {results.map((run) => <option key={run.run_id} value={run.run_id}>{run.run_id}</option>)}
                </select>
              </Field>
              <Field label="Gene database">
                <select name="db_name" value={selectedDb} onChange={(event) => setSelectedDb(event.target.value)} required>
                  <option value="">Select database</option>
                  {databases.map((db) => <option key={db.name}>{db.name}</option>)}
                </select>
              </Field>
              <Field label="Gene">
                <select name="gene" value={selectedGene} onChange={(event) => setSelectedGene(event.target.value)} required>
                  <option value="">Select gene</option>
                  {genes.map((gene) => <option key={gene.symbol} value={gene.symbol}>{gene.description || gene.symbol}</option>)}
                </select>
              </Field>
              <Field label="Position type">
                <select name="query_type" defaultValue="aa_pos">
                  <option value="genomic_pos">Genome coordinate</option>
                  <option value="cds_pos">CDS nucleotide</option>
                  <option value="aa_pos">Amino acid</option>
                </select>
              </Field>
              <Field label="Positions">
                <input name="query_value" required placeholder="436,437,540,581" />
              </Field>
              <Field label="Allele">
                <input name="alt" placeholder="optional; base, codon, or AA" />
              </Field>
              <Field label="Minimum alt frequency">
                <input name="min_alt_freq" type="number" min="0" max="1" step="0.01" defaultValue="0.05" />
              </Field>
              <input type="hidden" name="force" value="true" />
              <button type="button" onClick={() => refreshResults()}>Refresh Runs</button>
              <button className="primary" type="submit">Query Sites</button>
            </form>
          </section>
        )}

        {active === "results" && (
          <section className="panel">
            <div className="panelHeader">
              <h2>Results</h2>
              <button type="button" onClick={async () => {
                await refreshResults();
              }}>Load Latest</button>
            </div>
            <div className="resultMeta">
              <span>Sample name</span>
              <strong>{result?.run_id || "No sample loaded"}</strong>
            </div>
            <h3>Sample Summary</h3>
            <DataTable rows={result?.sample_summary || []} empty="No sample summary loaded." />
            <h3>Mutation Candidates</h3>
            <PaginatedSortableTable rows={result?.mutation_candidates || []} empty="No mutation candidates loaded." />
            <h3>Site Query Charts</h3>
            <SiteQueryCharts rows={result?.site_query || []} />
            <h3>Site Query</h3>
            <DataTable rows={result?.site_query || []} empty="No site query loaded." />
          </section>
        )}
      </section>
    </main>
  );
}

createRoot(document.getElementById("root")!).render(<App />);
