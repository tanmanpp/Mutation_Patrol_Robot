import React, { useEffect, useMemo, useRef, useState } from "react";
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
  LineChart,
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
  site_queries?: SiteQueryResult[];
  coverage_summary?: Record<string, string>[];
  coverage_regions?: Record<string, string>[];
  html_report?: string;
};

type SiteQueryResult = {
  query_id: string;
  path: string;
  updated_at: number;
  rows: Record<string, string>[];
  scan_summary?: Record<string, unknown>;
  no_call_regions?: Record<string, string>[];
  complete_table?: string;
  complete_rows?: Record<string, string>[];
};

type UploadStats = {
  total_bytes: number;
  file_count: number;
  dir_count: number;
};

type JobPayload<T = unknown> = {
  job_id: string;
  kind: string;
  label: string;
  status: "queued" | "running" | "completed" | "failed" | "cancelled";
  message: string;
  result: T | null;
  error_type?: string;
};

type SiteQueryPayload = {
  run_id: string;
  query_id: string;
  gene: string;
  query_type?: string;
};

type CoveragePayload = {
  run_id: string;
  coverage_summary: Record<string, string>[];
  coverage_regions: Record<string, string>[];
};

type CoveragePointPayload = {
  run_id: string;
  gene: string;
  minimum_call_depth: number;
  min_mapq: number;
  min_baseq: number;
  points: Record<string, string>[];
};

type IgvConfigPayload = {
  run_id: string;
  query_id: string;
  database: string;
  gene: string;
  chrom: string;
  strand: "+" | "-";
  start: number;
  end: number;
  locus: string;
  reference: {
    name: string;
    fasta_url: string;
    index_url: string;
  };
  alignment: {
    name: string;
    bam_url: string;
    index_url: string;
  };
};

const FIELD_HELP: Record<string, string> = {
  "Database name": "自訂基因庫的識別名稱。建議使用英文字母、數字、連字號或底線，例如 pf3d7_targets。",
  "Reference FASTA": "上傳與 BAM 或 FASTQ 分析所用版本相同的參考基因組 FASTA；contig 名稱與座標必須一致。",
  "Gene JSONL files": "上傳 NCBI Datasets Gene 的 Product report JSONL。可以一次選擇多個基因註解檔。",
  "Gene JSONL folder": "選擇包含多個 .jsonl 基因註解檔的資料夾；系統只會讀取其中的 JSONL 檔案。",
  "Gene symbols filter": "選填。只建立指定的 gene symbols，多個名稱請用逗號分隔；留白表示匯入全部上傳的註解。",
  "Run name": "這次樣本分析的識別名稱，也會成為結果資料夾名稱；請使用容易辨識且不與既有 run 重複的名稱。",
  "Gene database": "選擇本次分析或查詢使用的參考基因庫；它應與樣本比對所用的參考基因組版本一致。",
  FASTQ: "上傳一個 FASTQ 或 FASTQ.GZ 檔。FASTQ 與 BAM 必須擇一上傳，不可同時使用。",
  BAM: "上傳一個 BAM 檔。系統會進行排序並建立 index；BAM 的 contig 名稱必須與所選基因庫相符。",
  "Candidate source": "BAM pileup 會保留深度、allele frequency 與 MIXED_SIGNAL 資訊；VCF caller 適合標準變異呼叫，但可用的 allele QC 資訊可能較少。",
  "Minimum depth": "一個位置至少需要多少個合格 reads 才可分析；低於此值的點位會視為 NO_CALL／無法判斷。",
  Threads: "分析可使用的 CPU 執行緒數。數值較高通常較快，但也會使用更多電腦資源。",
  "Minimum alt count": "替代 allele 至少需要多少個 reads 支持才會列為候選；設為 1 對低頻訊號較敏感，也較容易納入雜訊。",
  "Minimum base quality": "只計算 Phred base quality 達到此門檻的鹼基；20 約代表 1% 的鹼基判讀錯誤率。",
  "Minimum alt frequency": "替代 allele 占合格 reads 的最低比例。0.05 代表 5%；降低門檻會增加低頻訊號，也可能增加雜訊。",
  "Analysis run": "選擇一個已完成、且仍保留 merged.sorted.bam 的分析結果，作為這次查詢或 coverage 計算的資料來源。",
  Gene: "選擇要掃描完整註解區段或查詢指定點位的基因。",
  "Position type": "Whole annotated gene 會掃描 gene_start 到 gene_end 的整個註解區段；其他模式則分別使用 genomic、CDS 或胺基酸座標。",
  Positions: "只有指定點位模式需要填寫。可輸入多個 1-based 位置並用逗號分隔，例如 436,437,540。",
  Allele: "選填，僅用於指定點位模式。可輸入想確認的鹼基、codon 或胺基酸；完整基因掃描會自動尋找所有合格變異。",
  "Minimum depth for a call": "點位至少需要多少個合格 reads 才能判讀；低於此值會回報 NO_CALL／無法判斷，而不是野生型。",
  "Minimum mapping quality": "只使用 mapping quality 達到此門檻的 reads。20 常作為基本篩選；提高門檻會更嚴格。",
  "Minimum allele support": "某個 allele 至少需要多少個 reads 支持，才能參與 VARIANT 或 MIXED_SIGNAL 判定；預設為 2。",
  "Minimum call depth": "Coverage 中判定某位置為 CALLABLE 的最低深度；低於此值的區域會列為 NO_CALL／無法判斷。",
};

function Field(props: { label: string; help?: string; children: React.ReactNode }) {
  const help = props.help ?? FIELD_HELP[props.label];
  return (
    <label className="field">
      <span className="fieldLabel">
        {props.label}
        {help && (
          <span
            className="fieldHelpIcon"
            aria-label={`${props.label} 說明：${help}`}
            tabIndex={0}
          >
            ?
          </span>
        )}
      </span>
      {props.children}
      {help && (
        <span className="fieldTooltip" role="tooltip">
          {help}
        </span>
      )}
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
  const columns = useMemo(
    () => (
      rows[0]
        ? Object.keys(rows[0]).filter((column) => rows.some((row) => row[column] !== ""))
        : []
    ),
    [rows],
  );
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
              {columns.map((column) => (
                <td key={column}>{displayTableValue(column, row[column])}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function formatBytes(bytes: number) {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value >= 10 || unitIndex === 0 ? value.toFixed(0) : value.toFixed(1)} ${units[unitIndex]}`;
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
  const columns = useMemo(
    () => (
      rows[0]
        ? Object.keys(rows[0]).filter((column) => rows.some((row) => row[column] !== ""))
        : []
    ),
    [rows],
  );
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
                {columns.map((column) => (
                  <td key={column}>{displayTableValue(column, row[column])}</td>
                ))}
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
  callStatus: string;
  alleleSpectrum: string;
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

function optionalNumber(value: string | undefined) {
  if (value === undefined || value === "") return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function compareBiologicalPositions(
  left: Record<string, string>,
  right: Record<string, string>,
) {
  const leftAa = optionalNumber(left.aa_pos);
  const rightAa = optionalNumber(right.aa_pos);
  if (leftAa !== null || rightAa !== null) {
    if (leftAa === null) return 1;
    if (rightAa === null) return -1;
    if (leftAa !== rightAa) return leftAa - rightAa;
  }

  const leftCds = optionalNumber(left.cds_pos);
  const rightCds = optionalNumber(right.cds_pos);
  if (leftCds !== null || rightCds !== null) {
    if (leftCds === null) return 1;
    if (rightCds === null) return -1;
    if (leftCds !== rightCds) return leftCds - rightCds;
  }

  const leftGenomic = optionalNumber(left.pos || left.genomic_pos);
  const rightGenomic = optionalNumber(right.pos || right.genomic_pos);
  if (leftGenomic !== null && rightGenomic !== null && leftGenomic !== rightGenomic) {
    return leftGenomic - rightGenomic;
  }
  return compareCellValues(left.alt || left.alt_codon, right.alt || right.alt_codon);
}

function sortByBiologicalPosition(rows: Record<string, string>[]) {
  return [...rows].sort(compareBiologicalPositions);
}

function displayTableValue(column: string, value: string | undefined) {
  if (column === "call_status" && value === "NO_CALL") {
    return "無法判斷 (NO_CALL)";
  }
  if (column === "call_status" && value === "MIXED_SIGNAL") {
    return "多 allele 訊號 (MIXED_SIGNAL)";
  }
  if (column === "call_status" && value === "REFERENCE") {
    return "與參考序列一致 (REFERENCE)";
  }
  if (column === "call_status" && value === "VARIANT") {
    return "觀察到變異 (VARIANT)";
  }
  if (column === "filter" && value === "LOW_DEPTH") {
    return "深度不足 (LOW_DEPTH)";
  }
  if (column === "filter" && value === "NO_MAPPING") {
    return "沒有可用 reads (NO_MAPPING)";
  }
  return value || "";
}

function siteQueryLabel(row: Record<string, string>) {
  const gene = row["Gene symbol"] || row.gene || "";
  const geneName = row["Gene name"] || gene;
  const aaPos = row.aa_pos || "";
  const aaChange = row.aa_change && row.aa_change !== "NA" ? row.aa_change : "";
  const position = aaPos || row.pos || "";
  const allele = row.alt_codon || row.alt || "";
  if (aaChange) return `${geneName} ${aaChange}`;
  if (position && allele) return `${geneName} ${position} ${allele}`;
  return `${geneName} ${position}`.trim();
}

function siteQueryChartData(rows: Record<string, string>[]): SiteChartDatum[] {
  return sortByBiologicalPosition(rows).map((row, index) => ({
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
    callStatus: row.call_status || "",
    alleleSpectrum: row.allele_spectrum || "",
  }));
}

function sitePositionKey(item: SiteChartDatum) {
  return `${item.gene}|${item.position}`;
}

function sitePositionLabel(item: SiteChartDatum) {
  const geneName = item.geneDescription || item.gene;
  return geneName && item.position ? `${geneName} ${item.position}` : item.label;
}

function siteQueryDisplayName(query: SiteQueryResult) {
  const firstRow = query.rows?.[0];
  const geneName = firstRow?.["Gene name"];
  if (query.scan_summary?.scan_type === "whole_gene") {
    return `${String(query.scan_summary.gene_name || geneName || query.query_id)} — whole gene`;
  }
  return geneName || query.query_id;
}

function scanSummaryRows(summary: Record<string, unknown> | undefined) {
  if (!summary) return [];
  const preferred = [
    "gene", "gene_name", "chrom", "start", "end", "gene_length",
    "callable_positions", "callable_percent", "no_call_positions",
    "no_call_regions", "reference_only_positions", "variant_sites",
    "variant_rows", "complete_table_rows", "mixed_signal_sites",
    "variant_type_counts",
  ];
  return preferred
    .filter((key) => summary[key] !== undefined)
    .map((key) => ({
      metric: key,
      value: (
        typeof summary[key] === "object"
          ? JSON.stringify(summary[key])
          : String(summary[key])
      ),
    }));
}

function aminoAcidDetailRows(rows: Record<string, string>[]) {
  return sortByBiologicalPosition(rows)
    .filter((row) => row.region_type === "CDS")
    .map((row) => ({
      gene: row["Gene symbol"] || row.gene || "",
      genomic_pos: row.pos || "",
      variant_type: row.variant_type || "",
      nucleotide_change: `${row.ref || ""}>${row.alt || ""}`,
      cds_pos: row.cds_pos || "",
      codon_pos: row.codon_pos || "",
      codon_change: row.codon_change || "",
      aa_pos: row.aa_pos || "",
      ref_aa: row.ref_aa || "",
      ref_aa_name: row.ref_aa_name || "",
      alt_aa: row.alt_aa || "",
      alt_aa_name: row.alt_aa_name || "",
      aa_change: row.aa_change || "",
      effect: row.effect || "",
      allele_freq: row.allele_freq || "",
      call_status: row.call_status || "",
    }));
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
    if (item.variantType !== "REF" && item.variantType !== "NO_CALL") {
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
  if (item.callStatus === "NO_CALL") return "#a3adb0";
  if (item.callStatus === "MIXED_SIGNAL") return "#d38b2f";
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
      {item.callStatus && <span>Call: {displayTableValue("call_status", item.callStatus)}</span>}
      {item.alleleSpectrum && <span>Alleles: {item.alleleSpectrum}</span>}
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

type CoverageChartPoint = {
  position: number;
  depth: number;
  minimumDepth: number;
};

function coverageChartData(
  points: Record<string, string>[],
  minimumDepth: number,
): CoverageChartPoint[] {
  if (!points.length) return [];
  const maxPoints = 1200;
  const step = Math.max(1, Math.ceil(points.length / maxPoints));
  const output: CoverageChartPoint[] = [];
  for (let index = 0; index < points.length; index += step) {
    const bucket = points.slice(index, index + step);
    const totalDepth = bucket.reduce((sum, point) => sum + parseNumber(point.depth), 0);
    output.push({
      position: parseNumber(bucket[Math.floor(bucket.length / 2)]?.pos),
      depth: totalDepth / bucket.length,
      minimumDepth,
    });
  }
  return output;
}

function CoverageChart({
  payload,
}: {
  payload: CoveragePointPayload | null;
}) {
  const data = useMemo(
    () => coverageChartData(payload?.points || [], payload?.minimum_call_depth || 10),
    [payload],
  );
  if (!payload || !data.length) {
    return <div className="empty">Select a gene to view position-level coverage.</div>;
  }
  return (
    <section className="chartPanel chartPanelWide">
      <h3>{payload.gene} Position Coverage</h3>
      <p className="helperText">
        Depth is displayed across the selected gene or amplicon-aligned region.
        This is coverage evidence only and is not a CNV estimate.
      </p>
      <ResponsiveContainer width="100%" height={330}>
        <LineChart data={data} margin={{ top: 14, right: 24, bottom: 42, left: 0 }}>
          <CartesianGrid strokeDasharray="3 3" vertical={false} />
          <XAxis
            dataKey="position"
            type="number"
            domain={["dataMin", "dataMax"]}
            tick={{ fontSize: 11 }}
          />
          <YAxis allowDecimals={false} />
          <Tooltip />
          <Legend />
          <Line
            dataKey="depth"
            name="Mean depth"
            stroke="#2f6f9f"
            dot={false}
            strokeWidth={1.5}
          />
          <Line
            dataKey="minimumDepth"
            name="Minimum call depth"
            stroke="#b85c38"
            dot={false}
            strokeDasharray="6 4"
          />
        </LineChart>
      </ResponsiveContainer>
    </section>
  );
}

function apiResourceUrl(path: string) {
  if (/^https?:\/\//i.test(path)) return path;
  return `${API_BASE}${path}`;
}

function IgvViewer({
  config,
  onError,
}: {
  config: IgvConfigPayload;
  onError: (message: string) => void;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    let browser: import("igv").Browser | null = null;
    let igvModule: typeof import("igv") | null = null;

    async function createViewer() {
      if (!containerRef.current) return;
      setLoading(true);
      try {
        igvModule = await import("igv");
        const options = {
          reference: {
            id: config.database,
            name: config.reference.name,
            fastaURL: apiResourceUrl(config.reference.fasta_url),
            indexURL: apiResourceUrl(config.reference.index_url),
          },
          locus: config.locus,
          showNavigation: true,
          tracks: [
            {
              name: "Reference sequence — three-frame translation",
              type: "sequence",
              frameTranslate: true,
              reversed: config.strand === "-",
              removable: false,
            },
            {
              name: config.alignment.name,
              type: "alignment",
              format: "bam",
              url: apiResourceUrl(config.alignment.bam_url),
              indexURL: apiResourceUrl(config.alignment.index_url),
              height: 500,
              displayMode: "EXPANDED",
            },
          ],
        } as unknown as import("igv").CreateOpt;
        const created = await igvModule.default.createBrowser(
          containerRef.current,
          options,
        );
        if (cancelled) {
          igvModule.default.removeBrowser(created);
          return;
        }
        browser = created;
        setLoading(false);
      } catch (error) {
        if (!cancelled) {
          setLoading(false);
          onError(
            error instanceof Error
              ? error.message
              : "IGV could not be initialized.",
          );
        }
      }
    }

    createViewer();
    return () => {
      cancelled = true;
      if (browser && igvModule) {
        igvModule.default.removeBrowser(browser);
      }
    };
  }, [config, onError]);

  return (
    <section className="igvPanel" id="igv-viewer-panel">
      <div className="sectionHeader">
        <div>
          <h3>IGV Genome Viewer</h3>
          <p className="igvLocus">
            {config.gene}: {config.locus}
          </p>
        </div>
      </div>
      <div className="infoBox">
        Reference FASTA and the analysis BAM are loaded automatically. The
        initial view is centered on the complete annotated gene interval.
        Three-frame translation is enabled on the reference sequence track and
        follows the annotated gene strand.
        Alignment mismatches and indels are visual evidence and are not an
        automatic clinical interpretation.
      </div>
      {loading && <div className="empty">Loading reference and alignments into IGV...</div>}
      <div className="igvContainer" ref={containerRef} />
    </section>
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
  const [selectedResultId, setSelectedResultId] = useState("");
  const [selectedSiteQueryId, setSelectedSiteQueryId] = useState("");
  const [siteQueryType, setSiteQueryType] = useState("gene_region");
  const [coverageGene, setCoverageGene] = useState("");
  const [coveragePoints, setCoveragePoints] = useState<CoveragePointPayload | null>(null);
  const [igvConfig, setIgvConfig] = useState<IgvConfigPayload | null>(null);
  const [igvLoading, setIgvLoading] = useState(false);
  const [igvError, setIgvError] = useState("");
  const [uploadStats, setUploadStats] = useState<UploadStats>({ total_bytes: 0, file_count: 0, dir_count: 0 });

  function pushStatus(message: string) {
    setStatus(message);
    setStatusLog((previous) => [message, ...previous].slice(0, 8));
  }

  const reportIgvError = React.useCallback((message: string) => {
    setIgvError(message);
    setStatus(`IGV error: ${message}`);
    setStatusLog((previous) => [`IGV error: ${message}`, ...previous].slice(0, 8));
  }, []);

  async function checkRuntimeHealth() {
    const response = await fetch(`${API_BASE}/api/health`);
    const payload = await response.json();
    if (!response.ok || !payload.api_ok) {
      throw new Error(payload.detail || "Backend health check failed.");
    }
    if (!payload.ok) {
      const missing = Object.entries(payload.tools || {})
        .filter(([, value]) => !(value as { ok?: boolean }).ok)
        .map(([name]) => name);
      throw new Error(`Environment is missing: ${missing.join(", ")}`);
    }
    pushStatus("Environment ready.");
  }

  async function waitForJob<T>(initialJob: JobPayload<T>, description: string): Promise<T> {
    let job = initialJob;
    while (job.status === "queued" || job.status === "running") {
      pushStatus(`${description}: ${job.status}...`);
      await new Promise((resolve) => window.setTimeout(resolve, 1000));
      const response = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(job.job_id)}`);
      const payload = await response.json();
      if (!response.ok) {
        throw new Error(payload.detail || `Could not read job ${job.job_id}.`);
      }
      job = payload;
    }
    if (job.status === "completed" && job.result) {
      return job.result;
    }
    throw new Error(job.message || `${description} did not complete.`);
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
    checkRuntimeHealth()
      .then(() => Promise.all([
        refreshDatabases(),
        refreshResults(),
        refreshUploadStats(),
      ]))
      .catch((error) => pushStatus(
        error instanceof Error ? error.message : "Backend is not reachable."
      ));
  }, []);

  const currentDb = databases.find((db) => db.name === selectedDb);
  const genes = currentDb?.genes || [];
  const selectedSiteQuery = useMemo(() => {
    const queries = result?.site_queries || [];
    if (!queries.length) return null;
    return queries.find((query) => query.query_id === selectedSiteQueryId) || queries[queries.length - 1];
  }, [result, selectedSiteQueryId]);
  const selectedSiteQueryRows = useMemo(
    () => sortByBiologicalPosition(
      selectedSiteQuery?.rows || result?.site_query || [],
    ),
    [selectedSiteQuery?.rows, result?.site_query],
  );
  const selectedScanSummaryRows = scanSummaryRows(selectedSiteQuery?.scan_summary);
  const selectedNoCallRegions = selectedSiteQuery?.no_call_regions || [];
  const selectedCompleteRows = selectedSiteQuery?.complete_rows || [];
  const selectedAminoAcidRows = aminoAcidDetailRows(selectedSiteQueryRows);
  const selectedCoverageRegions = (result?.coverage_regions || []).filter(
    (row) => row["Gene symbol"] === coverageGene,
  );

  useEffect(() => {
    if (active === "coverage" && result?.run_id && coverageGene) {
      loadCoverageGene(result.run_id, coverageGene);
    }
  }, [active, result?.run_id, coverageGene]);

  useEffect(() => {
    setIgvConfig(null);
    setIgvError("");
  }, [result?.run_id, selectedSiteQuery?.query_id]);

  function setLoadedResult(payload: ResultPayload) {
    setResult(payload);
    setSelectedResultId(payload.run_id);
    const queries = payload.site_queries || [];
    setSelectedSiteQueryId(queries[queries.length - 1]?.query_id || "");
    const coverageGenes = payload.coverage_summary || [];
    setCoverageGene((previous) => (
      coverageGenes.some((row) => row["Gene symbol"] === previous)
        ? previous
        : coverageGenes[0]?.["Gene symbol"] || ""
    ));
    setCoveragePoints(null);
  }

  async function submitDatabase(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      pushStatus("Uploading reference and annotation files...");
      const form = new FormData(event.currentTarget);
      const response = await fetch(`${API_BASE}/api/databases`, { method: "POST", body: form });
      const job = await response.json();
      if (!response.ok) {
        pushStatus(job.detail || "Database build failed.");
        return;
      }
      const payload = await waitForJob<Database & { name: string }>(job, "Building gene database");
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
    } catch (error) {
      pushStatus(error instanceof Error ? error.message : "Database build failed.");
    }
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
    try {
      pushStatus("Uploading sample file...");
      const form = new FormData(formElement);
      const response = await fetch(`${API_BASE}/api/analyze`, { method: "POST", body: form });
      const job = await response.json();
      if (!response.ok) {
        pushStatus(job.detail || "Sample analysis failed.");
        return;
      }
      const payload = await waitForJob<ResultPayload>(job, "Running sample analysis");
      setLoadedResult(payload);
      await refreshResults(payload.run_id);
      pushStatus(`Analysis complete: ${payload.run_id}`);
      setActive("results");
    } catch (error) {
      pushStatus(error instanceof Error ? error.message : "Sample analysis failed.");
    }
  }

  async function clearSampleUploads() {
    const confirmed = window.confirm(
      "Clear uploaded sample files? Existing analysis results will not be deleted."
    );
    if (!confirmed) return;

    pushStatus("Clearing uploaded sample files...");
    const response = await fetch(`${API_BASE}/api/uploads/samples`, { method: "DELETE" });
    const payload = await response.json();
    if (!response.ok) {
      pushStatus(payload.detail || "Failed to clear uploaded sample files.");
      return;
    }
    setUploadStats({
      total_bytes: payload.stats?.total_bytes || 0,
      file_count: payload.stats?.file_count || 0,
      dir_count: payload.stats?.dir_count || 0,
    });
    pushStatus(`Uploaded sample files cleared: ${payload.deleted_files || 0} files, ${formatBytes(payload.deleted_bytes || 0)} removed.`);
  }

  async function submitSiteQuery(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      const form = new FormData(event.currentTarget);
      const wholeGene = form.get("query_type") === "gene_region";
      pushStatus(wholeGene ? "Starting whole-gene scan..." : "Starting site query...");
      const response = await fetch(`${API_BASE}/api/site-query`, { method: "POST", body: form });
      const job = await response.json();
      if (!response.ok) {
        pushStatus(job.detail || "Site query failed.");
        return;
      }
      const payload = await waitForJob<SiteQueryPayload>(
        job,
        wholeGene ? "Scanning complete gene interval" : "Querying selected sites",
      );
      await loadResult(payload.run_id);
      pushStatus(
        `${wholeGene ? "Whole-gene scan" : "Site query"} complete: `
        + `${payload.run_id} / ${payload.gene || payload.query_id}`,
      );
      setActive("results");
    } catch (error) {
      pushStatus(error instanceof Error ? error.message : "Site query failed.");
    }
  }

  async function submitCoverage(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      const form = new FormData(event.currentTarget);
      pushStatus("Generating gene and region coverage...");
      const response = await fetch(`${API_BASE}/api/coverage`, { method: "POST", body: form });
      const job = await response.json();
      if (!response.ok) {
        pushStatus(job.detail || "Coverage generation failed.");
        return;
      }
      const payload = await waitForJob<CoveragePayload>(job, "Calculating coverage");
      await loadResult(payload.run_id);
      const firstGene = payload.coverage_summary?.[0]?.["Gene symbol"] || "";
      if (firstGene) await loadCoverageGene(payload.run_id, firstGene);
      pushStatus(`Coverage ready: ${payload.run_id}`);
    } catch (error) {
      pushStatus(error instanceof Error ? error.message : "Coverage generation failed.");
    }
  }

  async function loadCoverageGene(runId: string, gene: string) {
    setCoverageGene(gene);
    if (!runId || !gene) {
      setCoveragePoints(null);
      return;
    }
    const response = await fetch(
      `${API_BASE}/api/results/${encodeURIComponent(runId)}/coverage/${encodeURIComponent(gene)}`,
    );
    const payload = await response.json();
    if (!response.ok) {
      setCoveragePoints(null);
      pushStatus(payload.detail || `Could not load coverage for ${gene}.`);
      return;
    }
    setCoveragePoints(payload);
    pushStatus(`Coverage loaded: ${gene}`);
  }

  async function refreshResults(preferredRunId = selectedResultId) {
    pushStatus("Refreshing result runs...");
    const response = await fetch(`${API_BASE}/api/results`);
    const payload = await response.json();
    const nextResults = payload.results || [];
    setResults(nextResults);
    const nextSelection = preferredRunId && nextResults.some((run: ResultPayload) => run.run_id === preferredRunId)
      ? preferredRunId
      : nextResults[nextResults.length - 1]?.run_id || "";
    if (nextSelection) {
      setSelectedResultId(nextSelection);
      await loadResult(nextSelection);
    } else {
      setSelectedResultId("");
      setSelectedSiteQueryId("");
      setResult(null);
    }
    pushStatus(`Result list updated: ${nextResults.length} available.`);
  }

  async function loadResult(runId: string) {
    if (!runId) {
      setSelectedResultId("");
      setResult(null);
      return;
    }
    pushStatus(`Loading result: ${runId}...`);
    const response = await fetch(`${API_BASE}/api/results/${encodeURIComponent(runId)}`);
    const payload = await response.json();
    if (!response.ok) {
      pushStatus(payload.detail || `Failed to load result: ${runId}`);
      return;
    }
    setLoadedResult(payload);
    pushStatus(`Result loaded: ${payload.run_id}`);
  }

  function openHtmlReport() {
    if (!result?.run_id) {
      pushStatus("No result selected for HTML report.");
      return;
    }
    window.open(`${API_BASE}/api/results/${encodeURIComponent(result.run_id)}/report`, "_blank");
    pushStatus(`HTML report requested: ${result.run_id}`);
  }

  function downloadCompleteGeneTable() {
    if (!result?.run_id || !selectedSiteQuery?.query_id || !selectedSiteQuery.complete_table) {
      pushStatus("No complete whole-gene table is available.");
      return;
    }
    window.open(
      `${API_BASE}/api/results/${encodeURIComponent(result.run_id)}`
      + `/site-query/${encodeURIComponent(selectedSiteQuery.query_id)}/complete-table`,
      "_blank",
    );
    pushStatus(`Complete gene table requested: ${selectedSiteQuery.query_id}`);
  }

  async function openIgvViewer() {
    if (!result?.run_id || !selectedSiteQuery?.query_id) {
      pushStatus("Select a Whole Gene Scan result before opening IGV.");
      return;
    }
    setIgvLoading(true);
    setIgvError("");
    try {
      const response = await fetch(
        `${API_BASE}/api/results/${encodeURIComponent(result.run_id)}`
        + `/site-query/${encodeURIComponent(selectedSiteQuery.query_id)}/igv`,
      );
      const payload = await response.json();
      if (!response.ok) {
        throw new Error(payload.detail || "IGV configuration could not be loaded.");
      }
      setIgvConfig(payload);
      pushStatus(`IGV ready: ${payload.gene} / ${payload.locus}`);
      window.setTimeout(() => {
        document.getElementById("igv-viewer-panel")?.scrollIntoView({
          behavior: "smooth",
          block: "start",
        });
      }, 100);
    } catch (error) {
      const message = error instanceof Error ? error.message : "IGV could not be opened.";
      setIgvError(message);
      pushStatus(message);
    } finally {
      setIgvLoading(false);
    }
  }

  async function refreshUploadStats() {
    const response = await fetch(`${API_BASE}/api/uploads/samples`);
    const payload = await response.json();
    if (!response.ok) return;
    setUploadStats({
      total_bytes: payload.total_bytes || 0,
      file_count: payload.file_count || 0,
      dir_count: payload.dir_count || 0,
    });
  }

  return (
    <main className="shell">
      <aside className="sidebar">
        <div>
          <h1>Mutation Patrol Robot</h1>
          <p>Configurable gene database, mutation browsing, and coverage inspection.</p>
        </div>
        <nav>
          {[
            ["database", "Gene Database"],
            ["analysis", "Sample Analysis"],
            ["query", "Whole Gene Scan"],
            ["coverage", "Coverage"],
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
            <div className="panelHeader">
              <h2>Sample Analysis</h2>
              <button
                type="button"
                className="dangerButton"
                title={`${uploadStats.file_count} uploaded sample files`}
                onClick={clearSampleUploads}
              >
                Clear Sample Uploads ({formatBytes(uploadStats.total_bytes)})
              </button>
            </div>
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
              <Field label="Minimum base quality">
                <input name="min_baseq" type="number" min="0" max="93" defaultValue="20" />
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
            <h2>Whole Gene Scan / Site Query</h2>
            <div className="infoBox">
              Whole Gene Scan examines the complete annotated genomic interval
              of the selected gene and reports supported SNV/indel rows plus
              NO_CALL region summaries. MIXED_SIGNAL describes multiple
              supported alleles only; it does not diagnose mixed infection.
            </div>
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
                <select
                  name="query_type"
                  value={siteQueryType}
                  onChange={(event) => setSiteQueryType(event.target.value)}
                >
                  <option value="gene_region">Whole annotated gene</option>
                  <option value="genomic_pos">Genome coordinate</option>
                  <option value="cds_pos">CDS nucleotide</option>
                  <option value="aa_pos">Amino acid</option>
                </select>
              </Field>
              <Field label="Positions">
                <input
                  name="query_value"
                  required={siteQueryType !== "gene_region"}
                  disabled={siteQueryType === "gene_region"}
                  placeholder={
                    siteQueryType === "gene_region"
                      ? "Not required for whole-gene scan"
                      : "436,437,540,581"
                  }
                />
              </Field>
              <Field label="Allele">
                <input
                  name="alt"
                  disabled={siteQueryType === "gene_region"}
                  placeholder={
                    siteQueryType === "gene_region"
                      ? "All supported variants will be reported"
                      : "optional; base, codon, or AA"
                  }
                />
              </Field>
              <Field label="Minimum alt frequency">
                <input name="min_alt_freq" type="number" min="0" max="1" step="0.01" defaultValue="0.05" />
              </Field>
              <Field label="Minimum depth for a call">
                <input name="min_depth" type="number" min="1" step="1" defaultValue="10" />
              </Field>
              <Field label="Minimum mapping quality">
                <input name="min_mapq" type="number" min="0" max="255" step="1" defaultValue="20" />
              </Field>
              <Field label="Minimum base quality">
                <input name="min_baseq" type="number" min="0" max="93" step="1" defaultValue="20" />
              </Field>
              <Field label="Minimum allele support">
                <input name="min_allele_count" type="number" min="1" step="1" defaultValue="2" />
              </Field>
              <input type="hidden" name="force" value="true" />
              <button type="button" onClick={() => refreshResults()}>Refresh Runs</button>
              <button className="primary" type="submit">
                {siteQueryType === "gene_region" ? "Scan Whole Gene" : "Query Sites"}
              </button>
            </form>
          </section>
        )}

        {active === "coverage" && (
          <section className="panel">
            <div className="panelHeader">
              <h2>Gene / Amplicon Coverage</h2>
              <button type="button" onClick={() => refreshResults()}>Refresh Runs</button>
            </div>
            <div className="infoBox">
              This page displays observed read depth and callable positions only.
              It does not calculate or infer CNV. When explicit amplicon regions
              are not present in the database, exon/target regions are shown.
            </div>
            <form onSubmit={submitCoverage} className="gridForm">
              <Field label="Analysis run">
                <select
                  name="analysis_run"
                  value={selectedResultId}
                  onChange={(event) => loadResult(event.target.value)}
                  required
                >
                  <option value="">Select completed analysis</option>
                  {results.map((run) => (
                    <option key={run.run_id} value={run.run_id}>{run.run_id}</option>
                  ))}
                </select>
              </Field>
              <Field label="Gene database">
                <select
                  name="db_name"
                  value={selectedDb}
                  onChange={(event) => setSelectedDb(event.target.value)}
                  required
                >
                  <option value="">Select database</option>
                  {databases.map((db) => <option key={db.name}>{db.name}</option>)}
                </select>
              </Field>
              <Field label="Minimum call depth">
                <input name="minimum_call_depth" type="number" min="1" defaultValue="10" />
              </Field>
              <Field label="Minimum mapping quality">
                <input name="min_mapq" type="number" min="0" max="255" defaultValue="20" />
              </Field>
              <Field label="Minimum base quality">
                <input name="min_baseq" type="number" min="0" max="93" defaultValue="20" />
              </Field>
              <input type="hidden" name="force" value="true" />
              <button className="primary" type="submit">Generate Coverage</button>
            </form>

            <h3>Gene Coverage Summary</h3>
            <PaginatedSortableTable
              rows={result?.coverage_summary || []}
              empty="No coverage is available. Select a run and generate coverage."
            />

            <div className="sectionHeader">
              <h3>Position Coverage</h3>
              <select
                value={coverageGene}
                onChange={(event) => setCoverageGene(event.target.value)}
              >
                <option value="">Select gene</option>
                {(result?.coverage_summary || []).map((row) => (
                  <option key={row["Gene symbol"]} value={row["Gene symbol"]}>
                    {row["Gene name"] || row["Gene symbol"]}
                  </option>
                ))}
              </select>
            </div>
            <CoverageChart payload={coveragePoints} />

            <h3>Amplicon / Target-region Coverage</h3>
            <PaginatedSortableTable
              rows={selectedCoverageRegions}
              empty="No amplicon or target-region coverage is available for this gene."
            />
          </section>
        )}

        {active === "results" && (
          <section className="panel">
            <div className="panelHeader">
              <h2>Results</h2>
              <div className="resultControls">
                <select
                  value={selectedResultId}
                  onChange={(event) => loadResult(event.target.value)}
                >
                  <option value="">Select result</option>
                  {results.map((run) => <option key={run.run_id} value={run.run_id}>{run.run_id}</option>)}
                </select>
                <button type="button" onClick={() => refreshResults()}>
                  Refresh
                </button>
                <button type="button" onClick={() => {
                  const latest = results[results.length - 1];
                  if (latest) loadResult(latest.run_id);
                }}>
                  Load Latest
                </button>
                <button type="button" onClick={openHtmlReport} disabled={!result?.run_id}>
                  Export HTML Report
                </button>
              </div>
            </div>
            <div className="resultMeta">
              <span>Sample name</span>
              <strong>{result?.run_id || "No sample loaded"}</strong>
            </div>
            <h3>Sample Summary</h3>
            <DataTable rows={result?.sample_summary || []} empty="No sample summary loaded." />
            <h3>Mutation Candidates</h3>
            <PaginatedSortableTable rows={result?.mutation_candidates || []} empty="No mutation candidates loaded." />
            <div className="sectionHeader">
              <h3>Gene Scan / Site Query</h3>
              <div className="sectionControls">
                <select
                  value={selectedSiteQuery?.query_id || ""}
                  onChange={(event) => setSelectedSiteQueryId(event.target.value)}
                >
                  <option value="">Select site query</option>
                  {(result?.site_queries || []).map((query) => (
                    <option key={query.query_id} value={query.query_id}>{siteQueryDisplayName(query)}</option>
                  ))}
                </select>
                {selectedSiteQuery?.scan_summary?.scan_type === "whole_gene" && (
                  <button
                    className="primary"
                    type="button"
                    onClick={openIgvViewer}
                    disabled={igvLoading}
                  >
                    {igvLoading ? "Preparing IGV..." : "Open in IGV"}
                  </button>
                )}
              </div>
            </div>
            {igvError && <div className="errorBox">{igvError}</div>}
            {igvConfig && (
              <IgvViewer
                config={igvConfig}
                onError={reportIgvError}
              />
            )}
            {selectedScanSummaryRows.length > 0 && (
              <>
                <h3>Whole Gene Scan Summary</h3>
                <DataTable
                  rows={selectedScanSummaryRows}
                  empty="No whole-gene scan summary loaded."
                />
                <h3>NO_CALL Regions</h3>
                <PaginatedSortableTable
                  rows={selectedNoCallRegions}
                  empty="No NO_CALL regions were found."
                />
              </>
            )}
            <SiteQueryCharts rows={selectedSiteQueryRows} />
            <h3>Detected Variant Alleles</h3>
            <PaginatedSortableTable
              rows={selectedSiteQueryRows}
              empty="No supported variants were detected for this query."
            />
            <h3>Coding and Amino-acid Details</h3>
            <PaginatedSortableTable
              rows={selectedAminoAcidRows}
              empty="No coding-region amino-acid changes were detected."
            />
            {selectedSiteQuery?.complete_table && (
              <>
                <div className="sectionHeader">
                  <h3>Complete Gene Table</h3>
                  <button type="button" onClick={downloadCompleteGeneTable}>
                    Export Complete CSV
                  </button>
                </div>
                <div className="infoBox">
                  The preview shows up to 1,000 rows. Export Complete CSV
                  contains every genomic position and all reported allele rows,
                  including REFERENCE, VARIANT, MIXED_SIGNAL, and NO_CALL.
                </div>
                <PaginatedSortableTable
                  rows={selectedCompleteRows}
                  empty="No complete gene table preview is available."
                />
              </>
            )}
          </section>
        )}
      </section>
    </main>
  );
}

createRoot(document.getElementById("root")!).render(<App />);
