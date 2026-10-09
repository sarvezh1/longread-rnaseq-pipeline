#!/usr/bin/env python3
"""Build an offline Block 6 report and publication source tables from manifests."""

from __future__ import annotations

import argparse, csv, html, json, shutil
from collections import Counter
from pathlib import Path


def rows(path: Path):
    with path.open(newline="",encoding="utf-8") as handle: return list(csv.DictReader(handle,delimiter="\t"))


def write_tsv(path: Path, records: list[dict]):
    fields=list(records[0]) if records else []
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(records)


def table(records: list[dict], fields: list[str]) -> str:
    if not records: return '<p class="status">No records.</p>'
    header="".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body="".join("<tr>"+"".join(f"<td>{html.escape(str(row.get(field,'')))}</td>" for field in fields)+"</tr>" for row in records)
    return f"<div class=scroll><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>"


def status_svg(counts: Counter, destination: Path):
    labels=sorted(counts); maximum=max(counts.values(),default=1); width=820; height=90+55*len(labels)
    bars=[]
    for index,label in enumerate(labels):
        y=55+index*55; bar=int(540*counts[label]/maximum)
        bars.append(f'<text x="15" y="{y+20}" font-size="14">{html.escape(label)}</text><rect x="245" y="{y}" width="{bar}" height="28" fill="#356a8a"/><text x="{250+bar}" y="{y+20}" font-size="14">{counts[label]}</text>')
    destination.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><text x="15" y="28" font-size="18" font-weight="bold">Benchmark execution status</text>{"".join(bars)}</svg>\n',encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--sources",type=Path,required=True); parser.add_argument("--truth",type=Path,required=True); parser.add_argument("--runs",type=Path,required=True); parser.add_argument("--callers",type=Path,default=Path(__file__).resolve().parents[1]/"validation/manifests/caller_configurations.tsv"); parser.add_argument("--storage-audit",type=Path,required=True); parser.add_argument("--metrics-dir",type=Path); parser.add_argument("--output-dir",type=Path,required=True); args=parser.parse_args()
    output=args.output_dir; (output/"source_tables").mkdir(parents=True,exist_ok=True); (output/"figures").mkdir(exist_ok=True); (output/"tables").mkdir(exist_ok=True); (output/"report").mkdir(exist_ok=True)
    source_rows,truth_rows,run_rows,caller_rows=rows(args.sources),rows(args.truth),rows(args.runs),rows(args.callers); audit=json.loads(args.storage_audit.read_text(encoding="utf-8"))
    for source,name in [(args.sources,"datasets.tsv"),(args.truth,"truth_resources.tsv"),(args.runs,"benchmark_runs.tsv"),(args.callers,"caller_configurations.tsv")]: shutil.copyfile(source,output/"source_tables"/name)
    shutil.copyfile(args.callers,output/"tables"/"caller_configurations.tsv")
    metrics=[]
    if args.metrics_dir and args.metrics_dir.exists():
        for path in sorted(args.metrics_dir.glob("*.json")):
            item=json.loads(path.read_text(encoding="utf-8")); item["source_file"]=path.name; metrics.append(item)
    metric_rows=[]
    for item in metrics:
        for level in ("exact_all","exact_multi_exon","exact_single_exon","intron_chain","splice_junction"):
            if level in item: metric_rows.append({"dataset_id":item.get("dataset_id"),"platform":item.get("platform"),"caller":item.get("caller"),"level":level,**item[level]})
    if metric_rows: write_tsv(output/"tables"/"truth_metrics.tsv",metric_rows)
    counts=Counter(row["status"] for row in run_rows); status_svg(counts,output/"figures"/"benchmark_status.svg")
    write_tsv(output/"tables"/"experiment_status_summary.tsv",[{"status":status,"run_count":count} for status,count in sorted(counts.items())])
    dataset_counts=Counter((row["platform"],row["local_execution_status"]) for row in source_rows)
    write_tsv(output/"tables"/"dataset_summary.tsv",[{"platform":platform,"local_execution_status":status,"dataset_count":count} for (platform,status),count in sorted(dataset_counts.items())])
    overall="COMPLETE" if run_rows and all(row["status"]=="COMPLETE" for row in run_rows) else "PARTIAL_VALIDATION_DEFERRED_EXTERNAL_STORAGE"
    summary={"schema_version":"6.0","overall_status":overall,"production_defaults_changed":False,"source_dataset_records":len(source_rows),"benchmark_run_records":len(run_rows),"truth_metric_records":len(metric_rows),"storage_audit":audit,"run_status_counts":dict(sorted(counts.items())),"interpretation":"No global accuracy score is calculated. Deferred experiments are not treated as empty biological findings."}
    (output/"report"/"validation_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    metric_html=table(metric_rows,["dataset_id","platform","caller","level","tp","fp","fn","precision","recall","f1","status"]) if metric_rows else '<p class="status">No full external-data truth benchmark was executed locally; this is a deferred experiment, not a zero result.</p>'
    def experiment(name):
        selected=[row for row in run_rows if row["experiment"]==name]
        return table(selected,["run_id","dataset_id","caller","configuration","support_threshold","depth_fraction","truth_id","status"])
    document=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Block 6 validation report</title><style>body{{font:15px system-ui,sans-serif;max-width:1180px;margin:2rem auto;padding:0 1rem;color:#17242d}}h1,h2{{color:#174d6b}}.status{{padding:.8rem;background:#eef5f7;border-left:5px solid #356a8a}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd5da;padding:.4rem;text-align:left}}th{{background:#e8eff2}}.scroll{{overflow:auto}}code{{background:#eef1f2;padding:.1rem .25rem}}</style></head><body><h1>Block 6 scientific validation</h1><p class="status"><strong>Overall status:</strong> {overall}</p><p>This report keeps truth-bearing accuracy, caller concordance, platform/protocol effects, replicate reproducibility, depth, support thresholds, quantification, differential execution, simulation, and resources separate. It does not compute a global score.</p><h2>Feasibility gate</h2><p>Selected full WTC11 data: {audit['compressed_bytes']:,} compressed bytes. Estimated total local footprint: {audit['estimated_total_bytes']:,} bytes. Available at audit: {audit['available_bytes']:,} bytes. Decision: <code>{audit['decision']}</code>.</p><h2>Dataset manifest</h2>{table(source_rows,['dataset_id','experiment_accession','file_accession','platform','library_protocol','biological_replicate','benchmark_role','local_execution_status'])}<h2>Truth resources</h2>{table(truth_rows,['truth_id','resource_type','compressed_bytes','checksum_type','checksum','role'])}<h2>Caller configurations</h2>{table(caller_rows,['tool','version','role','configuration','support_semantics','container'])}<h2>All benchmark run statuses</h2><img src="../figures/benchmark_status.svg" alt="Benchmark status counts">{table(run_rows,['run_id','dataset_id','experiment','caller','configuration','support_threshold','depth_fraction','truth_id','status'])}<h2>Experimental SIRV truth</h2>{metric_html}{experiment('experimental_sirv_truth')}<h2>Simulation truth</h2>{experiment('simulation_truth')}<h2>Reference-annotation sensitivity</h2>{experiment('annotation_sensitivity')}<h2>Caller effect and concordance without truth</h2>{experiment('caller_effect')}<h2>Platform/protocol effect</h2>{experiment('platform_protocol_effect')}<h2>Replicate reproducibility</h2>{experiment('replicate_reproducibility')}<h2>Sequencing depth</h2>{experiment('depth')}<h2>Support thresholds</h2>{experiment('support_threshold')}<h2>Quantification validation</h2>{experiment('quantification')}<h2>Real differential-analysis execution</h2>{experiment('real_differential')}<h2>Runtime and resources</h2>{experiment('resources')}<p class="status">A status beginning with <code>PLANNED</code> or <code>DEFERRED</code> means unavailable locally, never zero effect.</p><h2>Scientific interpretation</h2><ul><li>Observed: reads, transcript structures, junctions, abundance, and resource measurements.</li><li>Evaluated against truth: only records tied to an explicit truth resource.</li><li>Concordance without truth: reproducibility or agreement, not correctness.</li><li>No production defaults were changed.</li></ul></body></html>'''
    (output/"report"/"index.html").write_text(document,encoding="utf-8")
    print(f"Wrote {output/'report'/'index.html'} ({overall})")


if __name__=="__main__": main()
