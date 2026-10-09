#!/usr/bin/env python3
"""Join explicit Block 4 evidence without producing narrative biological claims."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def rows(path,key):
    if not path.exists(): return {}
    with path.open(newline="") as h:return {r.get(key):r for r in csv.DictReader(h,delimiter="\t") if r.get(key)}
def matrix(path):
    with path.open(newline="") as h:
        rd=csv.DictReader(h,delimiter="\t"); return rd.fieldnames[1:],{r["feature_id"]:r for r in rd}
def primary_orfs(path):
    selected={}; absent={}; candidates={}
    if not path.exists(): return selected
    with path.open(newline="") as h:
        for row in csv.DictReader(h,delimiter="\t"):
            tx=row.get("transcript_id")
            if row.get("primary_orf","").lower()=="true": selected[tx]=row
            elif row.get("orf_status")=="NO_ORF_CANDIDATE": absent[tx]=row
            elif row.get("orf_status")=="ORF_CANDIDATE_AVAILABLE" and (tx not in candidates or int(row.get("peptide_length",0))>int(candidates[tx].get("peptide_length",0))): candidates[tx]=row
    result=dict(absent)
    for tx,row in candidates.items(): result[tx]={**row,"orf_id":"NA","primary_orf":"False"}
    result.update(selected); return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--dataset",type=Path,required=True);p.add_argument("--statistics",type=Path,required=True);p.add_argument("--dtu",type=Path,required=True);p.add_argument("--splicing",type=Path,required=True);p.add_argument("--ends",type=Path,required=True);p.add_argument("--orfs",type=Path,required=True);p.add_argument("--switches",type=Path,required=True);p.add_argument("--enrichment",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True, exist_ok=True)
    samples,counts=matrix(a.dataset/"transcript_counts.tsv");_,tpm=matrix(a.dataset/"transcript_tpm.tsv")
    meta=rows(a.dataset/"transcript_metadata.tsv","transcript_id"); dte=rows(a.statistics/"transcripts.differential_expression.tsv","transcript_id");dtu=rows(a.dtu/"dtu_results.tsv","transcript_id");ends=rows(a.ends/"transcript_boundaries.tsv","transcript_id")
    orf_rows=primary_orfs(a.orfs/"orf_predictions.tsv"); switches=rows(a.switches/"isoform_switches.tsv","transcript_id")
    events=defaultdict(list)
    with (a.splicing/"events"/"events.tsv").open(newline="") as h:
        for r in csv.DictReader(h,delimiter="\t"):
            for tx in (r.get("total_transcripts","")+","+r.get("alternative_transcripts","")).split(","):
                if tx:events[tx].append(r.get("event_id",""))
    fields=["gene_id","transcript_id","sample","sqanti3_class","known_novel","raw_count","tpm","dte_logFC","dte_p_value","dte_fdr","dtu_estimate","dtu_p_value","dtu_fdr","splicing_event_ids","putative_tss","putative_tes","tss_cluster","tes_cluster","orf_status","primary_orf_id","peptide_length","isoform_switch_dIF","isoform_switch_fdr","statistical_switch","functional_consequence_evidence","caller_benchmark_context","observed_evidence","statistical_evidence","computational_prediction","experimental_validation"]
    with (a.output/"integrated_biology.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for tx in sorted(meta):
            m=meta[tx];o=orf_rows.get(tx,{});sw=switches.get(tx,{});e=ends.get(tx,{})
            for sample in samples:w.writerow({"gene_id":m.get("gene_id","NA"),"transcript_id":tx,"sample":sample,"sqanti3_class":m.get("structural_category","NA"),"known_novel":m.get("known_novel","NA"),"raw_count":counts.get(tx,{}).get(sample,"0"),"tpm":tpm.get(tx,{}).get(sample,"0"),"dte_logFC":dte.get(tx,{}).get("logFC","NA"),"dte_p_value":dte.get(tx,{}).get("PValue","NA"),"dte_fdr":dte.get(tx,{}).get("FDR","NA"),"dtu_estimate":dtu.get(tx,{}).get("estimate","NA"),"dtu_p_value":dtu.get(tx,{}).get("p_value","NA"),"dtu_fdr":dtu.get(tx,{}).get("adjusted_p_value","NA"),"splicing_event_ids":";".join(sorted(set(events[tx]))) or "NA","putative_tss":e.get("putative_tss","NA"),"putative_tes":e.get("putative_tes","NA"),"tss_cluster":e.get("tss_cluster","NA"),"tes_cluster":e.get("tes_cluster","NA"),"orf_status":o.get("orf_status","NO_ORF_CANDIDATE"),"primary_orf_id":o.get("orf_id","NA"),"peptide_length":o.get("peptide_length","0"),"isoform_switch_dIF":sw.get("dIF","NA"),"isoform_switch_fdr":sw.get("adjusted_p_value","NA"),"statistical_switch":sw.get("statistical_switch","NA"),"functional_consequence_evidence":sw.get("functional_consequence_evidence","NOT_EVALUATED"),"caller_benchmark_context":"NOT_ATTACHED","observed_evidence":"production IsoQuant structure and abundance","statistical_evidence":"see method-specific status/results","computational_prediction":"ORF and enrichment outputs when available","experimental_validation":"NONE_PROVIDED"})
    (a.output/"integrated_status.json").write_text(json.dumps({"status":"COMPLETE","natural_language_claims_generated":False,"caller_benchmark_context_policy":"optional context only; never filters or changes significance","evidence_levels":["observed","statistically_inferred","computationally_predicted","experimentally_validated"]},indent=2)+"\n")
if __name__=="__main__":main()
