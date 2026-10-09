#!/usr/bin/env python3
"""Normalize TransDecoder candidate ORFs without equating absent predictions to non-coding."""
import argparse,csv,json,re
from pathlib import Path

def fasta(path):
    out={}; name=None; seq=[]
    if not path.exists(): return out
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name: out[name]=(header,"".join(seq))
            header=line[1:]; name=header.split()[0]; seq=[]
        else: seq.append(line.strip())
    if name: out[name]=(header,"".join(seq))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--metadata",type=Path,required=True);p.add_argument("--transdecoder-dir",type=Path,required=True);p.add_argument("--prediction-status",default="COMPLETE");p.add_argument("--output",type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True, exist_ok=True)
    with a.metadata.open(newline="") as h: models=list(csv.DictReader(h,delimiter="\t"))
    pep_files=list(a.transdecoder_dir.rglob("*.transdecoder.pep")); peps=fasta(pep_files[0]) if pep_files else {}; by_tx={}
    candidates=[]
    for oid,(header,sequence) in peps.items():
        tx=re.sub(r"\.p\d+$","",oid); length=len(sequence.rstrip("*")); complete=(re.search(r"type:([^\s]+)",header).group(1) if re.search(r"type:([^\s]+)",header) else "unknown")
        match=re.search(r"(\d+)-(\d+)\(([+-])\)",header); start,end,strand=(match.groups() if match else ("NA","NA","NA"))
        row={"transcript_id":tx,"orf_id":oid,"orf_status":"PREDICTED_ORF","orf_start":start,"orf_end":end,"orf_strand":strand,"frame":"NA","peptide_length":length,"completeness":complete,"primary_orf":False,"peptide_sequence":sequence}
        candidates.append(row); by_tx.setdefault(tx,[]).append(row)
    for tx,rows in by_tx.items(): max(rows,key=lambda x:x["peptide_length"])["primary_orf"]=True
    selected_ids={row["orf_id"] for row in candidates}
    candidate_files=list(a.transdecoder_dir.rglob("longest_orfs.pep")); raw_candidates=fasta(candidate_files[0]) if candidate_files else {}
    for oid,(header,sequence) in raw_candidates.items():
        if oid in selected_ids: continue
        tx=re.sub(r"\.p\d+$","",oid); length=len(sequence.rstrip("*")); complete=(re.search(r"type:([^\s]+)",header).group(1) if re.search(r"type:([^\s]+)",header) else "unknown")
        match=re.search(r"(\d+)-(\d+)\(([+-])\)",header); start,end,strand=(match.groups() if match else ("NA","NA","NA"))
        candidates.append({"transcript_id":tx,"orf_id":oid,"orf_status":"ORF_CANDIDATE_AVAILABLE","orf_start":start,"orf_end":end,"orf_strand":strand,"frame":"NA","peptide_length":length,"completeness":complete,"primary_orf":False,"peptide_sequence":sequence})
    candidate_transcripts={re.sub(r"\.p\d+$","",oid) for oid in raw_candidates}
    for model in models:
        if model["transcript_id"] not in by_tx and model["transcript_id"] not in candidate_transcripts: candidates.append({"transcript_id":model["transcript_id"],"orf_id":"NA","orf_status":"NO_ORF_CANDIDATE","orf_start":"NA","orf_end":"NA","orf_strand":"NA","frame":"NA","peptide_length":0,"completeness":"NA","primary_orf":False,"peptide_sequence":"NA"})
    fields=["transcript_id","orf_id","orf_status","orf_start","orf_end","orf_strand","frame","peptide_length","completeness","primary_orf","peptide_sequence"]
    with (a.output/"orf_predictions.tsv").open("w",newline="") as h:w=csv.DictWriter(h,fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(sorted(candidates,key=lambda x:(x["transcript_id"],not x["primary_orf"],x["orf_id"])))
    status="NO_TRANSCRIPTS" if not models else a.prediction_status
    (a.output/"orf_status.json").write_text(json.dumps({"status":status,"prediction_status":a.prediction_status,"transcripts":len(models),"predicted_orf_transcripts":len(by_tx),"candidate_orfs":len(raw_candidates),"candidate_transcripts":len(candidate_transcripts),"candidate_state":"ORF_CANDIDATE_AVAILABLE","no_candidate_state":"NO_ORF_CANDIDATE","primary_orfs_fabricated":False,"homology_search":"DISABLED","pfam_search":"DISABLED"},indent=2)+"\n")
if __name__=="__main__":main()
