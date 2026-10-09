#!/usr/bin/env python3
"""Audit Block 6 source-data storage before any download."""

import argparse, csv, json, shutil
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--manifest",type=Path,required=True); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--target",type=Path,default=Path(".")); parser.add_argument("--unpacked-multiplier",type=float,default=3.0); parser.add_argument("--work-multiplier",type=float,default=6.0); parser.add_argument("--result-multiplier",type=float,default=1.0); args=parser.parse_args()
    with args.manifest.open(newline="",encoding="utf-8") as handle: rows=list(csv.DictReader(handle,delimiter="\t"))
    known=[int(row["compressed_bytes"]) for row in rows if row["compressed_bytes"].isdigit() and row["local_execution_status"].startswith("FULL_DEFERRED")]
    compressed=sum(known); unpacked=int(compressed*args.unpacked_multiplier); work=int(compressed*args.work_multiplier); results=int(compressed*args.result_multiplier); required=compressed+unpacked+work+results
    free=shutil.disk_usage(args.target.resolve()).free
    report={"schema_version":"6.0","manifest":args.manifest.as_posix(),"datasets_with_exact_size":len(known),"compressed_bytes":compressed,"estimated_unpacked_bytes":unpacked,"estimated_work_bytes":work,"estimated_result_bytes":results,"estimated_total_bytes":required,"available_bytes":free,"fits":required<=free,"decision":"FULL_LOCAL_EXECUTION_ALLOWED" if required<=free else "FULL_LOCAL_EXECUTION_DEFERRED_STORAGE","estimation_policy":{"unpacked_multiplier":args.unpacked_multiplier,"work_multiplier":args.work_multiplier,"result_multiplier":args.result_multiplier}}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8"); print(json.dumps(report,indent=2))


if __name__=="__main__": main()
