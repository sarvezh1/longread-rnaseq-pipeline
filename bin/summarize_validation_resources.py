#!/usr/bin/env python3
"""Convert a Nextflow trace into stable Block 6 resource source tables."""

import argparse, csv, json, re
from pathlib import Path

UNITS={"B":1,"KB":1000,"MB":1000**2,"GB":1000**3,"TB":1000**4,"KIB":1024,"MIB":1024**2,"GIB":1024**3,"TIB":1024**4}


def bytes_value(text):
    if not text or text=="-": return None
    match=re.fullmatch(r"([0-9.]+)\s*([KMGT]?i?B)",text.strip(),re.I)
    return int(float(match.group(1))*UNITS[match.group(2).upper()]) if match else None


def seconds(text):
    if not text or text=="-": return None
    total=0.0
    for value,unit in re.findall(r"([0-9.]+)\s*(ms|s|m|h|d)",text,re.I): total+=float(value)*{"ms":.001,"s":1,"m":60,"h":3600,"d":86400}[unit.lower()]
    return total


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--trace",type=Path,required=True); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--summary",type=Path,required=True); args=parser.parse_args()
    with args.trace.open(newline="",encoding="utf-8") as handle: source=list(csv.DictReader(handle,delimiter="\t"))
    output=[]
    for row in source:
        output.append({"task_id":row.get("task_id",row.get("task_id","NA")),"process":row.get("process","NA"),"tag":row.get("tag","NA"),"status":row.get("status","NA"),"wall_seconds":seconds(row.get("realtime","")),"cpu_percent":row.get("%cpu","NA"),"peak_rss_bytes":bytes_value(row.get("peak_rss","")),"read_bytes":bytes_value(row.get("read_bytes","")),"write_bytes":bytes_value(row.get("write_bytes",""))})
    fields=["task_id","process","tag","status","wall_seconds","cpu_percent","peak_rss_bytes","read_bytes","write_bytes"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",newline="",encoding="utf-8") as handle: writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(output)
    complete=[row for row in output if row["status"] in {"COMPLETED","CACHED"}]
    summary={"schema_version":"6.0","tasks":len(output),"complete_or_cached":len(complete),"failed":sum(row["status"]=="FAILED" for row in output),"wall_seconds_sum":sum(row["wall_seconds"] or 0 for row in complete),"maximum_peak_rss_bytes":max((row["peak_rss_bytes"] or 0 for row in complete),default=None),"note":"summed task wall time is not end-to-end elapsed workflow time"}
    args.summary.write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__": main()
