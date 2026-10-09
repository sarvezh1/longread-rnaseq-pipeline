#!/usr/bin/env python3
"""Fetch one manifest-selected benchmark file with preflight and checksum recording."""

import argparse, csv, hashlib, json, shutil, urllib.request
from urllib.parse import urlparse
from pathlib import Path


def digest(path,algorithm):
    value=hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): value.update(block)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--manifest",type=Path,required=True); parser.add_argument("--dataset-id",required=True); parser.add_argument("--output-dir",type=Path,required=True); parser.add_argument("--minimum-free-multiplier",type=float,default=2.0); args=parser.parse_args()
    with args.manifest.open(newline="",encoding="utf-8") as handle: matches=[row for row in csv.DictReader(handle,delimiter="\t") if row.get("dataset_id")==args.dataset_id or row.get("truth_id")==args.dataset_id]
    if len(matches)!=1: raise SystemExit(f"Expected one manifest record for {args.dataset_id}; found {len(matches)}")
    row=matches[0]; size=int(row.get("compressed_bytes","0")) if row.get("compressed_bytes","0").isdigit() else 0
    args.output_dir.mkdir(parents=True,exist_ok=True); free=shutil.disk_usage(args.output_dir).free
    if size and free < size*args.minimum_free_multiplier: raise SystemExit(f"Insufficient free space: {free} bytes available; require {int(size*args.minimum_free_multiplier)}")
    url=row["source_url"]; filename=Path(urlparse(url).path).name or args.dataset_id; destination=args.output_dir/filename; partial=destination.with_suffix(destination.suffix+".partial")
    with urllib.request.urlopen(url) as response, partial.open("wb") as output: shutil.copyfileobj(response,output)
    if size and partial.stat().st_size!=size: partial.unlink(); raise SystemExit("Downloaded size does not match manifest")
    checksum_type=row.get("checksum_type","").lower(); expected=row.get("checksum","").lower(); observed_sha256=digest(partial,"sha256")
    if checksum_type=="md5" and expected and digest(partial,"md5")!=expected: partial.unlink(); raise SystemExit("MD5 mismatch")
    partial.replace(destination)
    record={"dataset_id":args.dataset_id,"source_url":url,"bytes":destination.stat().st_size,"sha256":observed_sha256,"source_checksum_type":row.get("checksum_type"),"source_checksum":row.get("checksum"),"path":destination.name}
    (args.output_dir/f"{args.dataset_id}.fetch.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8"); print(json.dumps(record,indent=2))


if __name__=="__main__": main()
