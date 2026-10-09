#!/usr/bin/env python3
"""Create or validate a complete Block 6 benchmark execution record."""

import argparse, hashlib, json
from pathlib import Path

REQUIRED = ["run_id","dataset_id","source_url","source_checksum","subset_method","random_seed","platform","library_protocol","biological_replicate","reference","annotation","caller","caller_version","parameters","support_threshold","read_depth","container","container_digest","runtime","metric_definitions"]


def sha256(path):
    value=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): value.update(block)
    return value.hexdigest()


def validate(record):
    missing=[field for field in REQUIRED if field not in record or record[field] in (None,"")]
    if missing: raise ValueError("benchmark provenance lacks: "+", ".join(missing))
    if record.get("production_defaults_changed") is not False: raise ValueError("benchmark record must state production_defaults_changed=false")


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--validate",type=Path); parser.add_argument("--output",type=Path); parser.add_argument("--run-id"); parser.add_argument("--dataset-id"); parser.add_argument("--source-url"); parser.add_argument("--source-checksum"); parser.add_argument("--subset-method",default="none; full source population"); parser.add_argument("--random-seed",default="not applicable"); parser.add_argument("--platform"); parser.add_argument("--library-protocol"); parser.add_argument("--biological-replicate"); parser.add_argument("--reference",type=Path); parser.add_argument("--annotation",type=Path); parser.add_argument("--caller"); parser.add_argument("--caller-version"); parser.add_argument("--parameters",default="production defaults"); parser.add_argument("--support-threshold",default="production default"); parser.add_argument("--read-depth",default="100%"); parser.add_argument("--container"); parser.add_argument("--container-digest"); parser.add_argument("--runtime"); parser.add_argument("--metric-definitions",type=Path); args=parser.parse_args()
    if args.validate:
        record=json.loads(args.validate.read_text(encoding="utf-8")); validate(record); print(f"Valid benchmark provenance: {args.validate}"); return
    if not args.output: raise SystemExit("--output is required when creating a record")
    record={"schema_version":"6.0","run_id":args.run_id,"dataset_id":args.dataset_id,"source_url":args.source_url,"source_checksum":args.source_checksum,"subset_method":args.subset_method,"random_seed":args.random_seed,"platform":args.platform,"library_protocol":args.library_protocol,"biological_replicate":args.biological_replicate,"reference":str(args.reference),"reference_sha256":sha256(args.reference),"annotation":str(args.annotation),"annotation_sha256":sha256(args.annotation),"caller":args.caller,"caller_version":args.caller_version,"parameters":args.parameters,"support_threshold":args.support_threshold,"read_depth":args.read_depth,"container":args.container,"container_digest":args.container_digest,"runtime":args.runtime,"metric_definitions":str(args.metric_definitions),"metric_definitions_sha256":sha256(args.metric_definitions),"production_defaults_changed":False}
    validate(record); args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__": main()
