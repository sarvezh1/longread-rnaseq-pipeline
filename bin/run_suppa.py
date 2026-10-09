#!/usr/bin/env python3
"""Run pinned SUPPA event/PSI/differential workflows with explicit statuses."""

import argparse, csv, json, subprocess
from pathlib import Path


def write_status(path, status, **extra):
    path.write_text(json.dumps({"status": status, "engine": "SUPPA", **extra}, indent=2) + "\n")


def empty(path, header): path.write_text("\t".join(header) + "\n")


def run(command, log):
    with log.open("a") as handle: subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def write_suppa_expression(source, destination, selected_samples=None):
    """Write SUPPA's unusual format: sample-only header, ID-prefixed rows."""
    table = read_table(source)
    samples = table[0][1:]
    selected = samples if selected_samples is None else selected_samples
    indexes = [samples.index(sample) + 1 for sample in selected]
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(selected)
        for row in table[1:]:
            writer.writerow([row[0], *[row[index] for index in indexes]])


def write_suppa_psi(source, destination, selected_samples):
    table = read_table(source)
    # Native SUPPA PSI has a sample-only header and event-ID-prefixed rows.
    samples = table[0]
    indexes = [samples.index(sample) + 1 for sample in selected_samples]
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(selected_samples)
        for row in table[1:]:
            writer.writerow([row[0], *[row[index] for index in indexes]])


def bh_adjust(values):
    adjusted = [None] * len(values)
    valid = [(index, value) for index, value in enumerate(values) if value is not None]
    previous = 1.0
    for rank_from_end, (index, value) in enumerate(sorted(valid, key=lambda item: item[1], reverse=True), 1):
        rank = len(valid) - rank_from_end + 1
        previous = min(previous, value * len(valid) / rank)
        adjusted[index] = previous
    return adjusted


def main():
    p=argparse.ArgumentParser(); p.add_argument("--dataset",type=Path,required=True); p.add_argument("--output",type=Path,required=True); p.add_argument("--suppa",default="suppa.py"); p.add_argument("--numerator",default=""); p.add_argument("--denominator",default="")
    a=p.parse_args(); a.output.mkdir(parents=True, exist_ok=True); log=a.output/"suppa.log"; log.write_text("")
    events=a.output/"events.tsv"; psi=a.output/"psi.tsv"; diff=a.output/"differential_splicing.tsv"
    empty(events,["event_id","gene_id","event_type","alternative_transcripts","total_transcripts"]); empty(psi,["event_id","event_type"]); empty(diff,["event_id","event_type","delta_psi","p_value","adjusted_p_value","contrast"])
    gtf=a.dataset/"production_transcripts.gtf"
    if not any(line and not line.startswith("#") for line in gtf.read_text().splitlines()):
        write_status(a.output/"events_status.json","COMPLETE_EMPTY_TRANSCRIPTOME"); write_status(a.output/"psi_status.json","NOT_RUN_NO_EVENTS"); write_status(a.output/"differential_status.json","NOT_RUN_NO_EVENTS"); return
    try:
        prefix=a.output/"production"
        run([a.suppa,"generateEvents","-i",str(gtf),"-o",str(prefix),"-f","ioe","-e","SE","SS","MX","RI","FL"],log)
        ioe_files=sorted(a.output.glob("production_*.ioe"))
        combined=a.output/"all_events.ioe"; rows=[]; header=None
        for source in ioe_files:
            lines=source.read_text().splitlines()
            if lines and header is None: header=lines[0]
            rows.extend(lines[1:]);
        combined.write_text((header+"\n" if header else "seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts\n")+"\n".join(rows)+("\n" if rows else ""))
        with events.open("w",newline="") as handle:
            w=csv.writer(handle,delimiter="\t",lineterminator="\n"); w.writerow(["event_id","gene_id","event_type","alternative_transcripts","total_transcripts"])
            for line in rows:
                f=line.split("\t"); event=f[2]; etype=event.split(";")[-1].split(":")[0] if ";" in event else "UNKNOWN"; w.writerow([event,f[1],etype,f[3],f[4]])
        write_status(a.output/"events_status.json","COMPLETE" if rows else "COMPLETE_EMPTY_RESULT",event_count=len(rows),native_event_types=["SE","A3","A5","MX","RI","AF","AL"])
        if not rows: write_status(a.output/"psi_status.json","NOT_RUN_NO_EVENTS"); write_status(a.output/"differential_status.json","NOT_RUN_NO_EVENTS"); return
        suppa_tpm=a.output/"transcript_tpm.suppa.tsv"
        write_suppa_expression(a.dataset/"transcript_tpm.tsv",suppa_tpm)
        run([a.suppa,"psiPerEvent","-i",str(combined),"-e",str(suppa_tpm),"-o",str(a.output/"production")],log)
        native=a.output/"production.psi"; native_psi=read_table(native)
        event_details={row[0]:row[1:] for row in read_table(events)[1:]}
        with psi.open("w",newline="") as handle:
            writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
            writer.writerow(["event_id","gene_id","event_type","alternative_transcripts","total_transcripts",*native_psi[0]])
            for row in native_psi[1:]: writer.writerow([row[0],*event_details[row[0]],*row[1:]])
        write_status(a.output/"psi_status.json","COMPLETE",missing_values_preserved_as_na=True,event_count=len(native_psi)-1)
        design=json.loads((a.dataset/"design_status.json").read_text())
        if design["status"]!="READY": write_status(a.output/"differential_status.json",design["status"],reasons=design["reasons"])
        else:
            metadata=list(csv.DictReader((a.dataset/"sample_metadata.tsv").open(newline=""),delimiter="\t"))
            numerator_samples=[row["sample"] for row in metadata if row["condition"]==a.numerator]
            denominator_samples=[row["sample"] for row in metadata if row["condition"]==a.denominator]
            inputs=[]
            for label,samples in (("denominator",denominator_samples),("numerator",numerator_samples)):
                condition_psi=a.output/f"{label}.psi"; condition_tpm=a.output/f"{label}.tpm"
                write_suppa_psi(native,condition_psi,samples)
                write_suppa_expression(a.dataset/"transcript_tpm.tsv",condition_tpm,samples)
                inputs.append((condition_psi,condition_tpm))
            differential_prefix=a.output/"differential"
            run([a.suppa,"diffSplice","--method","empirical","--psi",str(inputs[0][0]),str(inputs[1][0]),"--tpm",str(inputs[0][1]),str(inputs[1][1]),"--input",str(combined),"--area","1000","--lower-bound","0.05","--output",str(differential_prefix)],log)
            native_diff=a.output/"differential.dpsi"
            if not native_diff.exists():
                # SUPPA 2.4 leaves the single pairwise result at its temporary
                # name when there is only one comparison.
                candidates=sorted(a.output.glob("differential.dpsi.temp.*"))
                if len(candidates)!=1: raise RuntimeError("SUPPA did not produce an unambiguous differential result")
                native_diff=candidates[0]
            native_rows=read_table(native_diff)
            parsed=[]
            for row in native_rows[1:]:
                try: p_value=float(row[2]) if row[2] not in {"NA","nan",""} else None
                except (IndexError,ValueError): p_value=None
                try: delta=float(row[1]) if row[1] not in {"NA","nan",""} else None
                except (IndexError,ValueError): delta=None
                # Inputs are denominator then numerator; expose the requested
                # numerator-vs-denominator sign in the standardized contract.
                parsed.append((row[0],None if delta is None else -delta,p_value))
            adjusted=bh_adjust([row[2] for row in parsed])
            with diff.open("w",newline="") as handle:
                writer=csv.writer(handle,delimiter="\t",lineterminator="\n"); writer.writerow(["event_id","gene_id","event_type","alternative_transcripts","total_transcripts","delta_psi","p_value","adjusted_p_value","contrast"])
                for (event_id,delta,p_value),padj in zip(parsed,adjusted):
                    writer.writerow([event_id,*event_details[event_id],"NA" if delta is None else delta,"NA" if p_value is None else p_value,"NA" if padj is None else padj,f"{a.numerator}_vs_{a.denominator}"])
            status="COMPLETE" if parsed else "COMPLETE_EMPTY_RESULT"
            write_status(a.output/"differential_status.json",status,contrast=f"{a.numerator}_vs_{a.denominator}",event_count=len(parsed),method="SUPPA empirical",multiple_testing="Benjamini-Hochberg applied by workflow")
    except subprocess.CalledProcessError as error:
        write_status(a.output/"events_status.json","TOOL_ERROR",exit_code=error.returncode); write_status(a.output/"psi_status.json","NOT_RUN_UPSTREAM_ERROR"); write_status(a.output/"differential_status.json","NOT_RUN_UPSTREAM_ERROR"); raise

if __name__=="__main__": main()
