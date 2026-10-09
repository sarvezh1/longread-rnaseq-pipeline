#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "${script_dir}/../.." && pwd)
test_root=$(mktemp -d)
container_test=".block4_contract_tmp_$$"
trap 'rm -rf "${test_root}" "${project_root}/${container_test}"' EXIT
cd "${project_root}"

python3 -m py_compile bin/prepare_sample_abundance.py bin/build_abundance_matrices.py \
  bin/expression_qc.py bin/run_suppa.py bin/transcript_ends.py \
  bin/parse_transdecoder.py bin/integrate_biology.py bin/build_biology_contract.py

# Structural identity, not run-local IsoQuant ID, defines cross-sample rows.
python3 - "${test_root}" <<'PY'
import csv, json, pathlib, subprocess, sys
root=pathlib.Path(sys.argv[1]); fields=['transcript_id','gene_id','chrom','strand','transcript_start','transcript_end','transcript_length','exon_count','exon_class','exon_coordinates','intron_chain','structural_category','associated_gene','associated_transcript','supporting_read_count']
for sample,native,condition,replicate in [('s1','run1.1','control','1'),('s2','run2.9','treatment','1')]:
    bundle=root/sample; bundle.mkdir()
    (bundle/'sample_metadata.json').write_text(json.dumps({'sample':sample,'condition':condition,'replicate':replicate,'platform':'ont','library_type':'cdna'}))
    for name,feature,value,column in [('gene_counts','G1',10,'count'),('transcript_counts',native,10,'count'),('gene_tpm','G1',100,'TPM'),('transcript_tpm',native,100,'TPM')]:
        (bundle/f'{name}.tsv').write_text(f'feature_id\t{column}\n{feature}\t{value}\n')
    row=dict.fromkeys(fields,'NA'); row.update(transcript_id=native,gene_id='G1',chrom='chr1',strand='+',transcript_start='100',transcript_end='299',transcript_length='100',exon_count='2',exon_class='multi-exon',exon_coordinates='100-149,250-299',intron_chain='150-249',structural_category='novel_in_catalog',associated_gene='G1',supporting_read_count='10')
    with (bundle/'transcript_structure.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(row)
    (bundle/'transcripts.gtf').write_text(f'chr1\ttest\ttranscript\t100\t299\t.\t+\t.\tgene_id "G1"; transcript_id "{native}";\nchr1\ttest\texon\t100\t149\t.\t+\t.\tgene_id "G1"; transcript_id "{native}";\nchr1\ttest\texon\t250\t299\t.\t+\t.\tgene_id "G1"; transcript_id "{native}";\n')
    (bundle/'transcripts.fasta').write_text(f'>{native}\n'+'ATG'*33+'TAA\n')
    (bundle/'sqanti3_classification.tsv').write_text('isoform\tstructural_category\n')
    (bundle/'block2_contract.json').write_text('{}\n')
out=root/'matrices'
subprocess.run([sys.executable,'bin/build_abundance_matrices.py','--bundles',str(root/'s1'),str(root/'s2'),'--output',str(out)],check=True)
rows=list(csv.DictReader((out/'transcript_id_map.tsv').open(),delimiter='\t'))
assert len(rows)==2 and len({r['standard_transcript_id'] for r in rows})==1
assert {r['caller_transcript_id'] for r in rows}=={'run1.1','run2.9'}
assert json.loads((out/'design_status.json').read_text())['status']=='NOT_RUN_NO_EXPLICIT_CONTRAST'
PY

python3 bin/expression_qc.py --dataset tests/data/biology/dataset --output "${test_root}/qc"
python3 bin/transcript_ends.py --dataset tests/data/biology/dataset --output "${test_root}/ends" --window 50
docker run --rm -v "${project_root}:/work" -w /work longread-rnaseq-pipeline/biology:bioc3.23-r4.6 \
  Rscript bin/run_edger.R tests/data/biology/dataset "${container_test}/edger" treatment control 0.05 '~0 + condition'
docker run --rm -v "${project_root}:/work" -w /work longread-rnaseq-pipeline/biology:bioc3.23-r4.6 \
  Rscript bin/run_saturn.R tests/data/biology/dataset "${container_test}/dtu" treatment control 0.05 '~0 + condition'
docker run --rm -v "${project_root}:/work" -w /work longread-rnaseq-pipeline/biology:bioc3.23-r4.6 \
  Rscript bin/run_isoform_switch.R tests/data/biology/dataset "${container_test}/switch" treatment control 0.05 0.1
docker run --rm -v "${project_root}:/work" -w /work longread-rnaseq-pipeline/suppa:2.4 \
  python3 bin/run_suppa.py --dataset tests/data/biology/dataset --output "${container_test}/suppa" --numerator treatment --denominator control
docker run --rm -v "${project_root}:/work" -w /work trinityrnaseq/transdecoder:6.0.0 \
  sh -lc "mkdir -p '${container_test}/orfs/native' && TransDecoder -t tests/data/biology/dataset/production_transcripts.fasta -O '${container_test}/orfs/native' -m 100 >/dev/null && python3 bin/parse_transdecoder.py --metadata tests/data/biology/dataset/transcript_metadata.tsv --transdecoder-dir '${container_test}/orfs/native' --output '${container_test}/orfs'"
docker run --rm -v "${project_root}:/work" -w /work trinityrnaseq/transdecoder:6.0.0 \
  sh -lc "mkdir -p '${container_test}/orfs_sparse/native' && /usr/local/bin/util/TransDecoder.LongOrfs -t tests/data/biology/dataset/production_transcripts.fasta -O '${container_test}/orfs_sparse/native' -m 100 >/dev/null && python3 bin/parse_transdecoder.py --metadata tests/data/biology/dataset/transcript_metadata.tsv --transdecoder-dir '${container_test}/orfs_sparse/native' --prediction-status PREDICTION_NOT_AVAILABLE_SPARSE_INPUT --output '${container_test}/orfs_sparse'"

python3 - "${container_test}" "${test_root}" <<'PY'
import csv,json,pathlib
import sys
root=pathlib.Path(sys.argv[1]); test_root=pathlib.Path(sys.argv[2])
genes={r['gene_id']:r for r in csv.DictReader((root/'edger/genes.differential_expression.tsv').open(),delimiter='\t')}
assert float(genes['G1']['FDR'])<0.05 and float(genes['G3']['FDR'])>=0.05
dtu=list(csv.DictReader((root/'dtu/dtu_results.tsv').open(),delimiter='\t'))
assert any(r['transcript_id'] in {'T2a','T2b'} and float(r['adjusted_p_value'])<0.05 for r in dtu)
filtering={r['transcript_id']:r for r in csv.DictReader((root/'dtu/dtu_filtering.tsv').open(),delimiter='\t')}
assert filtering['T1']['reason']=='GENE_HAS_ONE_TRANSCRIPT' and filtering['TZERO']['reason']=='GENE_HAS_ONE_TRANSCRIPT'
events=list(csv.DictReader((root/'suppa/events.tsv').open(),delimiter='\t')); assert any(r['event_type']=='SE' for r in events)
diff=list(csv.DictReader((root/'suppa/differential_splicing.tsv').open(),delimiter='\t')); assert abs(float(diff[0]['delta_psi'])-0.8)<0.01
switch=list(csv.DictReader((root/'switch/isoform_switches.tsv').open(),delimiter='\t')); assert any(r['statistical_switch']=='TRUE' for r in switch)
orfs=list(csv.DictReader((root/'orfs/orf_predictions.tsv').open(),delimiter='\t'))
assert any(r['orf_status']=='ORF_CANDIDATE_AVAILABLE' for r in orfs) and any(r['transcript_id']=='TZERO' and r['orf_status']=='NO_ORF_CANDIDATE' for r in orfs)
sparse_status=json.loads((root/'orfs_sparse/orf_status.json').read_text())
sparse_orfs=list(csv.DictReader((root/'orfs_sparse/orf_predictions.tsv').open(),delimiter='\t'))
assert sparse_status['prediction_status']=='PREDICTION_NOT_AVAILABLE_SPARSE_INPUT'
assert sparse_status['predicted_orf_transcripts']==0 and sparse_status['candidate_orfs']>0
assert sparse_status['primary_orfs_fabricated'] is False
assert any(r['orf_status']=='ORF_CANDIDATE_AVAILABLE' for r in sparse_orfs)
assert not any(r['primary_orf'].lower()=='true' for r in sparse_orfs)
ends={r['transcript_id']:r for r in csv.DictReader((test_root/'ends/transcript_boundaries.tsv').open(),delimiter='\t')}
assert ends['T3']['putative_tss']=='2399' and ends['T3']['putative_tes']=='2000'
PY

nextflow run . -profile test_stub -stub-run -w "${test_root}/work" --outdir "${test_root}/results" \
  --mode biology --input tests/data/structural/ont_valid.csv
test -s "${test_root}/results/biology/biology_contract.json"
test -s "${test_root}/results/biology/abundance/transcript_counts.tsv"
test -s "${test_root}/results/biology/integrated/integrated_biology.tsv"

echo 'Block 4 structural and deterministic statistical contracts passed.'
