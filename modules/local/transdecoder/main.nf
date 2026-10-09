process TRANSDECODER {
    tag 'TransDecoder ORF prediction'
    label 'process_high'
    container 'trinityrnaseq/transdecoder:6.0.0@sha256:78895668eaaa35d00659a98b30bd7341b512979547c0c517bb5ac49573975511'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'orfs'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.transdecoder.versions.yml'
    input: path dataset
    output: path 'orfs', emit: results; path 'block4.transdecoder.versions.yml', emit: versions
    script:
    """
    mkdir -p orfs/native
    prediction_status='NO_TRANSCRIPTS'
    if test -s ${dataset}/production_transcripts.fasta; then
        /usr/local/bin/util/TransDecoder.LongOrfs -t ${dataset}/production_transcripts.fasta -O orfs/native -m ${params.transdecoder_min_peptide_length} > orfs/native/transdecoder_longorfs.log 2>&1
        prediction_status='PREDICTION_COMPLETE'
        set +e
        /usr/local/bin/util/TransDecoder.Predict -t ${dataset}/production_transcripts.fasta -O orfs/native ${params.transdecoder_single_best_only ? '--single_best_only' : ''} > orfs/native/transdecoder_predict.log 2>&1
        predict_exit=\$?
        set -e
        if test \$predict_exit -ne 0; then
            candidate_count=\$(grep -c '^>' orfs/native/*transdecoder_dir/longest_orfs.pep 2>/dev/null || true)
            if test \$candidate_count -gt 0 && grep -Eq 'and pwm length = 0|no lines available in input' orfs/native/transdecoder_predict.log; then
                prediction_status='PREDICTION_NOT_AVAILABLE_SPARSE_INPUT'
            else
                cat orfs/native/transdecoder_predict.log >&2
                exit \$predict_exit
            fi
        fi
    fi
    python3 ${projectDir}/bin/parse_transdecoder.py --metadata ${dataset}/transcript_metadata.tsv --transdecoder-dir orfs/native --prediction-status "\$prediction_status" --output orfs
    printf '{"TransDecoder":"%s","container":"trinityrnaseq/transdecoder:6.0.0@sha256:78895668eaaa35d00659a98b30bd7341b512979547c0c517bb5ac49573975511"}\n' "\$(TransDecoder --version 2>&1 | awk 'NF{print \$NF}' | tail -1)" > orfs/runtime_versions.json
    printf 'TRANSDECODER:\n  transdecoder: "%s"\n  container: "trinityrnaseq/transdecoder:6.0.0@sha256:78895668eaaa35d00659a98b30bd7341b512979547c0c517bb5ac49573975511"\n  minimum_peptide_length: %s\n  homology: "disabled unless explicitly supplied in a future hook"\n  pfam: "disabled unless explicitly supplied in a future hook"\n' "\$(TransDecoder --version 2>&1 | awk 'NF{print \$NF}' | tail -1)" '${params.transdecoder_min_peptide_length}' > block4.transdecoder.versions.yml
    """
    stub:
    """
    mkdir -p orfs/native; printf 'transcript_id\torf_id\torf_status\torf_start\torf_end\torf_strand\tframe\tpeptide_length\tcompleteness\tprimary_orf\tpeptide_sequence\n' > orfs/orf_predictions.tsv; printf '{"status":"NO_TRANSCRIPTS","prediction_status":"NO_TRANSCRIPTS","homology_search":"DISABLED","pfam_search":"DISABLED"}' > orfs/orf_status.json; printf '{"TransDecoder":"6.0.0"}' > orfs/runtime_versions.json; printf 'TRANSDECODER:\n  transdecoder: "6.0.0"\n' > block4.transdecoder.versions.yml
    """
}
