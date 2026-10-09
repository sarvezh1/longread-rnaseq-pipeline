#!/usr/bin/env python3
"""Fast deterministic Block 6 metric, manifest, and reporting tests."""

import csv, json, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PY=sys.executable


def run(*args):
    return subprocess.run([PY,*map(str,args)],cwd=ROOT,check=True,text=True,capture_output=True)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.context=tempfile.TemporaryDirectory(); self.temp=Path(self.context.name)
        self.data=ROOT/"tests/data/validation"
    def tearDown(self): self.context.cleanup()

    def truth(self,prediction="prediction.gtf",truth="truth.gtf",tolerance=0,name="metrics"):
        output=self.temp/f"{name}.json"; matches=self.temp/f"{name}.matches.tsv"
        run("bin/validation_metrics.py","truth-match","--prediction",self.data/prediction,"--truth",self.data/truth,"--end-tolerance",tolerance,"--dataset-id","deterministic_fixture","--caller","fixture","--platform","synthetic","--truth-type","deterministic test fixture","--output",output,"--matches",matches)
        return json.loads(output.read_text()),matches

    def test_truth_metrics_and_tolerance(self):
        zero,_=self.truth(); ten,matches=self.truth(tolerance=10,name="tol10")
        self.assertEqual(zero["exact_all"]|{"ignore":0},{"tp":2,"fp":2,"fn":1,"precision":0.5,"recall":2/3,"f1":0.5714285714285715,"status":"COMPLETE","ignore":0})
        self.assertEqual(zero["intron_chain"]["tp"],2); self.assertEqual(zero["splice_junction"]["fp"],1)
        self.assertEqual(ten["exact_all"]["tp"],3); self.assertEqual(ten["exact_all"]["fn"],0)
        with matches.open(newline="") as handle: rows=list(csv.DictReader(handle,delimiter="\t"))
        self.assertEqual(len(rows),2); self.assertIn("absolute_tss_error_bp",rows[0])

    def test_empty_truth_and_zero_detection(self):
        empty_truth,_=self.truth(truth="empty.gtf",name="empty_truth")
        zero_detection,_=self.truth(prediction="empty.gtf",name="zero_detection")
        self.assertEqual(empty_truth["exact_all"]["status"],"EMPTY_TRUTH"); self.assertIsNone(empty_truth["exact_all"]["precision"])
        self.assertEqual(zero_detection["exact_all"]["status"],"ZERO_DETECTION"); self.assertEqual(zero_detection["exact_all"]["recall"],0.0)

    def test_downsampling_is_deterministic_and_nested(self):
        paths={}
        for fraction in (0.25,0.75):
            output=self.temp/f"d{fraction}.fastq"; manifest=self.temp/f"d{fraction}.json"
            run("bin/validation_metrics.py","downsample","--input",self.data/"reads.fastq","--output",output,"--fraction",fraction,"--seed",61703,"--manifest",manifest); paths[fraction]=output
            self.assertEqual(json.loads(manifest.read_text())["input_reads"],8)
        low=set(paths[0.25].read_text().splitlines()[::4]); high=set(paths[0.75].read_text().splitlines()[::4]); self.assertTrue(low<=high)
        repeat=self.temp/"repeat.fastq"; run("bin/validation_metrics.py","downsample","--input",self.data/"reads.fastq","--output",repeat,"--fraction",0.75,"--seed",61703,"--manifest",self.temp/"repeat.json")
        self.assertEqual(repeat.read_bytes(),paths[0.75].read_bytes())
        run("bin/validation_metrics.py","fastq-stats","--input",self.data/"reads.fastq","--output",self.temp/"stats.json")
        stats=json.loads((self.temp/"stats.json").read_text()); self.assertEqual((stats["reads"],stats["total_bases"],stats["n50"]),(8,32,4))

    def test_depth_manifest_generation(self):
        output=self.temp/"depth.tsv"
        run("bin/validation_metrics.py","depth-manifest","--source-manifest",ROOT/"validation/manifests/source_datasets.tsv","--dataset-id","wtc11_ont_cdna_rep1","--fractions",0.1,0.5,1.0,"--seed",61703,"--output",output)
        with output.open(newline="") as handle: records=list(csv.DictReader(handle,delimiter="\t"))
        self.assertEqual([row["fraction"] for row in records],["0.10","0.50","1.00"]); self.assertTrue(all(row["biological_replicate_claim"]=="false" for row in records))

    def test_support_replicates_and_quantification(self):
        run("bin/validation_metrics.py","support-threshold","--input",self.data/"structures_r1.tsv","--support-semantics","fixture read count","--thresholds",1,2,5,10,"--output-dir",self.temp/"thresholds","--summary",self.temp/"thresholds.tsv")
        with (self.temp/"thresholds.tsv").open(newline="") as handle: thresholds=list(csv.DictReader(handle,delimiter="\t"))
        self.assertEqual([row["retained_models"] for row in thresholds],["2","1","1","0"])
        run("bin/validation_metrics.py","support-benchmark","--input",self.data/"structures_r1.tsv","--truth",self.data/"truth.gtf","--support-semantics","fixture read count","--thresholds",1,5,10,"--output",self.temp/"support_metrics.tsv")
        with (self.temp/"support_metrics.tsv").open(newline="") as handle: support=list(csv.DictReader(handle,delimiter="\t"))
        exact=[row for row in support if row["level"]=="exact_all"]; self.assertEqual([row["tp"] for row in exact],["2","1","0"])
        run("bin/validation_metrics.py","replicates","--structure",f"r1={self.data/'structures_r1.tsv'}","--structure",f"r2={self.data/'structures_r2.tsv'}","--output",self.temp/"replicates.tsv","--membership",self.temp/"membership.json")
        membership=json.loads((self.temp/"membership.json").read_text()); self.assertEqual(membership["exact_structures_by_replicate_support"],{"1":2,"2":1})
        run("bin/validation_metrics.py","quantification","--truth",self.data/"quant_truth.tsv","--predicted",self.data/"quant_predicted.tsv","--output",self.temp/"quant.json","--table",self.temp/"quant.tsv")
        quant=json.loads((self.temp/"quant.json").read_text()); self.assertEqual(quant["n_union"],5); self.assertEqual(quant["zero_policy"],"union of truth and prediction; missing values represented as zero")

    def test_provenance_and_publication_report(self):
        record=self.temp/"run.json"
        run("bin/benchmark_provenance.py","--output",record,"--run-id","fixture","--dataset-id","fixture","--source-url","https://example.org/fixture","--source-checksum","sha256:fixture","--platform","synthetic","--library-protocol","fixture","--biological-replicate","1","--reference",self.data/"truth.gtf","--annotation",self.data/"truth.gtf","--caller","fixture","--caller-version","1.0","--container","fixture:1.0","--container-digest","sha256:fixture","--runtime","1 second","--metric-definitions",ROOT/"validation/manifests/metric_definitions.json")
        run("bin/benchmark_provenance.py","--validate",record)
        report=self.temp/"report"
        run("bin/build_validation_report.py","--sources",ROOT/"validation/manifests/source_datasets.tsv","--truth",ROOT/"validation/manifests/truth_resources.tsv","--runs",ROOT/"validation/manifests/benchmark_runs.tsv","--storage-audit",ROOT/"validation/resources/storage_audit.json","--output-dir",report)
        self.assertTrue((report/"report/index.html").exists()); self.assertTrue((report/"figures/benchmark_status.svg").exists()); self.assertIn("PARTIAL_VALIDATION_DEFERRED_EXTERNAL_STORAGE",(report/"report/index.html").read_text())
        run("bin/summarize_validation_resources.py","--trace",self.data/"trace.tsv","--output",self.temp/"resources.tsv","--summary",self.temp/"resources.json")
        resources=json.loads((self.temp/"resources.json").read_text()); self.assertEqual(resources["complete_or_cached"],2); self.assertEqual(resources["maximum_peak_rss_bytes"],2000000000)
        figure=self.temp/"depth.svg"; run("bin/plot_validation.py","--input",self.data/"plot.tsv","--x","depth","--y","f1","--group","caller","--title","Deterministic fixture","--x-label","Depth fraction","--y-label","F1","--output",figure)
        self.assertIn("IsoQuant",figure.read_text()); self.assertIn("Depth fraction",figure.read_text())

    def test_authoritative_manifest_pairing(self):
        with (ROOT/"validation/manifests/source_datasets.tsv").open(newline="",encoding="utf-8") as handle: manifest=list(csv.DictReader(handle,delimiter="\t"))
        wtc=[row for row in manifest if row["source_material"]=="WTC11"]
        self.assertEqual(len(wtc),6); self.assertEqual({row["biological_replicate"] for row in wtc},{"1","2","3"})
        for replicate in ("1","2","3"):
            self.assertEqual({row["platform"] for row in wtc if row["biological_replicate"]==replicate},{"ONT","PacBio"})
        with (ROOT/"validation/manifests/benchmark_runs.tsv").open(newline="",encoding="utf-8") as handle: runs={row["run_id"]:row for row in csv.DictReader(handle,delimiter="\t")}
        expected={row["dataset_id"] for row in wtc}
        for run_id in ("platform_protocol_triplicates","replicate_reproducibility"):
            self.assertEqual(set(runs[run_id]["dataset_id"].split(";")),expected)
        with (ROOT/"validation/manifests/differential_candidate.tsv").open(newline="",encoding="utf-8") as handle: differential={row["sample"] for row in csv.DictReader(handle,delimiter="\t")}
        self.assertEqual(set(runs["real_differential"]["dataset_id"].split(";")),differential)


if __name__=="__main__": unittest.main(verbosity=2)
