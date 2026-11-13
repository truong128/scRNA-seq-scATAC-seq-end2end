import csv
import gzip
import io
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from prep_fastq import unpack
from download import get
from verify_local_inputs import check

class TestInputHelpers(unittest.TestCase):
    def test_gse_manifest(self):
        with (ROOT/'data'/'geo_pilot.csv').open() as f:
            rows=list(csv.DictReader(f))
        self.assertEqual({r['sample_id'] for r in rows},{'AL05','AL06'})
        for r in rows:
            self.assertIn(r['gsm'],r['filtered_url'])
            self.assertIn(r['gsm'],r['raw_url'])
            self.assertIn('TCR_filtered_contig_annotations.csv.gz',r['tcr_url'])
            self.assertTrue(r['tcr_url'].startswith('https://ftp.ncbi.nlm.nih.gov/geo/samples/'))

    def test_untrusted_tar_path_rejected(self):
        import tarfile
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)
            name=path/'bad.tar'
            with tarfile.open(name,'w') as f:
                payload=b'bad'
                m=tarfile.TarInfo('../../escape.txt')
                m.size=len(payload)
                f.addfile(m,io.BytesIO(payload))
            with self.assertRaises(ValueError):
                unpack(name,path/'out')
            self.assertFalse((path/'escape.txt').exists())

    def test_processed_inputs_are_not_fastqs(self):
        with (ROOT/'data'/'geo_pilot.csv').open() as f:
            rows=list(csv.DictReader(f))
        self.assertTrue(all('.h5' in r['filtered_url'] for r in rows))


class TestOfflineInputCheck(unittest.TestCase):
    def _write_package(self, root, files, samples=("AL05",)):
        inputs = root / "inputs"
        inputs.mkdir()
        rows = []
        for name, payload in files.items():
            (inputs / name).write_bytes(payload)
            rows.append(f"{name}\t{len(payload)}\t{__import__('hashlib').sha256(payload).hexdigest()}")
        geo = [f"{line.split(chr(9))[0]}\tGEO_{line.split(chr(9))[0]}\t{line.split(chr(9))[1]}\t{line.split(chr(9))[2]}" for line in rows]
        (inputs / "GEO_SOURCE.tsv").write_text("short_name\tgeo_filename\tbytes\tsha256\n" + "\n".join(geo) + "\n")
        sheet = root / "samples.csv"
        sheet.write_text("sample_id,patient_id,condition\n" + "".join(f"{s},{s},WT_LUAD\n" for s in samples))
        return sheet

    def test_package_without_geo_supplement_passes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sheet = self._write_package(root, {
                "AL05_filtered.h5": b"filtered",
                "AL05_raw.h5": b"raw",
                "AL05_tcr.csv.gz": b"tcr",
            })
            message = check(root, sheet)
            self.assertIn("Workflow input checksums OK (3 files)", message)
            self.assertIn("GSE274934_suppl/manifest.tsv absent", message)
            self.assertFalse((root / "GSE274934_suppl").exists())

    def test_missing_input_manifest_names_the_required_file(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "inputs").mkdir()
            sheet = root / "samples.csv"
            sheet.write_text("sample_id,patient_id,condition\nAL05,AL05,WT_LUAD\n")
            with self.assertRaises(ValueError) as caught:
                check(root, sheet)
            self.assertIn("GEO_SOURCE.tsv", str(caught.exception))
            self.assertIn("not required", str(caught.exception))

    def test_present_geo_supplement_is_size_checked(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sheet = self._write_package(root, {
                "AL05_filtered.h5": b"filtered",
                "AL05_raw.h5": b"raw",
                "AL05_tcr.csv.gz": b"tcr",
            })
            suppl = root / "GSE274934_suppl"
            suppl.mkdir()
            (suppl / "extra.bin").write_bytes(b"abc")
            (suppl / "manifest.tsv").write_text("extra.bin\t3\n")
            self.assertIn("GEO supplement size check OK (1 files)", check(root, sheet))
            (suppl / "manifest.tsv").write_text("extra.bin\t9\n")
            with self.assertRaises(ValueError) as caught:
                check(root, sheet)
            self.assertIn("GEO size check failed", str(caught.exception))

    def test_half_an_atac_pair_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sheet = self._write_package(root, {
                "AL05_filtered.h5": b"filtered",
                "AL05_raw.h5": b"raw",
                "AL05_tcr.csv.gz": b"tcr",
                "AL05_atac_peaks.h5": b"peaks",
            })
            with self.assertRaises(ValueError) as caught:
                check(root, sheet)
            self.assertIn("only one of the two ATAC files", str(caught.exception))

if __name__=='__main__': unittest.main()
