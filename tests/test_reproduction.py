"""Tests for the new reproduction layer and a small scientific regression set."""
from __future__ import annotations
import importlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import reproduce
from tools.compare_reference import compare_json, compare_csv
from tools import build_release


class RunnerTests(unittest.TestCase):
    def test_source_integrity(self):
        result=reproduce.check_manifest(ROOT/'SOURCE_MANIFEST.sha256')
        self.assertTrue(result['passed'],result['failures'])
        self.assertEqual(result['checked_files'],73)

    def test_paper_metadata(self):
        data=json.loads((ROOT/'SOURCE_PROVENANCE.json').read_text())
        self.assertEqual(data['paper'], {'title': 'A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design', 'author': 'G. Calafiore'})
        self.assertEqual(data['package_version'], build_release.VERSION)
        self.assertEqual(len(data['files']),73)
        for item in data['files']:
            self.assertEqual(set(item),{'path','sha256'})
            self.assertEqual(reproduce.sha256(ROOT/item['path']),item['sha256'])

    def test_missing_manifest_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=reproduce.check_manifest(Path(tmp)/'SOURCE_MANIFEST.sha256')
            self.assertFalse(result['passed'])
            self.assertEqual(result['checked_files'],0)
            self.assertIn('missing required manifest',result['failures'][0]['reason'])

    def test_release_requires_core_files(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(build_release,'ROOT',Path(tmp)):
            with self.assertRaisesRegex(FileNotFoundError,'SOURCE_MANIFEST'):
                build_release.release_files()

    def test_release_contains_required_manifests(self):
        names={p.relative_to(ROOT).as_posix() for p in build_release.release_files()}
        self.assertIn('SOURCE_MANIFEST.sha256',names)
        self.assertIn('SHA256SUMS',names)
        self.assertIn('SOURCE_PROVENANCE.json',names)
        self.assertIn('.github/workflows/ci.yml',names)
        self.assertFalse(any(n.startswith('runs/') for n in names))

    def test_metadata_uses_repository_relative_paths(self):
        data=json.loads((ROOT/'SOURCE_PROVENANCE.json').read_text())
        for item in data['files']:
            path=Path(item['path'])
            self.assertFalse(path.is_absolute())
            self.assertNotIn('..',path.parts)
            self.assertTrue((ROOT/path).is_file())

    def test_full_starts_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            workspace=Path(tmp)/'run'
            record=reproduce.prepare_workspace(workspace,False)
            self.assertEqual(record['seeded_reference_inputs'],[])
            self.assertTrue(record['starts_from_empty_outputs'])
            self.assertFalse(any((workspace/'data').iterdir()))
            self.assertFalse(any((workspace/'checks').iterdir()))
            self.assertEqual(len(list((workspace/'code').glob('*.py'))),16)

    def test_seeded_inputs_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            workspace=Path(tmp)/'run'
            record=reproduce.prepare_workspace(workspace,True)
            self.assertFalse(record['starts_from_empty_outputs'])
            self.assertGreater(len(record['seeded_reference_inputs']),20)
            self.assertTrue((workspace/'data/mechanical_results.json').is_file())
            self.assertFalse((workspace/'checks/reported_results_verification.json').exists())
            self.assertFalse(any((workspace/'figures').iterdir()))

    def test_existing_outputs_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            workspace=Path(tmp)/'run';workspace.mkdir()
            saved=workspace/'keep.txt';saved.write_text('keep this')
            with self.assertRaises(ValueError):reproduce.prepare_workspace(workspace,False)
            self.assertEqual(saved.read_text(),'keep this')

    def test_protected_paths_rejected(self):
        for path in [ROOT,ROOT.parent,ROOT/'reference/new',ROOT/'code/new',ROOT/'validation/new']:
            with self.subTest(path=path):
                with self.assertRaises(ValueError):reproduce.prepare_workspace(path,False)

    def test_child_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            env=reproduce.child_environment(Path(tmp))
            self.assertEqual(env['MPLBACKEND'],'Agg')
            self.assertNotIn('PYTHONOPTIMIZE',env)
            self.assertNotIn('PYTHONPATH',env)
            self.assertTrue(all(env[k]=='1' for k in reproduce.THREAD_VARIABLES))

    def test_five_source_verifiers(self):
        self.assertEqual(len(reproduce.VERIFY),5)
        self.assertEqual(reproduce.FULL[-1],'verify_reported_results')


class ComparisonTests(unittest.TestCase):
    def test_elapsed_time_is_excluded(self):
        self.assertEqual(compare_json({'value':1,'elapsed_seconds':2.},{'value':1,'elapsed_seconds':999.},1e-9,1e-11),[])

    def test_integer_counts_are_exact(self):
        self.assertTrue(compare_json({'failures':34},{'failures':35},1.,1.))

    def test_float_tolerance(self):
        self.assertEqual(compare_json({'risk':.05},{'risk':.05+1e-13},1e-9,1e-11),[])

    def test_missing_json_field_fails(self):
        self.assertTrue(compare_json({'n':1,'k':0},{'n':1},1e-9,1e-11))

    def test_conditional_certificate_status_is_checked(self):
        self.assertTrue(compare_json({'certificate_status':'conditional_on_almost_sure_support_reconstruction'},
                                     {'certificate_status':'unconditional'},1e-9,1e-11))

    def test_csv_counts_and_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a.csv';b=Path(tmp)/'b.csv'
            a.write_text('replicate,count,risk\n0,34,0.05\n')
            b.write_text('replicate,count,risk\n0,34,0.05000000000001\n')
            self.assertTrue(compare_csv(a,b)['passed'])
            b.write_text('replicate,count,risk\n0,35,0.05\n')
            self.assertFalse(compare_csv(a,b)['passed'])
            b.write_text('replicate,wrong,risk\n0,34,0.05\n')
            self.assertFalse(compare_csv(a,b)['passed'])


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.workspace=Path(cls.temp.name)/'science'
        reproduce.prepare_workspace(cls.workspace,False)
        sys.path.insert(0,str(cls.workspace/'code'))
        cls.cert=importlib.import_module('certificates')
        cls.horizon=importlib.import_module('horizon_baselines')
        cls.selection=importlib.import_module('selection_audit')

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(cls.workspace/'code'))
        for name,module in list(sys.modules.items()):
            file=getattr(module,'__file__',None)
            if file and Path(file).is_relative_to(cls.workspace):sys.modules.pop(name,None)
        cls.temp.cleanup()

    def test_k_equals_n(self):
        self.assertEqual(self.cert.fixed_bound(5,5,.05),1.)
        self.assertEqual(self.cert.AnytimeCertificate(.05)(5,5),1.)

    def test_zero_complexity_root(self):
        u=self.cert.fixed_bound(100,0,.01)
        self.assertAlmostEqual(u,1-.01**.01,places=14)

    def test_root_is_conservative(self):
        for n,k in [(100,1),(100,5),(1000,20)]:
            u=self.cert.fixed_bound(n,k,.01)
            self.assertLessEqual(self.cert.log_envelope(n,k,u),math.log(.01)+1e-10)

    def test_invalid_inputs(self):
        for args in [(0,0,.05),(10,11,.05),(10,-1,.05),(10,2,0.),(10,2,1.)]:
            with self.subTest(args=args):
                with self.assertRaises(ValueError):self.cert.fixed_bound(*args)

    def test_uniform_complexity_weights(self):
        cert=self.cert.AnytimeCertificate(.05,d=4)
        self.assertAlmostEqual(sum(cert.weight(k) for k in range(5)),1.)
        with self.assertRaises(ValueError):cert.weight(5)

    def test_stricter_confidence(self):
        loose=self.cert.AnytimeCertificate(.05,d=2)(400,2)
        strict=self.cert.AnytimeCertificate(.01,d=2)(400,2)
        self.assertGreaterEqual(strict,loose)

    def test_paper_endpoint(self):
        self.assertAlmostEqual(100000*self.cert.AnytimeCertificate(.05)(100000,2),21.38364321847749,places=10)

    def test_comparator_independent_roots(self):
        for n,k in [(10,0),(50,2),(200,5)]:
            self.assertLess(abs(self.horizon.cg23(n,k,.05)-self.selection.cg23_upper(n,k,.05)),5e-12)

    def test_summable_weights(self):
        n=1000
        partial=math.fsum(self.horizon.time_weight(i) for i in range(1,n+1))
        self.assertAlmostEqual(partial,1-1/math.log(n+math.e),places=14)

    def test_fixed_sample_size(self):
        self.assertEqual(self.cert.classical_sample_size(2,.05,.05),93)

    def test_gcc_stages(self):
        self.assertEqual(self.cert.gcc_refined_stages(2,.05,.05),[92,132,162])

    def test_small_qp_path(self):
        import numpy as np
        from experiments import lmi_path
        n=128;seed=26090800;rng=np.random.default_rng(seed)
        angles=rng.uniform(0,2*np.pi,n);a=rng.uniform(-.7,.7,n)
        fast=self.selection.fast_path(angles,a)
        slow=lmi_path(np.random.default_rng(seed),n,check=True)
        self.assertLess(float(np.max(np.abs(fast[0]-slow['gain']))),2e-10)
        self.assertTrue(np.array_equal(fast[1],slow['k']))


if __name__=='__main__':unittest.main()
