import csv
import gzip
import importlib
import io
import sys
import tempfile
import unittest
from collections import Counter
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

import pysam

from scCirRL.assign_bc import assign_ref_bc
from scCirRL.cluster_umi import umi_clustering
from scCirRL.parameter import _five_ada, scrl_para_class
from scCirRL.scCirRL_assign_bc_from_tsv import assign_bc_from_tsv
from scCirRL.scCirRL_merge_bc_umi import merge_bc_umi_tsvs
from scCirRL.umi_stats import BC_UMI_FIELDS, BC_UMI_HEADER, format_umi_clustering_stats


BASE_UMI = 'ACGTACGTACGT'
BRIDGE_UMI = 'ACGTTCGTACGT'
DISTANT_UMI = 'ACGTTCGAACGT'
BC1, BC2 = 'ACGT' * 4, 'TGCA' * 4


class UMIClusteringStatsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # Repeated sequences occur across reads, transcripts, and cell barcodes.
        observations = (
            [(BC1, BASE_UMI, 'T1', 'chr1')] * 4
            + [(BC1, BRIDGE_UMI, 'T1', 'chr1')] * 2
            + [(BC1, DISTANT_UMI, 'T1', 'chr1')]
            + [(BC1, BASE_UMI, 'T2', 'chr1')] * 3
            + [(BC2, BASE_UMI, 'T1', 'chr2')] * 2
        )
        self.reads = [
            ('r{}'.format(i), bc, umi, tx, chrom)
            for i, (bc, umi, tx, chrom) in enumerate(observations)
        ]
        self.read_to_trans = {r: [tx] for r, bc, umi, tx, chrom in self.reads}
        self.read_to_cate = {r: 'FSM' for r, bc, umi, tx, chrom in self.reads}
        self.genes = {'T1': {'id': 'G1', 'name': 'Gene1'},
                      'T2': {'id': 'G2', 'name': 'Gene2'}}

    def write_reads_tsv(self, reads=None):
        path = self.root / 'reads.tsv'
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(['#qname', 'chrom', 'start', 'end', 'bc', 'umi',
                             'bu', 'is_perfect', 'ts_tag'])
            for read, bc, umi, tx, chrom in self.reads if reads is None else reads:
                writer.writerow([read, chrom, 100, 200, bc, umi,
                                 bc + umi + 'TTTTT', 1, '-'])
        return path

    def run_tsv_pipeline(self, reads=None, name='chunk.tsv', **kwargs):
        input_path = self.write_reads_tsv(reads)
        output_path = self.root / name
        log = io.StringIO()
        with redirect_stderr(log):
            result = assign_bc_from_tsv(
                str(input_path), [BC1, BC2], BC1 + 'NNN' + BC2,
                self.genes, self.read_to_trans, self.read_to_cate,
                str(output_path), 16, 2, 12, 1, **kwargs)
        return output_path, result, log.getvalue()

    def read_clusters(self, path):
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def write_clusters(self, path, rows, legacy=False):
        opener = gzip.open if str(path).endswith('.gz') else open
        fields = BC_UMI_FIELDS[:-1] if legacy else BC_UMI_FIELDS
        with opener(path, 'wt', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t',
                                    lineterminator='\n', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

    def assert_fixture_counts(self, rows):
        counts = {(row['#CellBarcode'], row['TranscriptID']):
                  (int(row['ReadCount']), int(row['RawUMICount'])) for row in rows}
        self.assertEqual(counts, {(BC1, 'T1'): (7, 3), (BC1, 'T2'): (3, 1),
                                  (BC2, 'T1'): (2, 1)})

    def test_graph_counts_sequences_not_reads_or_nodes(self):
        observations = [(r, bc, umi, True) for r, bc, umi, tx, chrom in self.reads]
        assigned, clusters, to_top, to_parent = umi_clustering(
            '', self.read_to_trans, self.genes, observations, 1)
        self.assertEqual(len(assigned), 12)
        self.assertEqual(sorted((len(row[2]), row[6]) for row in clusters),
                         [(2, 1), (3, 1), (7, 3)])
        self.assertEqual(assigned['r6'][-2:], (2, 1))
        self.assertEqual(sum(to_top.values()), 12)
        self.assertEqual(sum(to_parent.values()), 12)

    def test_tsv_pipeline_reports_per_cluster_and_run_counts(self):
        output, result, log = self.run_tsv_pipeline()
        self.assertEqual(result, (12, 12))
        self.assert_fixture_counts(self.read_clusters(output))
        self.assertIn(format_umi_clustering_stats(Counter({1: 2, 3: 1})), log)
        self.assertIn('1 (33.33% of counted clusters)', log)
        self.assertIn('Additional raw UMIs merged        : 2', log)

    def test_empty_tsv_reports_zero_clusters(self):
        output, result, log = self.run_tsv_pipeline(reads=[])
        self.assertEqual(result, (0, 0))
        self.assertEqual(output.read_text(), BC_UMI_HEADER)
        self.assertIn(format_umi_clustering_stats(Counter()), log)

    def test_bam_pipeline_aggregates_regions_and_matches_tsv(self):
        bam_path = self.root / 'input.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': 2000}, {'SN': 'chr2', 'LN': 2000}]}
        with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as bam:
            for i, (name, bc, umi, tx, chrom) in enumerate(self.reads):
                prefix = _five_ada + bc + umi + 'T' * 10
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = prefix + 'CG' * 50
                read.flag = 0
                read.reference_id = 0 if chrom == 'chr1' else 1
                read.reference_start = 100 + i
                read.mapping_quality = 60
                read.cigartuples = [(4, len(prefix)), (0, 100)]
                read.set_tag('ts', '-')
                bam.write(read)
        pysam.index(str(bam_path))
        para = scrl_para_class()
        para.umi_len = 12
        para.long_bam = str(bam_path)
        para.out_bu_fn = str(self.root / 'main.tsv')
        para.out_bu_bam = str(self.root / 'main.bam')
        para.log_fn = str(self.root / 'main.log')
        with redirect_stderr(io.StringIO()):
            assign_ref_bc([('chr1', 0, 1000), ('chr2', 0, 1000)], 12,
                          [BC1, BC2], BC1 + 'NNN' + BC2,
                          self.genes, self.read_to_trans, self.read_to_cate, para)
        rows = self.read_clusters(para.out_bu_fn)
        self.assert_fixture_counts(rows)
        self.assertIn(format_umi_clustering_stats(Counter({1: 2, 3: 1})),
                      Path(para.log_fn).read_text())
        tsv, _, _ = self.run_tsv_pipeline(
            in_bam_fn=str(bam_path), out_bam_fn=str(self.root / 'chunk.bam'))
        self.assertEqual(rows, self.read_clusters(tsv))
        for path in [para.out_bu_bam, str(self.root / 'chunk.bam')]:
            with pysam.AlignmentFile(path) as bam:
                reads = list(bam.fetch())
            self.assertEqual(len(reads), 12)
            self.assertTrue(all(read.has_tag('CB') and read.has_tag('UB') for read in reads))

    def test_merge_sums_cluster_histograms(self):
        source, _, _ = self.run_tsv_pipeline()
        rows = self.read_clusters(source)
        first, second = self.root / 'first.tsv.gz', self.root / 'second.tsv'
        self.write_clusters(first, rows[:2])
        self.write_clusters(second, rows[2:])
        output = self.root / 'merged.tsv.gz'
        log = io.StringIO()
        with redirect_stderr(log):
            result = merge_bc_umi_tsvs([str(first), str(second)], str(output))
        self.assertEqual(result, (2, 3))
        self.assertEqual(rows, self.read_clusters(output))
        self.assertIn(format_umi_clustering_stats(Counter({1: 2, 3: 1})), log.getvalue())

    def test_merge_does_not_guess_counts_for_legacy_chunks(self):
        source, _, _ = self.run_tsv_pipeline()
        rows = self.read_clusters(source)
        legacy = self.root / 'legacy.tsv'
        self.write_clusters(legacy, rows, legacy=True)
        output = self.root / 'mixed.tsv'
        log = io.StringIO()
        with redirect_stderr(log):
            merge_bc_umi_tsvs([str(legacy), str(source)], str(output))
        merged = self.read_clusters(output)
        self.assertEqual([row['RawUMICount'] for row in merged[:3]], ['NA'] * 3)
        self.assertEqual(merged[3:], rows)
        self.assertIn(format_umi_clustering_stats(Counter({1: 2, 3: 1}), 3), log.getvalue())

    def test_merge_rejects_invalid_raw_counts(self):
        source, _, _ = self.run_tsv_pipeline()
        row = self.read_clusters(source)[0]
        for invalid in ['0', '-1', 'not-a-count', str(int(row['ReadCount']) + 1)]:
            with self.subTest(invalid=invalid):
                path = self.root / 'invalid.tsv'
                self.write_clusters(path, [dict(row, RawUMICount=invalid)])
                with self.assertRaises(ValueError):
                    merge_bc_umi_tsvs([str(path)], str(self.root / 'invalid-merged.tsv'))

    def test_merge_preserves_long_read_name_fields(self):
        names = ['read{:032d}'.format(i) for i in range(4000)]
        row = dict(zip(BC_UMI_FIELDS, [
            BC1, BASE_UMI, '4000', 'T1', 'G1', 'Gene1', ','.join(names),
            ','.join(['FSM'] * 4000), ','.join(['N'] * 4000), '1']))
        self.assertGreater(len(row['ReadNames']), 131072)
        source, output = self.root / 'deep.tsv', self.root / 'deep-merged.tsv'
        self.write_clusters(source, [row])
        with redirect_stderr(io.StringIO()):
            merge_bc_umi_tsvs([str(source)], str(output))
        self.assertEqual(self.read_clusters(output), [row])

    def test_expression_readers_accept_old_and_new_columns(self):
        # Table parsing does not call R; isolate the unrelated R import dependency.
        r_modules = {name: mock.MagicMock() for name in [
            'rpy2', 'rpy2.robjects', 'rpy2.robjects.packages']}
        with mock.patch.dict(sys.modules, r_modules):
            matrix = importlib.import_module('scCirRL.scCirRL_make_expression_matrix')
        source, _, _ = self.run_tsv_pipeline()
        legacy = self.root / 'legacy.tsv'
        self.write_clusters(legacy, self.read_clusters(source), legacy=True)
        self.assertEqual(matrix.get_bu_cmpt(str(source), [], False, {}),
                         matrix.get_bu_cmpt(str(legacy), [], False, {}))
        self.assertEqual(matrix.get_bu_hap_cmpt(str(source), [], False, {}, {}),
                         matrix.get_bu_hap_cmpt(str(legacy), [], False, {}, {}))
        old_cells, old_paths, old_trans = matrix.get_bu_cmpt_streaming(
            str(legacy), [], False, {}, str(self.root / 'old_stream'))
        new_cells, new_paths, new_trans = matrix.get_bu_cmpt_streaming(
            str(source), [], False, {}, str(self.root / 'new_stream'))
        self.assertEqual((old_cells, old_trans), (new_cells, new_trans))
        self.assertEqual({k: Path(v).read_text() for k, v in old_paths.items()},
                         {k: Path(v).read_text() for k, v in new_paths.items()})


if __name__ == '__main__':
    unittest.main()
