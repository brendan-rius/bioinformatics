from unittest import TestCase

from median_string import median_string


class TestMedian_string(TestCase):
    def test_median_string(self):
        self.assertEqual(median_string(3, ['AAATTGACGCAT',
                                           'GACGACCACGTT',
                                           'CGTCAGCGCCTG',
                                           'GCTGAGCACCGG',
                                           'AGTTCGGGACAG']), 'GAC')
