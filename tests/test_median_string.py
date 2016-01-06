from unittest import TestCase

from median_string import median_string, median_strings


class TestMedianString(TestCase):
    def test_median_string(self):
        dna = dna = ['AAATTGACGCAT',
                     'GACGACCACGTT',
                     'CGTCAGCGCCTG',
                     'GCTGAGCACCGG',
                     'AGTTCGGGACAG']
        self.assertEqual(median_string(3, dna), 'GAC')

    def test_median_strings(self):
        dna = ['AAATTGACGTTT',
               'GACGAAACCTTT',
               'CAAAAGCGTTTG',
               'GCTTTGACAAAC',
               'AGTTTCGGAAAG']
        self.assertEqual(median_strings(3, dna), ['TTT', 'AAA'])
