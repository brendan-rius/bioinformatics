from unittest import TestCase

from kmer import all_kmers


class Test_kmers(TestCase):
    def test_all_1mers(self):
        kmers = set(all_kmers(1))
        self.assertEqual(kmers, {'A', 'T', 'G', 'C'})

    def test_all_2mers(self):
        kmers = set(all_kmers(2))
        self.assertEqual(kmers, {'AA', 'AT', 'AG', 'AC',
                                 'TA', 'TT', 'TG', 'TC',
                                 'GA', 'GT', 'GG', 'GC',
                                 'CA', 'CT', 'CG', 'CC'})
