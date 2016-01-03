from unittest import TestCase

from kmer import all_kmers, kmers


class TestKmers(TestCase):
    def test_all_1mers(self):
        """
        Try to generate all possible 1-mers
        """
        kmers = set(all_kmers(1))
        self.assertEqual(kmers, {'A', 'T', 'G', 'C'})

    def test_all_2mers(self):
        """
        Try to generate all possible 2-mers
        """
        kmers = set(all_kmers(2))
        self.assertEqual(kmers, {'AA', 'AT', 'AG', 'AC',
                                 'TA', 'TT', 'TG', 'TC',
                                 'GA', 'GT', 'GG', 'GC',
                                 'CA', 'CT', 'CG', 'CC'})

    def test_kmers(self):
        """
        Iterate over all the possible 2-mers of a string
        """
        s = "ATTGCTCA"
        possible_2mers = list(kmers(s, 2))
        self.assertEqual(possible_2mers, ['AT', 'TT', 'TG', 'GC', 'CT', 'TC', 'CA'])
