from unittest import TestCase

from hash_kmer import hash_nucleotide, unhash_nucleotide, hash_kmer, unhash_kmer


class TestHash(TestCase):
    def test_hash_nucleotide(self):
        self.assertEqual(hash_nucleotide('A'), 0)
        self.assertEqual(hash_nucleotide('C'), 1)
        self.assertEqual(hash_nucleotide('G'), 2)
        self.assertEqual(hash_nucleotide('T'), 3)

    def test_unhash_nucleotide(self):
        self.assertEqual(unhash_nucleotide(0), 'A')
        self.assertEqual(unhash_nucleotide(1), 'C')
        self.assertEqual(unhash_nucleotide(2), 'G')
        self.assertEqual(unhash_nucleotide(3), 'T')

    def test_hash_kmer(self):
        self.assertEqual(hash_kmer("TCTGATACTCTGTGGG"), 3727810282)

    def test_unhash_kmer(self):
        self.assertEqual(unhash_kmer(3727810282, 16), "TCTGATACTCTGTGGG")
