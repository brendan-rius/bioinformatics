from unittest import TestCase

from motifs import profile, probability_from_profile, most_probable_kmer_from_profile


class TestMotifs(TestCase):
    def test_profile(self):
        """
        Try to create a profile matrix from a list of motifs
        """
        motifs = ['TCGGGGGTTTTT',
                  'CCGGTGACTTAC',
                  'ACGGGGATTTTC',
                  'TTGGGGACTTTT',
                  'AAGGGGACTTCC',
                  'TTGGGGACTTCC',
                  'TCGGGGATTCAT',
                  'TCGGGGATTCCT',
                  'TAGGGGAACTAC',
                  'TCGGGTATAACC']
        profile_matrix = profile(motifs)
        expected_profile_matrix = [
            [0.2, 0.1, 0, 0.7],
            [0.2, 0.6, 0, 0.2],
            [0, 0, 1, 0],
            [0, 0, 1, 0],
            [0, 0, 0.9, 0.1],
            [0, 0, 0.9, 0.1],
            [0.9, 0, 0.1, 0],
            [0.1, 0.4, 0, 0.5],
            [0.1, 0.1, 0, 0.8],
            [0.1, 0.2, 0, 0.7],
            [0.3, 0.4, 0, 0.3],
            [0.0, 0.6, 0, 0.4]
        ]
        self.assertEqual(profile_matrix, expected_profile_matrix)

    def test_probability_from_profile(self):
        profile = [
            [0.2, 0.1, 0, 0.7],
            [0.2, 0.6, 0, 0.2],
            [0, 0, 1, 0],
            [0, 0, 1, 0],
            [0, 0, 0.9, 0.1],
            [0, 0, 0.9, 0.1],
            [0.9, 0, 0.1, 0],
            [0.1, 0.4, 0, 0.5],
            [0.1, 0.1, 0, 0.8],
            [0.1, 0.2, 0, 0.7],
            [0.3, 0.4, 0, 0.3],
            [0.0, 0.6, 0, 0.4]
        ]
        kmer = "TCGGGGATTTCC"
        self.assertEqual(probability_from_profile(kmer, profile), 0.020575296)

    def test_most_probable_kmer_from_profile(self):
        profile_matrix = [
            (0.2, 0.4, 0.3, 0.1),
            (0.2, 0.3, 0.3, 0.2),
            (0.3, 0.1, 0.5, 0.1),
            (0.2, 0.5, 0.2, 0.1),
            (0.3, 0.1, 0.4, 0.2)
        ]
        k = 5
        sequence = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
        self.assertEqual(most_probable_kmer_from_profile(sequence, k, profile_matrix), "CCGAG")
