from unittest import TestCase

from motifs import profile


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
            [0.2, 0.7, 0, 0.1],
            [0.2, 0.2, 0, 0.6],
            [0, 0, 1, 0],
            [0, 0, 1, 0],
            [0, .1, 0.9, 0],
            [0, .1, 0.9, 0],
            [0.9, 0, 0.1, 0],
            [0.1, 0.5, 0, 0.4],
            [0.1, 0.8, 0, 0.1],
            [0.1, 0.7, 0, 0.2],
            [0.3, 0.3, 0, 0.4],
            [0.0, 0.4, 0, 0.6]
        ]
        self.assertEqual(profile_matrix, expected_profile_matrix)
