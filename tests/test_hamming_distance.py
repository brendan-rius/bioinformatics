from unittest import TestCase

from hamming_distance import hamming_distance


class TestHamming_distance(TestCase):
    def test_hamming_distance(self):
        """
        Test the case when the strings are the same size
        :return:
        """
        self.assertEqual(hamming_distance("TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC",
                                          "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"), 50)

    def test_hamming_distance_bigger(self):
        """
        Test the case when the second string is bigger than the first one
        """
        self.assertEqual(hamming_distance("GATTCTCA", "GCAAAGACGCTGACCAA"), 3)
