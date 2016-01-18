from unittest import TestCase

from motifs import profile, probability_from_profile, most_probable_kmer_from_profile, greedy_motifs_search, motifs, \
    randomized_motifs_search, consensus, score, biased_random, profile_random_kmer, gibbs_motifs_search


class TestMotifs(TestCase):
    def test_profile_no_cromwell(self):
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

    def test_profile_cromwell(self):
        """
        Try to create a profile matrix from a list of motifs
        """
        motifs = ['ACCT', 'ATGT', 'ACGG', 'ACGA']
        profile_matrix = profile(motifs, True)
        expected_profile_matrix = [
            [0.625, 0.125, 0.125, 0.125],
            [0.125, 0.5, 0.125, 0.25],
            [0.125, 0.25, 0.5, 0.125],
            [0.25, 0.125, 0.25, 0.375]
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

    def test_greedy_motif_search_no_cromwell(self):
        sequences = [
            'GGCGTTCAGGCA',
            'AAGAATCAGTCA',
            'CAAGGAGTTCGC',
            'CACGTCAATCAC',
            'CAATAATATTCG'
        ]
        k = 3
        motifs = greedy_motifs_search(sequences, k, False)
        self.assertEqual(motifs, ['CAG', 'CAG', 'CAA', 'CAA', 'CAA'])

    def test_greedy_motif_search_cromwell(self):
        sequences = """CCGTGGGGGGTATGTTCCCATTTCAAATGCAGACTAGCCGCGATACCCATCCTGTAACCATGGCTGCAACGTCCCTATCGTGCGGTTTGAATAGGAGAGTGCTATCTAATCGGGGAATGTTGGTTTTATGGTGTGCAACCGCACATTTGAAACCGT
GGTGTGTCGTTCATTGAATGCCGCAACTTCCTGTATGGGGTATTGTTTACACAGTAGATAGGCGAGATCCCCAAGGGCTCGCGTCGGCCGGCTTGTTATGTCTCTTCCTTTCAACCCGGCAGCTGTCGCCGAAGACACTGGACAGCGTATACCAGC
CCTATACTCGGGACCTTTTAGGTACGACTATTCTTGGTACCTATATTAAGCTCAATCGCAATTTTGGGGTGAGTCGAGATTCCCGGCGTAAAGCGATATAACGCATACTATACTAGCAGAAAAGTCCGAACCTCGCGTGATAGTGCAACTAAATTA
GACGCAATTCCCCATATCAGAGCATCCATTGAGCCAACGAACCTGACACAGATAATGTATCATAACGCCTCTTTCACGTCCTACAGCCCTATATTGAACTACGCCTGACACCCGCCAGGAAAACTCAGGAAGTCCCTGTGCCATACGCCACAAGGT
ACTTAATAGTTACTGTCTGTTCGAGGATGAGCACCGTACAAAATGAGCGCTTTACCAACAGGCCCAGCAGATATAGATCAAGCCTAGATGGACAGCTCACCAATGGGCGTCGCCATTCCCGCGTTTGCAGTCAGCAAATGGCTGCAGCAGTAGTGT
TCCTCCCAGGCGAACCAGCTAAGTGGCGAGAAGCCGCATTGCTCTAATTCCTACCCAAGGAAGAGCATCCACCAAAAGCTAAAAAACCCGCCCGAGTGAACATCGCCGGCCGCGATCCCCACAGAAGAGGGCTAGATCTTTATTTAAATTAGTACC
TGGCATACGTAATCTGCACTTCCGACGTCATGGTCCCTTAACCGTTTGGTATACAGGCGGATCCATGGCGGAGCCGTGATGCCCCGGCACACCCTGCCCCCGGACATCTCGCAATACCCCGACATTTTATGCAATAACGCAGCGCATCAGCTTTTC
GTTCACTAAATATACCCCCGGGGTGAGTGGTGGGGGAGGCGCTTCAAGCAGGGTCTGGAGGTGGCGGTCAGTCCCGTTATATGATACCACCTTATAGTCCCCCGGGCCGGTATGACGTAATAACTGCTGCTCGACGGAATGCCCGGTCTCATCGCT
ACCTGGAAGGTCACCTTCAAATGTCCAGTACAGTTGAGCGAGACCTCAGGCTTTATGCGTGCCGGTATGCCCTACCGCGACCACGGACTCAAGATACGCGTATATGAGAGGGTCCTAAACCAATAAATAAAGTGTGGGACAGCTACCTGACAGAAG
GTAAAGCACGATAGGGCTTACTACGAATCGCGTCGTGGCGTTATCCCCCGCTCTATTGGTCGGTAACTCAGACACGCGAGGTGTCAAATACGAGGCCTCTGTAGCTCCAAGTATGATACGGAGGCCTGTTCTCTGCGGCATGATAACTAGATTGCG
ATTTTTGACTGGACCCACGATACAGTCGAGATACCCAAAGCCCTAGAAGTTTTGTGCTTGAGGTGTCGTACATCCATCCGCACCTACATGGCTACAAACACGACGAACACTTAAGTACTCAATGGCCAATAGATTCTTGGGGGTAGGGCCCTGTTC
AGTTCTGCGCCTGAAGAGCCCTAATGTTACTTTGATCAGACGATAGCGCATTTCACCTGAGACAGCGCGGAGGGGCTCACCAAGTGCTTTAAGGATACGGCTCCGACGGACGAAATCCCCAGGTACGGTTCATCATTACGTCGGTGATTACGCGTC
TCGCAGGGCCCAACTATATAATTGTGATTCCGGAACCGACAACAAAGCGTGTATGCTGGTCATCACACTACAATCGGGCGGTACGTGAGTTAGAGAGCAATCGGATCGAAACTTTCCTCTGCCGTTATGCCCGCGATTTGGCCACAACTCGCCTAA
CGAATAAACATTATAAATGCTGAAAGCGTTAGGACGCGAACTCACCCCTTATCGTCTTGAGACTTGTTGCACTCCAACGGGATGGCCGTCATCCCCGGCTGTAAAGAAGCCAGATCGATAGTATGTGCAATGTTAAGGGCCCCTCTTACACAACCG
CTCGATGGGTAGAGATGCCAGTATAAAATGGGTCCTTCGCTTCTGCTAAGGTCACTGCATTGGCGTGAAAATGGTTTCGACTATACGGGATGGTTGGCCGACATCCCCGGCGGCTCTCTTAAGTTTCTTTGTTTTCGAGCGTTGCTAGCGCTCGCA
AAATCCAGAAGAGCGCAGCAGCCCTCCGCTGAGAATAGCTCCTGTTTGCCCTCCCAAGTTTACCCTGCGGCACTCCATCAGGACACCAACATTGATGTCAAATTCAACGACGCGATACCCTATTTGTTATGTTTTTAGAAGCGCTGTAGTCTCAGC
AATGACCAATGATTAACCATGCGAACTTGCTATGTCGACGACATTCCCTAATATGGTCGGAGATGTGGGCAGTAAATGTGGTCAGCCTTCCGATGTATCCAGAACGGATGAATTGGCGTAACTGTGTAGACATCAGTCGGAATGTACCGTTCTGTT
TGGTGCGAACCGTGACCGTCCGGCAAAAAGAAGCCTCCCCCAGTTCAACTGTGATAGCCAGCCGAGATGCCCTCTTGTGCGATTCCCACGTTACGAGCTCATTCTAAATTCGATGAGGGGGGGAATCGCGATCCTCGTAAAGCGAGCGAAGTTAAT
TGGTTAAACAGGAAACGACCAACCAGATTTTCGTCTGTCGGAATTCCCCTTAGTAGAATTAGGAATCGGCCAGACACCAATTTGAAGCGTGATCGCGTATCACCCATGGTATCCTGTAAATTCAATGTCGCGTCCAGAGCACGCGACGTCAATACA
TACGTCTCCAGGATCGCTGTCGACGGTCCCCATATCGCCGCGATGCCCGGTCGGAGCGTAGGAGTGTGGAGAAATATAATGTACTTTGCTGAGCCACGAGCTGTCGTAACCTCCCCTTTAGGACCTCACGGGTTGTTTGGTACCGATTATTTCACC
ATACTTCGGCCGATGCCGCACGACACAAAATAAGGCTCGCGCAGGTCACTGGAAGTATAGCTCGCACAAGAAGCCGGTATCCCCAGCGAGTACGAACTCTTCTATCTGTCATCGAACGGAATGACGTCATTACGCAACAAAGGAATACTCATGGTA
CCGGCGAGATTGACCGAAATCGAAGTGAAGGTCGGGTTGAGGTCAACAAACGGCAAAGACGACGCAATGCCCTCCTGTCGGATCTCTGGGGGAGCTTCCCATCCACCATGCAGATTCGGGACCTCGCAATACTATCGATGTGGTGAATCTGGGCTA
TTCAACCCGAACGAGCCACCGACCCTCGTCTTAGCATAGATTAAACCGGGATCGATCGTGGTACGTCCTACGCGAAATAGGCAAGACCGGTTGAAAGGCCGTCGACTGCCCAGCGGTAAAGGCGCAATCCCCAGTGTTGGTGCCGAGGCGGCGGGC
TTCATGATGACTGCGCCCTGCCACGGCACGGTCCACCCGGGTTAAATGCTACGGGAAGCGGGACAGAAAGTCCCTTCTCTATTTTAGACTTTCCTAACGGTTACTCGGGGCGTTATCCCCTCTTACGACAGCCTGCGCTTTCAGGCATGCTATCAG
TCATATTAATTGACCTCAACCACATCTCTTTAAGACGATTACCATAGCGGTCAAGTCAACGACGAGATTCCCCGGCGTCATGTGTCATTAATAGTAAGGGCCTACTTAGGATGGAAACATCACAGAAACAGACCAGCGTTCTATCGCGAGCAACCA
""".upper().split()
        k = 12
        motifs = greedy_motifs_search(sequences, k)
        self.assertEqual(motifs, ['GCCGCGATACCC',
                                  'GGCGAGATCCCC',
                                  'GTCGAGATTCCC',
                                  'GACGCAATTCCC',
                                  'GTCGCCATTCCC',
                                  'GCCGCGATCCCC',
                                  'GCCGTGATGCCC',
                                  'GACGGAATGCCC',
                                  'GCCGGTATGCCC',
                                  'GGCGTTATCCCC',
                                  'GTCGAGATACCC',
                                  'GACGAAATCCCC',
                                  'GCCGTTATGCCC',
                                  'GCCGTCATCCCC',
                                  'GCCGACATCCCC',
                                  'GACGCGATACCC',
                                  'GACGACATTCCC',
                                  'GCCGAGATGCCC',
                                  'GTCGGAATTCCC',
                                  'GCCGCGATGCCC',
                                  'GCCGGTATCCCC',
                                  'GACGCAATGCCC',
                                  'GGCGCAATCCCC',
                                  'GGCGTTATCCCC',
                                  'GACGAGATTCCC'])

    def test_motifs(self):
        # If we test with a matrix filled with zeros, it wont find 'ATA' because there would be a zero probability for
        # the 'T' to appear. Instead we follow Cromwell's rule, and avoid zeros and use small numbers instead
        profile = [
            [0.97, 0.01, 0.01, 0.01],
            [0.97, 0.01, 0.01, 0.01],
            [0.97, 0.01, 0.01, 0.01],
        ]
        sequences = [
            'ATGCAAATGC',
            'ATGCATATGC',
            'ATGCAAATGC'
        ]
        expected_motifs = [
            'AAA',
            'ATA',
            'AAA',
        ]
        self.assertEqual(motifs(profile, sequences), expected_motifs)

    def test_randomized_motifs_search(self):
        """
        Might fail because of randomness, so run it mutliple time if not sure
        """
        sequences = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                     'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                     'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

        k = 8
        self.assertEqual(randomized_motifs_search(sequences, k, 5000),
                         ['GTAAACGG', 'GTGTAAGT', 'GTATACAG', 'GTGCACGT', 'GTGCAATG'])

    def test_consensus(self):
        sequences = [
            'AAATACAGACAGCGT',
            'AAAAAATAGCAGGGT',
            'TAAAATAAACAGCGG',
            'ACAGAAAAAAAGGGG',
            'AAAATAAAACTGCGA'
            'ATAGACGAACACGGT',
            'CAAAAAGAGAAGGGG',
            'ATAGAAAAGGAAGGG',
            'AAGAAAAAAGAGAGG',
            'CATAATGAACTGTGA'
        ]
        self.assertEqual(consensus(sequences), 'AAAAAAAAACAGGGG')

    def test_score(self):
        sequences = [
            'AAATACAGACAGCGT',
            'AAAAAATAGCAGGGT',
            'TAAAATAAACAGCGG',
            'ACAGAAAAAAAGGGG',
            'AAAATAAAACTGCGA',
            'ATAGACGAACACGGT',
            'CAAAAAGAGAAGGGG',
            'ATAGAAAAGGAAGGG',
            'AAGAAAAAAGAGAGG',
            'CATAATGAACTGTGA'
        ]
        self.assertEqual(score(sequences), 43)

    def test_biased_random(self):
        sequence = tuple(range(10))  # The sequence is (0, 1, ..., 9)
        picks = [biased_random(sequence) for i in range(10000)]
        counts = [picks.count(i) for i in sequence]  # We count the occurrences of each number of the sequence
        # It should be sorted since 1 should appear less frequently than 2, and 2 than 3, ...
        self.assertEqual(sorted(counts), counts)

    def test_profile_random_kmer(self):
        profile = [
            [0.97, 0.01, 0.01, 0.01],
            [0.97, 0.01, 0.01, 0.01],
            [0.49, 0.01, 0.01, 0.49],
            [0.97, 0.01, 0.01, 0.01],
        ]  # Profile for "AAAA" or "AATA", according to Cromwell's rule
        sequence1 = "ATGCATGCATGCAATAATGCATGCGCTAGATC"  # AATA is contained in the string
        sequence2 = "ATGCATGCATGCAAAATGCATGCGCTAGATC"  # AAAA is contained in the string
        random_kmer1 = profile_random_kmer(profile, sequence1)
        random_kmer2 = profile_random_kmer(profile, sequence2)
        self.assertEqual(random_kmer1, "AATA")
        self.assertEqual(random_kmer2, "AAAA")

    def test_gibbs_motifs_search(self):
        sequences = [
            'CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
            'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
            'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
            'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
            'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
        ]
        motifs = gibbs_motifs_search(sequences, k=8, n=100, N=20)
        print('\n'.join(motifs))
        self.assertLessEqual(score(motifs), 9)
