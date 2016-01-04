import math

from kmer import kmers
from neighbors import neighbors


def motif_entropy(motifs):
    """
    Compute the entropy of a motifs matrix
    :param motifs: the motif matrix
    :return: its entropy
    """

    def distribution(vector):
        """
        Return a distribution vector of the nucleotide vector
        :param vector: the vector
        :return: (a, t, g, c) where a is the frequency of A's, t the frequency of T's, and so on
        """
        multiplier = 1 / float(len(vector))
        return (
            (vector.count('a') + vector.count('A')) * multiplier,
            (vector.count('t') + vector.count('T')) * multiplier,
            (vector.count('g') + vector.count('G')) * multiplier,
            (vector.count('c') + vector.count('C')) * multiplier)

    def entropy(vector):
        """
        Compute the entropy of a distribution vector
        :param vector: the vector
        :return: the entropy (between 0 and 2)
        """
        vector_entropy = 0
        for x in vector:
            if x != 0:  # We avoid to compute log_2(0), and we consider it to be zero
                vector_entropy += x * math.log(x, 2)
        return -vector_entropy

    def transpose(matrix):
        """
        Transpose a matrix
        :param matrix: the matrix
        :return: the transposed matrix
        """
        return zip(*matrix)

    def motifs_list_to_matrix(motifs):
        """
        Transform a motif list to a matrix
        :return:
        """
        return [list(motif) for motif in motifs]

    matrix = motifs_list_to_matrix(motifs)
    transposed = transpose(matrix)
    return sum(entropy(distribution(v)) for v in transposed)


def motif_enumeration(sequences, k, d):
    """
    Check if a motif of length k appears in each sequence in strings with at most d mismatches
    :param sequences: the array of sequences
    :param k: the length of the motif
    :param d: the maximum number of mismatches
    :return: the (k, d)-motifs in string as a set
    """
    motifs = set()
    for kmer in kmers(sequences, k):
        neighborhood = neighbors(kmer, d)
        for neighbor in neighborhood:
            neighborhood2 = neighbors(neighbor, d)
            if all(any(neighbor2 in seq for neighbor2 in neighborhood2) for seq in sequences):
                motifs.add(neighbor)
    return motifs


def __main__():
    sequences = """CAATTTACGGAGCCGTCTGATGTTT
CCGACAGGCGCTACCCCTGAGACGT
TTCCCTTCCATCTGAGAAACGCGCA
ACTGACAAATCAGCCAAGTCTCTAG
TTGGCGCAAGAAGAATCTGACATAT
TCTTAGTTTTCCAGGTTTCCCCTGA""".split('\n')
    k = 5
    d = 1
    print(' '.join(list(motif_enumeration(sequences, k, d))))


def __main2__():
    motifs = ['TCGGGGgTTTtt',
              'cCGGtGAcTTaC',
              'aCGGGGATTTtC',
              'TtGGGGAcTTtt',
              'aaGGGGAcTTCC',
              'TtGGGGAcTTCC',
              'TCGGGGATTcat',
              'TCGGGGATTcCt',
              'TaGGGGAacTaC',
              'TCGGGtATaaCC']
    print(motif_entropy(motifs))


if __name__ == '__main__':
    __main__()
