import math

from kmer import kmers
from neighbors import neighbors


def profile(motifs):
    """
    Compute the profile matrix of a list of motifs
    :param motifs: the list of motifs (a list of strings)
    :return: its profile matrix with each row being a motif, and each cell in the row being the probability of
    A, T, G, C (respectively)
    """

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

    def profile_vector(vector):
        """
        Return a distribution vector of the nucleotide vector
        :param vector: the vector
        :return: (a, t, g, c) where a is the frequency of A's, t the frequency of T's, and so on
        """
        len_vec = float(len(vector))
        return [
            vector.count('A') / len_vec,
            vector.count('T') / len_vec,
            vector.count('G') / len_vec,
            vector.count('C') / len_vec]

    matrix = motifs_list_to_matrix(motifs)
    transposed = transpose(matrix)
    profile_matrix = [profile_vector(v) for v in transposed]
    return profile_matrix


def motif_entropy(motifs):
    """
    Compute the entropy of a motifs matrix
    :param motifs: the motif matrix
    :return: its entropy
    """

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

    profile_matrix = profile(motifs)  # We get the profile matrix for the motifs
    return sum(entropy(v) for v in profile_matrix)  # We return the sum of entropy for all the motifs


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
    print(motif_entropy(motifs))


def __main2__():
    sequences = """CAATTTACGGAGCCGTCTGATGTTT
CCGACAGGCGCTACCCCTGAGACGT
TTCCCTTCCATCTGAGAAACGCGCA
ACTGACAAATCAGCCAAGTCTCTAG
TTGGCGCAAGAAGAATCTGACATAT
TCTTAGTTTTCCAGGTTTCCCCTGA""".split('\n')
    k = 5
    d = 1
    print(' '.join(list(motif_enumeration(sequences, k, d))))


if __name__ == '__main__':
    __main__()
