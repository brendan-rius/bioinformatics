import math

from hash_kmer import hash_nucleotide
from kmer import kmers
from neighbors import neighbors


def profile(motifs):
    """
    Compute the profile matrix of a list of motifs
    :param motifs: the list of motifs (a list of strings)
    :return: its profile matrix with each row being a motif, and each cell in the row being the probability of
    A, C, G, T (respectively)
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
        :return: (a, c, g, t) where a is the frequency of A's, t the frequency of T's, and so on
        """
        # Could be faster to compute 1.0 / len(vector) and to multiply instead of dividing, but gives floating point
        # imprecision. Probably automatically optimized anyway.
        len_vec = float(len(vector))

        return [
            vector.count('A') / len_vec,
            vector.count('C') / len_vec,
            vector.count('G') / len_vec,
            vector.count('T') / len_vec
        ]

    matrix = motifs_list_to_matrix(motifs)
    transposed = transpose(matrix)
    profile_matrix = [profile_vector(v) for v in transposed]
    return profile_matrix


def probability_from_profile(sequence, profile_matrix):
    """
    Compute the probability for a sequence to be generated knowing a profile matrix
    :param sequence: the sequence to compute the probability for
    :param profile_matrix: the profile matrix
    :return: the probability that the given profile matrix generates this sequence
    """
    probability = 1
    for nucleotide, profile_vector in zip(sequence, profile_matrix):
        probability_of_nucleotide = profile_vector[hash_nucleotide(nucleotide)]
        probability *= probability_of_nucleotide
    return probability


def motifs_entropy(motifs):
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


def motifs_enumeration(sequences, k, d):
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


def most_probable_kmer_from_profile(sequence, k, profile_matrix):
    """
    Find the most profile-probable k-mer in a sequence
    :param sequence: the sequence
    :param k: the size of the k-mer
    :param profile_matrix: the profile matrix
    :return: the most probable k-mer in sequence according to the given profile matrix
    """
    most_probable = (-1, None)
    for kmer in kmers(sequence, k):
        probability = probability_from_profile(kmer, profile_matrix)
        if probability > most_probable[0]:
            most_probable = (probability, kmer)
    return most_probable[1]


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
    print(motifs_entropy(motifs))


def __main2__():
    sequences = """CAATTTACGGAGCCGTCTGATGTTT
CCGACAGGCGCTACCCCTGAGACGT
TTCCCTTCCATCTGAGAAACGCGCA
ACTGACAAATCAGCCAAGTCTCTAG
TTGGCGCAAGAAGAATCTGACATAT
TCTTAGTTTTCCAGGTTTCCCCTGA""".split('\n')
    k = 5
    d = 1
    print(' '.join(list(motifs_enumeration(sequences, k, d))))


if __name__ == '__main__':
    __main__()
