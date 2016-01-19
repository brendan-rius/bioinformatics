import math
import random
from bisect import bisect
from itertools import accumulate

from hamming_distance import hamming_distance
from hash_kmer import hash_nucleotide, unhash_nucleotide
from kmer import kmers
from neighbors import neighbors


def profile(motifs, cromwell=False):
    """
    Compute the profile matrix of a list of motifs
    :param motifs: the list of motifs (a list of strings)
    :param cromwell: whether the probabilities of the profile matrix should follow Cromwell's rule (to use when
    generating profile of small set of sequences). Will avoid to have probabilities equal to zero, but instead use small
    numbers
    :return: its profile matrix with each row being a motif, and each cell in the row being the probability of
    A, C, G, T (respectively)
    """

    def transpose(matrix):
        """
        Transpose a matrix

        Example:
        [['A', 'T', 'G', 'C', 'T"],
         ['A', 'A', 'G', 'C', 'C']]
        becomes
        [['A', 'A'],
         ['T', 'A'],
         ['G', 'G'],
         ['C', 'C'],
         ['T', 'C']]

        :param matrix: the matrix
        :return: the transposed matrix
        """
        return list(zip(*matrix))

    def motifs_list_to_matrix(motifs):
        """
        Transform a motif list to a matrix.

        Example:
        ['ATGCT',
         'AAGCC']
        becomes
        [['A', 'T', 'G', 'C', 'T'],
         ['A', 'A', 'G', 'C', 'C']]

        :return: the matrix
        """
        return [list(motif) for motif in motifs]

    def profile_vector(vector, cromwell):
        """
        Return a distribution vector of the nucleotide vector

        Example:
        [['A', 'T', 'G', 'C', 'T"],
         ['A', 'A', 'G', 'C', 'C']]
        becomes
        [(1/5, 1/5, 1/5, 2/5),
         (1/5, 2/5, 1/5, 1/5)]  (without cromwell)
        or
        [(2/9, 2/9, 2/9, 3/9),
         (2/9, 3/9, 2/9, 2/9)]  (with cromwell)

        :param vector: the vector
        :param cromwell: whether the probabilities of the profile matrix should follow Cromwell's rule (to use when
        generating profile of small set of sequences). Will avoid to have probabilities equal to zero, but instead use small
        numbers
        :return: (a, c, g, t) where a is the frequency of A's, t the frequency of T's, and so on
        """

        if cromwell:
            # Here we add 4 to the length of the vector, because we fake that every nucleotide appeared at least once
            # (to avoid null probabilities).
            len_vec = float(len(vector) + 4)
            return [
                (vector.count('A') + 1) / len_vec,
                (vector.count('C') + 1) / len_vec,
                (vector.count('G') + 1) / len_vec,
                (vector.count('T') + 1) / len_vec
            ]
        else:
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
    profile_matrix = [profile_vector(v, cromwell) for v in transposed]
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


def consensus(motifs):
    """
    Return the consensus string from a list of sequences
    :param motifs: the lost of sequences
    :return: the consensus string
    """
    result = ""
    profile_motifs = profile(motifs)
    for distribution in profile_motifs:
        index = distribution.index(max(distribution))
        result += unhash_nucleotide(index)
    return result


def score(motifs):
    """
    Score a list of motifs found in sequences. The lower the score the better.
    This is done by computing the consensus string of the motifs, and computing the hamming distance between this
    string and the list of sequences.
    :param motifs: the motifs to score found in sequences
    :param sequences: the sequences
    :return: a score in [0; +inf[, the lower being the better
    """
    motifs_consensus = consensus(motifs)
    return hamming_distance(motifs_consensus, motifs)


def greedy_motifs_search(sequences, k, cromwell=True):
    """
    Tries to find a collection of motifs in a collection of sequences of DNA
    :param sequences: the collection of sequences
    :param k: the size of the motifs to search for
    :param cromwell: should we use Cromwell's rule when generating the profile matrix?
    :return: a collection of the most probable motifs (one motif for each sequence)
    """
    best_motifs = None
    for motif1 in kmers(sequences[0], k):
        motifs = [motif1]
        for sequence in sequences[1:]:
            profile_matrix = profile(motifs, cromwell)
            motifs.append(most_probable_kmer_from_profile(sequence, k, profile_matrix))
        if not best_motifs or motifs_entropy(motifs) < motifs_entropy(best_motifs):
            best_motifs = motifs
    return best_motifs


def random_motifs(sequences, k):
    """
    Generate a random list of motifs of length k from a list of sequences
    :param sequences: the sequences
    :param k: the length of the desired motifs
    :return: a random list of motifs
    """
    result = []
    for sequence in sequences:
        index = random.randint(0, len(sequence) - k)
        result.append(sequence[index:index + k])
    return result


def randomized_motifs_search(sequences, k, n=1, cromwell=True):
    """
    Tries to find a list of motifs in a list of sequences of DNA. Tries to smooth the results by running multiple times
    This is a Monte Carlo algorithm.
    :param sequences: the list of sequences
    :param k: the size of the motifs to search for
    :param n: the number of times to run the search (set higher for better results)
    :param cromwell: should we use Cromwell's rule when generating the profile matrix?
    :return: a probable list of the most-probable motifs
    """

    def single_randomized_motifs_search(sequences, k, cromwell=True):
        """
        Tries to find a list of motifs in a list of sequences of DNA.
        This is a Monte Carlo algorithm.
        :param sequences: the list of sequences
        :param k: the size of the motifs to search for
        :param cromwell: should we use Cromwell's rule when generating the profile matrix?
        :return: a probable list of the most-probable motifs
        """
        motifs_list = random_motifs(sequences, k)  # The first list of motifs is randomly sampled from the sequences
        best_motifs = (motifs_entropy(motifs_list), motifs_list)  # (entropy_of_motifs, list_of_motifs)
        while True:  # This algorithm will terminate when the entropy stops improving
            motifs_profile = profile(best_motifs[1], cromwell)  # We get the profile of the current motifs
            motifs_list = motifs(motifs_profile, sequences)  # We generate a new list of motifs based on the profile
            entropy = motifs_entropy(motifs_list)  # We compute the entropy of the new list
            # If the entropy is better, we have a new best motifs list
            if entropy < best_motifs[0]:
                best_motifs = (entropy, motifs_list)
            # If the entropy does not get better, we stop because we do not want to run into an infinite loop
            else:
                return best_motifs

    best_motifs = single_randomized_motifs_search(sequences, k, cromwell)
    for i in range(0, n - 1):  # n - 1, because we already performed a first search
        next_motifs = single_randomized_motifs_search(sequences, k, cromwell)
        if next_motifs[0] < best_motifs[0]:  # If the entropy is better, we have a new best motifs list
            best_motifs = next_motifs
    return best_motifs[1]  # We return only the list of motifs, not the entropy


def motifs(profile, sequences):
    """
    Get the profile-most probable motifs in a list of sequences
    :param profile: the profile matrix
    :param sequences: some sequences of DNA
    :return: a list of the profile most probable motifs in the sequences according
    """

    k = len(profile)  # The profile matrix has the same length as the researched k-mer
    result = []
    for sequence in sequences:
        result.append(most_probable_kmer_from_profile(sequence, k, profile))
    return result


def biased_random(sequence):
    """
    Return a random index from the sequence, biased by the numbers in the sequence
    :param sequence: the sequence of probabilities. If it does not sum to one, it will be adjusted accordingly.
    :return: an index in the sequence
    """

    def sequence_to_distribution(sequence):
        """
        Transform a list of integers in [0; +inf[ to a probability distribution that sums to 1
        :param sequence:
        :return:
        """
        multiplier = 1.0 / sum(sequence)  # We use multiplication instead of division
        return [x * multiplier for x in sequence]

    if sequence not in biased_random.cache:  # If the accumulated distribution is stored in cache
        distribution = sequence_to_distribution(sequence)  # We normalize the sequence
        biased_random.cache[sequence] = list(accumulate(distribution))  # We store the accumulated distribution
    dist = biased_random.cache[sequence]
    return bisect(dist, random.uniform(0, dist[-1]))


biased_random.cache = {}  # Cache for biased random


def profile_random_kmer(profile, sequence):
    """
    Generate a profile-randomly chosen k-mer in a sequence.
    The size of the kmer will be deducted from profile.
    :param profile: the profile
    :param sequence: the sequence to search in
    :return:
    """
    k = len(profile)
    distribution = tuple(probability_from_profile(kmer, profile) for kmer in kmers(sequence, k))
    index = biased_random(distribution)
    return sequence[index:index + k]


def gibbs_motifs_search(sequences, k, n, N=1, cromwell=True):
    """
    Tries to find a list of repeated motifs in a list of sequences.
    This is a Monte Carlo algorithm, but it is different from randomized_motifs_search because it uses
    Gibbs sampling.
    :param sequences: ths list of sequences
    :param k: the wanted size of the motifs
    :param n: the number of iterations
    :param N: the number of time ro run the algorithm
    :return: a list of repeated motifs
    """

    def single_gibbs_motifs_search(sequences, k, n, cromwell=True):
        """
        Tries to find a list of repeated motifs in a list of sequences.
        This is a Monte Carlo algorithm, but it is different from randomized_motifs_search because it uses
        Gibbs sampling.
        :param sequences: ths list of sequences
        :param k: the wanted size of the motifs
        :param n: the number of iterations
        :return: a list of repeated motifs
        """
        motifs_list = random_motifs(sequences, k)  # The first list of motifs is randomly sampled from the sequences
        number_of_motifs = len(motifs_list)
        best_motifs = (motifs_entropy(motifs_list), motifs_list)  # (entropy_of_motifs, list_of_motifs)
        for j in range(n):
            index = random.randint(0, number_of_motifs - 1)  # We choose one of the motifs
            # All the motifs, except the chosen one
            motifs_except_chosen = [m for i, m in enumerate(motifs_list) if i != index]
            # Compute the profile of the motifs (without the chosen motif to avoid biasing the randomness
            # to choose it again)
            motifs_profile = profile(motifs_except_chosen, cromwell)
            # We replace the chosen motif bu a new one which is "randomly" selected but biased by the profile
            motifs_list[index] = profile_random_kmer(motifs_profile, sequences[index])
            entropy = motifs_entropy(motifs_list)
            if entropy < best_motifs[0]:  # If the entropy of the new motifs list is better than the current best
                best_motifs = (entropy, list(motifs_list))
        return best_motifs

    best_motifs = single_gibbs_motifs_search(sequences, k, n, cromwell)
    for i in range(0, N - 1):  # N - 1, because we already performed a first search
        next_motifs = single_gibbs_motifs_search(sequences, k, n, cromwell)
        if next_motifs[0] < best_motifs[0]:  # If the entropy is better, we have a new best motifs list
            best_motifs = next_motifs
    return best_motifs[1]  # We return only the list of motifs, not the entropy
