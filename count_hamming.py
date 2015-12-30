from hamming_distance import hamming_distance
from reverse_complement import reverse_complement


def count_hamming(genome, kmer, distance):
    """
    Return the number of occurrences of a given sequence and its similar sequences in a genome
    :param genome: the genome
    :param kmer: the sequence
    :param distance: the maximum hamming distance for 2 sequence to be similar
    :return: the number of occurrences
    """
    k = len(kmer)
    len_genome = len(genome)
    if k > len_genome:
        return 0
    count = 0
    max_index = len_genome - k
    for i in range(0, max_index + 1):
        word = genome[i:i + k]
        if hamming_distance(kmer, word) <= distance:
            count += 1
    return count


def __main__():
    genome = "CATGCCATTCGCATTGTCCCAGTGA"
    kmer = "CCC"
    distance = 2

    print(count_hamming(genome, kmer, distance))


if __name__ == '__main__':
    __main__()
