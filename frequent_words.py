from collections import Counter
from count_hamming import count_hamming
from hash_kmer import hash_kmer
from kmer import kmer
from neighbors import neighbors
from reverse_complement import reverse_complement


def frequent_words_mismatch(genome, k, distance, reverse=False):
    """
    Return a list of the most frequent words in the genome with mismatches
    :param genome:the genome
    :param k: the length of the k-mer
    :param distance: the maximum distance
    :param reverse: should we take account of the reverse complements too?
    :return:
    """
    kmers = set()
    max_index = len(genome) - k
    result = Counter()
    for i in range(0, max_index + 1):
        kmer = genome[i:i + k]
        kmers = kmers.union(neighbors(kmer, distance))
    for kmer in kmers:
        result[kmer] += count_hamming(genome, kmer, distance)
        if reverse:
            rev_kmer = reverse_complement(kmer)
            result[kmer] += count_hamming(genome, rev_kmer, distance)
    maxCount = max(result.values())
    frequent_kmers = {word: frequency for word, frequency in result.items() if frequency == maxCount}
    return frequent_kmers


def frequent_words(text, k):
    """
    Return the most frequent k-mers in a text
    :param text: the text
    :param k: the length of the k-ners
    :param with_count: flag to ask the function to return the nu;ber of apparitions of the wor or not
    :return: a set of most frequent k-mers
    """
    frequency_map = frequency_kmer(text, k)
    if len(frequency_map) == 0:
        return {}
    maxCount = max(frequency_map.values())
    frequent_kmers = {word: frequency for word, frequency in frequency_map.items() if frequency == maxCount}
    return frequent_kmers


def frequency_kmer(sequence, k):
    """
    Returns the frequency of all k-mers found in a sequence
    :param sequence: the sequence
    :param k: the size of k-mers
    :return: a counter with k-mer as keys and number of apparitions as value
    """
    if k == 0:
        raise ValueError("k has to be positive")
    len_text = len(sequence)
    max_index = len_text - k
    if len_text == 0:
        return {}
    frequency_map = Counter()
    for i in range(0, max_index + 1):
        foundKmer = kmer(sequence, position=i, k=k)
        frequency_map[foundKmer] += 1
    return frequency_map


def __main__():
    text = "atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccgacccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaagggggggatgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccgagctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggagatcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttataggtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcacccaaattcagtgtgggcgagcgcaacggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcataacttgagttaaaaaaaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgtattggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcataaaaaaaagggggggaccgaaagggaagctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga"
    k = 15
    print(' '.join(frequent_words(text, k)))


def __main2__():
    genome = "AGGCCAGGAACCACACTGCAGGTGCAGGTGCTGCAGGAAAGGTGCCCAGGAGGAGGAAACACACACCCTGCAATGCAGGCCACAAAGGAACCCCACAACCACAAACAAAGGACAAAGGAATGCAGGAAACAAACAAAGGACAGGAGGACACAAAGGCCAGGAACCACAGGCCACCCCCCCCCAAACAGGCCACTGCAGGAATGCCCTGCAGGTGCAC"
    k = 6
    distance = 2
    print(' '.join(frequent_words_mismatch(genome, k, distance, True).keys()))


def __main3__():
    text = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
    kstart = 3
    kend = 9
    for k in range(kstart, kend + 1):
        print((k, frequent_words(text, k)))


if __name__ == '__main__':
    __main__()
