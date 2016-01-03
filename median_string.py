from hamming_distance import hamming_distance
from kmer import all_kmers


def median_string(k, dna):
    """
    Find the a k-mer x that minimizes the HammingDistance(x, dna)
    :param k: the size of the k-mer to find
    :param dna: a list of DNA sequences
    :return: a k-mer that minimizes the distance between itself and the list of sequences. If multiple k-mers are found,
    return only a single one.
    """
    min_distance = (float("inf"), None)
    for pattern in all_kmers(k):
        distance = hamming_distance(pattern, dna)
        if distance < min_distance[0]:
            min_distance = (distance, pattern)
    return min_distance[1]


def __main__():
    dna = """AGGTACTATCAGTAAAAAGGTTGAACATTTGGTGCGATAGCG
TCAACTCGGTACGCGCATAGCAGGGCAGTTGTACGGGGTGCC
CGTTTAATGACGTTAGCAGGACACCGGTACAGGCTAAAACTA
TGCCAAGTCAACGTTTCCTGGTACGGCTACCGTGTTGGACTC
CGGCCCTGGTACACATGCTTGTGATATCTGGCAATTCTTCAC
CATGGGGGGTACGACAGAAGACATCGATATCCAACTGCGATC
TAACTCAGGTACGTACAGTATTACCGCGGTCTAGTACAAGCT
TTGCCGTGATTTACCCGACACGTGAGGTACTGCTTGCGGAAT
AAACGCGACTGGCCCCACGCTTCAAACAACCTAGGAAGGTAC
TGTTAGGAGAAAACCCATACGCTGAATCGATCGGCCCGGTAC"""
    k = 6
    print(median_string(k, dna.split('\n')))

if __name__ == '__main__':
    __main__()