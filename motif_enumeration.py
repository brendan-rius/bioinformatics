from kmer import kmers
from neighbors import neighbors


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


if __name__ == '__main__':
    __main__()
