from kmer import kmers


def motif_enumeration(strings, k, d):
    """
    Check if a motif of length k appears in each sequence in strings with at most d mismatches
    :param strings: the array of sequences
    :param k: the length of the motif
    :param d: the maximum number of mismatches
    :return: the (k, d)-motifs in string as a set
    """
    for sequence in strings:
        for kmer in kmers(sequence, k):
            print(kmer)


def __main__():
    motif_enumeration([
        "ATTTGGC",
        "TGCCTTA",
        "CGGTATC",
        "GAAAATT"], 3, 1)


if __name__ == '__main__':
    __main__()
