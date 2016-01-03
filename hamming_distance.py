from kmer import kmers


def hamming_distance(seq1, seq2):
    """
    Compute the hamming distance between two strings.
    If the second string is bigger than the first one, this will return the minimum hamming distance it found between
    the first k-mer, and all the possible k-mers in the second string.
    :param seq1: the first string
    :param seq2: the second string (can be bigger than the first one)
    :return: the hamming distance between the two (an integer), or the minimum hamming distance if the second sequence
    is bigger than the first one.
    """

    def hamming_distance_same_size(s1, s2):
        """
        Compute the hamming distance of two strings that have the same size
        :param s1: the fist string
        :param s2: the second string
        :return: the hamming distance between the two strings
        """
        distance = 0
        for c1, c2 in zip(s1, s2):
            if c1 != c2:  # We compare the two string char-by-char
                distance += 1
        return distance

    k = len(seq1)
    # We compute the hamming distance between seq1 and s for s being all the possibles strings the same size as seq1
    # in seq2.
    distances = [hamming_distance_same_size(seq1, s) for s in kmers(seq2, k)]
    return min(distances)  # We return the minimum found hamming distance


def __main__():
    s1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
    s2 = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"
    print(hamming_distance(s1, s2))


if __name__ == '__main__':
    __main__()
