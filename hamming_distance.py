def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("The k-mers must have the same length")
    distance = 0
    for idx, nucleotide in enumerate(seq1):
        if seq2[idx] != nucleotide:
            distance += 1
    return distance


def __main__():
    s1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
    s2 = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"
    print(hamming_distance(s1, s2))


if __name__ == '__main__':
    __main__()
