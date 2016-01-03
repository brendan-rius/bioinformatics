import math


def motif_entropy(motifs):
    """
    Compute the entropy of a motifs matrix
    :param motifs: the motif matrix
    :return: its entropy
    """

    def distribution(vector):
        """
        Return a distribution vector of the nucleotide vector
        :param vector: the vector
        :return: (a, t, g, c) where a is the frequency of A's, t the frequency of T's, and so on
        """
        multiplier = 1 / float(len(vector))
        return (
            (vector.count('a') + vector.count('A')) * multiplier,
            (vector.count('t') + vector.count('T')) * multiplier,
            (vector.count('g') + vector.count('G')) * multiplier,
            (vector.count('c') + vector.count('C')) * multiplier)

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

    matrix = motifs_list_to_matrix(motifs)
    transposed = transpose(matrix)
    return sum(entropy(distribution(v)) for v in transposed)


def __main__():
    motifs = ['TCGGGGgTTTtt',
              'cCGGtGAcTTaC',
              'aCGGGGATTTtC',
              'TtGGGGAcTTtt',
              'aaGGGGAcTTCC',
              'TtGGGGAcTTCC',
              'TCGGGGATTcat',
              'TCGGGGATTcCt',
              'TaGGGGAacTaC',
              'TCGGGtATaaCC']
    print(motif_entropy(motifs))


if __name__ == '__main__':
    __main__()
