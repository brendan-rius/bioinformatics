def probability_kmer(s, n, k):
    """
    Compute the probability for a random k-mer to appear at least n times in a random n-long DNA sequence
    :param s: the length of the sequence
    :param n: the minimum number of apparitions of the k-mer
    :param k: the siwe of the k-mer
    :return: the probability
    """
    # See http://math.stackexchange.com/questions/1548725/probability-for-a-random-k-mer-to-repeat-at-least-n-times-in-a-random-s-char
    res_s = 0
    for j in range(n, s // k + 1):
        res_p = 1
        for i in range(1, j + 1):
            res_p *= s - i * k + 1
        res_p *= 4 ** (k - k * j)
        res_s += res_p
    return res_s


def __main__():
    s = 500
    n = 3
    k = 9

    print(1 / probability_kmer(s, n, k))


if __name__ == '__main__':
    __main__()
