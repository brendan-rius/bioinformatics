from frequent_words import frequency_kmer


def find_clumps(sequence, L, t, k):
    """
    Find (L, t)-clumps of k-mers in a sequence
    :param sequence: the sequence
    :param L: the length of a window
    :param t: the minimum number of repetitions for a k-mer to be a clump in a window
    :param k: the size of the k-mer
    :return: a set of k-mers
    """
    len_sequence = len(sequence)
    if len_sequence < L:
        raise ValueError("the sequence has to be larger than the window")
    clumps = set()
    end = len_sequence - L + 1
    for i in range(0, len_sequence - L + 1):
        print("\r{}%".format(round(i * 100.0 / end, 3)), end='')
        frequencies = frequency_kmer(sequence[i:i + L], k)
        clumps.update({clump for clump, freq in frequencies.items() if freq >= t})
    return clumps


def find_clumps_2(sequence, L, t, k):
    """
    Find (L, t)-clumps of k-mers in a sequence
    :param sequence: the sequence
    :param L: the length of a window
    :param t: the minimum number of repetitions for a k-mer to be a clump in a window
    :param k: the size of the k-mer
    :return: a set of k-mers
    """
    frequencies = _frequency_kmer_stupid(sequence, k)
    res = set()
    print("start")
    for name, positions in frequencies.items():
        for idx1, pos1 in enumerate(positions):
            tab2 = positions[idx1 + t - 1:]
            if len(tab2) == 0:
                continue
            pos2 = tab2[0]
            if pos2 - pos1 + k <= L:
                res.add(name)
    return res


def _frequency_kmer_stupid(sequence, k):
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
    frequency_map = {}
    for i in range(0, max_index + 1):
        foundKmer = sequence[i:i + k]
        if foundKmer in frequency_map:
            frequency_map[foundKmer].append(i)
        else:
            frequency_map[foundKmer] = [i]
    return frequency_map


def __main__():
    file = open("ecoli_genome.txt", "r")
    sequence = file.read()
    file.close()
    t = 3
    k = 9
    L = 500

    clumps = find_clumps_2(sequence, L, t, k)
    print(len(clumps))
    print(list(clumps)[:10])


if __name__ == '__main__':
    __main__()
