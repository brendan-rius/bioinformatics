from hamming_distance import hamming_distance


def immediate_neighbors(kmer):
    nucleotides = {'A', 'T', 'G', 'C'}
    neighbors = set()
    for i, c in enumerate(kmer):
        other_nucleotides = nucleotides.copy()
        other_nucleotides.remove(c)
        for n in other_nucleotides:
            neighbors.add(kmer[:i] + n + kmer[i + 1:])
    return neighbors


def neighbors(kmer, distance):
    """
    Return a set of all the neighbors of a k-mer
    :param kmer: the k-mer
    :param distance: the maximum distance between two k-mers for them to be neighbors
    :return: the set of all neighbors
    """
    k = len(kmer)
    if k == 0:
        return {}
    if distance == 0:
        return {kmer}
    if k == 1:
        return {'A', 'T', 'C', 'G'}
    suffix = kmer[1:]
    suffix_neighbors = neighbors(suffix, distance)
    result = set()
    for suffix_neighbor in suffix_neighbors:
        if hamming_distance(suffix_neighbor, suffix) == distance:
            result.add(kmer[0] + suffix_neighbor)
        else:
            result.add('A' + suffix_neighbor)
            result.add('T' + suffix_neighbor)
            result.add('C' + suffix_neighbor)
            result.add('G' + suffix_neighbor)
    return result


def __main__():
    kmer = "ACGT"
    distance = 3
    print('\n'.join(neighbors(kmer, distance)))


if __name__ == '__main__':
    __main__()
