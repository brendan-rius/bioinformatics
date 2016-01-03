import itertools


def kmer(text, position, k):
    """
    Find the k-mer at a certain position in a text
    :param text: the text
    :param position: the position in the text
    :param k: the size of the k-mer
    :return: the k-mer
    """
    if position < 0:
        raise ValueError("Position has to be positive")
    if k < 0:
        raise ValueError("k has to be positive")

    len_text = len(text)
    if position + k > len_text:
        raise ValueError(
                "Text is not long enough to find {}-mer starting at position {} in a string of length {}".format(k,
                                                                                                                 position,
                                                                                                                 len_text))
    return text[position:position + k]


def all_kmers(k):
    """
    Generate all the possible k-mers
    :param k: the size of the k-mers
    """
    base = ['A', 'T', 'G', 'C']
    for p in itertools.product(base, repeat=k):
        yield ''.join(p)


def kmers(text, k):
    """
    Iterate over text k-mer after k-mer
    :param text: the text
    :param k: the size of the kmers
    :return:
    """

    if isinstance(text, str):
        stop = len(text) - k
        for i in range(0, stop + 1):
            yield text[i:i + k]
    else:
        for seq in text:
            stop = len(seq) - k
            for i in range(0, stop + 1):
                yield seq[i:i + k]
