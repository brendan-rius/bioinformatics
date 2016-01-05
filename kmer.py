import itertools


def kmer(text, position, k):
    """
    Find the k-mer at a certain position in a text
    :param text: the text
    :param position: the position in the text
    :param k: the size of the k-mer
    :return: the k-mer
    """
    return text[position:position + k]


def all_kmers(k):
    """
    Generate all the possible k-mers.

    Generate 4^k elements

    :param k: the size of the k-mers
    """
    base = ['A', 'T', 'G', 'C']
    for p in itertools.product(base, repeat=k):
        yield ''.join(p)


def kmers(text, k):
    """
    Generator to iterate over a text k-mer after k-mer.
    - CASE 1: the fist param is a string, iterates through all the possible k-mers in this string.
    Generates n-k+1 elements with n being the size of the text.
    - CASE 2: the second param if a list of strings: iterate through all the possible k-mers in this list of strings.
    Generates sum(n_i-k+1) elements with n_i being the size of the i-th text.

    :param text: the text
    :param k: the size of the k-mers
    """

    if isinstance(text, str):  # If the text is a string
        stop = len(text) - k
        for i in range(0, stop + 1):
            yield text[i:i + k]
    else:
        for seq in text:
            stop = len(seq) - k
            for i in range(0, stop + 1):
                yield seq[i:i + k]
