def _hash_nucleotide(nucleotide):
    return "acgt".find(nucleotide.lower())


def _unhash_nucleotide(number):
    return "acgt"[number]


def hash_kmer(kmer):
    """
    Hash a kmer as a number. The size of the k-mer will be lost while hashing, so you will need to store it to
    unhash it after.
    :param kmer: the kmer
    :return: a unique number
    """
    len_kmer = len(kmer)
    if len_kmer == 0:
        raise ValueError("kmer cannot be null")
    elif len_kmer == 1:
        return "acgt".find(kmer[0].lower())
    return 4 * hash_kmer(kmer[:-1]) + hash_kmer(kmer[-1])


def unhash_kmer(hash, size):
    """
    Unhash a kmer of a certain size
    :param hash: the hashed kmer
    :param size: the size of the target kmer
    :return:
    """
    kmer = ""
    while hash != 0:
        kmer = _unhash_nucleotide(hash % 4) + kmer
        hash //= 4
    len_kmer = len(kmer)
    kmer = "a" * (size - len_kmer) + kmer  # fill with "a" until the size is correct
    return kmer


def __main__():
    kmer = "TCTGATACTCTGTGGG"
    hash = hash_kmer(kmer)  # Hash
    print(unhash_kmer(hash, len(kmer)).upper())  # Unhash
