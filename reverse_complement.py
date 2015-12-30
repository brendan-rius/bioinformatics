def reverse_complement(text):
    """
    Find the reverse complement of a DNA sequence
    :param text: the dna sequence
    :return: the reverse complement
    """
    complements = {"a": "t", "g": "c", "t": "a", "c": "g", "A": "T", "G": "C", "T": "A", "C": "G"}
    result = ""
    for c in text:
        result = complements[c] + result
    return result


def __main__():
    text = "GCTAGCT"
    print(reverse_complement(text))


if __name__ == '__main__':
    __main__()
