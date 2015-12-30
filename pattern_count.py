from kmer import kmer


def pattern_count(text, pattern):
    """
    Count the number of times a pattern repeats in a string
    :param text: the string
    :param pattern: the pattern
    :return: the number of times
    """
    len_pattern = len(pattern)
    len_text = len(text)
    if len_text == 0 or len_text < len_pattern:
        return 0
    if len_pattern == 0:
        raise ValueError("The pattern cannot be empty")
    max_index = len(text) - len_pattern
    count = 0
    for i in range(0, max_index + 1):
        if kmer(text, position=i, k=len_pattern) == pattern:
            count += 1
    return count


def __main__():
    text = "CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC"
    pattern = "CGCG"
    print(pattern_count(text, pattern))


if __name__ == '__main__':
    __main__()
