def to_profile_matrix(text):
    """
    Convert a profile matrix from coursera (string) to a usable profile matrix
    """
    matrix = [[float(n) for n in row.split()] for row in text.split('\n')]
    return list(zip(*matrix))  # We transpose the matrix


def to_sequences_list(text):
    return text.upper().split()
