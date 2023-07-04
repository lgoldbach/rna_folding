

def bp_to_dotbracket(bp: list, l: int) -> str:
    """Turn list of base-pairs into dot-bracket notation

    Args:
        bp (list): List of tuples defining base-pairs
        l (int): Length of sequence

    Returns:
        db (str): Secondary structure in dot-bracket notation

    """
    db = ["."] * l
    for pair in bp:
        db[pair[0] - 1] = "("
        db[pair[1] - 1] = ")"
    db = ''.join(db)
    return db