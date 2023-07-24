

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


def dotbracket_to_bp(db: str) -> set:
    """Convert RNA secondary structure from dot-bracket format to list of 
    base-pairs format

    Args:
        db (str): dot-bracket string

    Returns:
        bp (str): list of base-pair tuples (numeric).

    """
    opening_stack = []
    bp = []
    for i, site in enumerate(db):
        if site == "(":
            opening_stack.append(i)
        elif site == ")":
            bp.append((opening_stack.pop(), i))
    return set(bp)
