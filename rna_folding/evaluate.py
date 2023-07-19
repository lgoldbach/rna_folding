from rna_folding.utils import dotbracket_to_bp


def f1_score(A: str, B: str, abstract_level: int = 0) -> float:
    """Compute F1 score for two RNA secondary structures

    Args:
        A (str): First sequence
        B (str): Second sequence
        abstract_level (int, optional): Create sequences on abstract shape 
            level [0-5]. 0 is full dot-bracket notation, 5 is maximum 
            abstraction. Defaults to 0.

    Returns:
        bool: True if perfect match, False otherwise

    """
    A_ = dotbracket_to_bp(A)
    B_ = dotbracket_to_bp(B)

    tp = A_.intersection(B_)
    fp = B_.difference(tp)
    fn = A_.difference(tp)

    sensitivity = len(tp) / (len(tp) + len(fn))
    precision =  len(tp) / (len(tp) + len(fp))

    if sensitivity and precision:  # avoid division by zero
        f1 = (2*sensitivity*precision) / (sensitivity+precision)
    else:
        f1 = 0

    return f1


def compare_db(A: str, B: str, abstract_level: int = 0) -> bool:
    """Compare two RNA secondary structrues in dot-bracket notation

    Args:
        A (str): First sequence
        B (str): Second sequence
        abstract_level (int, optional): Create sequences on abstract shape 
            level [0-5]. 0 is full dot-bracket notation, 5 is maximum 
            abstraction. Defaults to 0.

    Returns:
        bool: True if perfect match, False otherwise

    """
    ############
    ### Placeholder for RNA shape abstraction code
    ############

    if A == B:
        return True
    else:
        False

