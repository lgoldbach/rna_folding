import argparse
import numpy as np
from itertools import product


def gp_map_sample(l: int, a: str, s: int):
    """Sample from a g-p map. For efficiency this creates samples with replacement. It is meant for large g-p maps.
    
    Args:
        l (int): Sequence length.
        a (str): Alphabet as a continuous string, e.g. "AUGC".
        s (int): Sample size

    returns:
        gp_map_sample (Generator)

    """
    counter = 0
    while counter < s:
        yield(random_sequence(l, a))
        counter += 1


def random_sequence(l, a):
    """Create a random sequence of length l from alphabet a

    Args:
        l (int): Sequence length.
        a (str): Alphabet as a continuous string, e.g. "AUGC".

    Returns:
        seq (str): A sequence (str).

    """
    seq_ = []
    for i in range(l):
        letter_numeric = np.random.choice(len(a), size=1)[0]
        letter = a[letter_numeric]
        seq_.append(letter)
    seq = "".join(seq_)
    return seq


def random_sample_generator(generator, generator_len, sample_size):
    """Create a random sample from a generator of known size (fails for very large generators due to numpy.random.choice())

    Args:
        generator (generator): Generator object
        generator_len (int): Size of the generator
        sample_size (int): Sample szie.

    Returns:
        sample (list): List of sampled objects
    """
    if sample_size > generator_len:
        raise ValueError("sample_size can't be larger than generator_len")
    
    rand_idx = np.sort(np.random.choice(generator_len, size=sample_size, replace=False))
    sample = []
    count = 0
    for idx in rand_idx:
        while count != idx:
            next(generator)
            count += 1
        sample.append(next(generator))
        count += 1

    return sample


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-l", "--seq_len", help="Sequence length", type=int)
    parser.add_argument("-a", "--alphabet", help="Define the alphabet as an unseparated string, e.g: AUGC", type=str)
    parser.add_argument("-s", "--sample_size", help="Number of samples to return", type=int)

    args = parser.parse_args()

    with open(args.output, "w") as f:
        for g in gp_map_sample(l=args.seq_len, a=args.alphabet, s=args.sample_size):
            f.write(g + "\n")
        f.close()
