import pytest

from rna_folding.utils import dotbracket_to_genotype

def test_dotbracket_to_genotype():
    dbs = ["(((...)))", ".....", "((()))"]
    genotype_gc = ["GGGGGGCCC", "GGGGG", "GGGCCC"]
    genotype_au = ["AAAAAAUUU", "AAAAA", "AAAUUU"]
    genotype_gc_rand = ["GCGCCGCGC", "GCGCC", "GCGCGC"]
    seed = 1996

    for db, gt_ref in zip(dbs, genotype_gc):
        gt = dotbracket_to_genotype(db, base_pair="GC", random=False)
        assert gt == gt_ref

    for db, gt_ref in zip(dbs, genotype_au):
        gt = dotbracket_to_genotype(db, base_pair="AU", random=False)
        assert gt == gt_ref

    for db, gt_ref in zip(dbs, genotype_gc_rand):
        gt = dotbracket_to_genotype(db, base_pair="GC", random=True, 
                                    seed=seed)
        assert gt == gt_ref
