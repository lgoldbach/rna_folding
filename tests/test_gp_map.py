import pytest

from rna_folding.gp_map import GenotypePhenotypeMap


def example_gp():
    genotypes = ['AAA', 'AAU', 'AUA', 'AUU', 'UAA', 'UAU', 'UUA', 'UUU']
    phenotypes = [')..', ')..', ')((', ')..', ')((', ').(', ')..', ').(']
    alphabet = ["A", "U"]

    gpm = GenotypePhenotypeMap(genotypes=genotypes,
                         phenotypes=phenotypes,
                         alphabet=alphabet)
    gpm.add_hamming_edges()
    
    return gpm


def test_connected_components():
    # define the correct neutral components per phenotype
    ref_nc = {')..': [{'AAA', 'AUU', 'AAU'}, {'UUA'}], 
              ').(': [{'UAU', 'UUU'}], 
              ')((': [{'AUA'}, {'UAA'}]}
    
    gpm = example_gp()

    for ph in gpm.phenotypes:
        nc_ = gpm.neutral_components(phenotypes=[ph])[0]
        nc = []
        for c in nc_:
            nc.append(set(c))  # unpack generator objects

        for c in ref_nc[ph]:
            try:
                nc.remove(c)
            except ValueError:
                assert True == False, f"Neutral component {c} missing."

        assert len(nc) == 0, f"Neutral components {nc} are incorrect."
        

def test_robustness():
    # define correct robustness values
    ref_robustness = {tuple(sorted({'UUU', 'UAU'})): 0.333333333333333,
                      tuple(sorted({'AAU', 'AAA', 'AUU'})): 0.4444444444444444, 
                      tuple(sorted({'UUA'})): 0.0,
                      tuple(sorted({'AUA'})): 0.0,
                      tuple(sorted({'UAA'})): 0.0}
    
    gpm = example_gp()
    
    for ph in gpm.phenotypes:
        nc = gpm.neutral_components(phenotypes=[ph])[0]
        for component in nc:
            r = gpm.phenotype_robustness(component)
            assert r == pytest.approx(ref_robustness[tuple(sorted(component))])
