'''
Testing the functions for distance score calculations
'''

# test cases taken from WhatsHap (version 2.8)

from giggles.align import edit_distance as ed
from giggles.ext import WFAWrapper


STRING_PAIRS = [
    ("", ""),
    ("", "A"),
    ("A", "A"),
    ("AB", ""),
    ("AB", "ABC"),
    ("TGAATCCC", "CCTGAATC"),
    ("ANANAS", "BANANA"),
    ("SISSI", "MISSISSIPPI"),
    ("GGAATCCC", "TGAGGGATAAATATTTAGAATTTAGTAGTAGTGTT"),
    ("TCTGTTCCCTCCCTGTCTCA", "TTTTAGGAAATACGCC"),
    (
        "TGAGACACGCAACATGGGAAAGGCAAGGCACACAGGGGATAGG",
        "AATTTATTTTATTGTGATTTTTTGGAGGTTTGGAAGCCACTAAGCTATACTGAGACACGCAACAGGGGAAAGGCAAGGCACA",
    ),
    (
        "TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "TTTTAGGAAATACGCCTGGTGGGGTTTGGAGTATAGTGAAAGATAGGTGAGTTGGTCGGGTG",
    ),
    ("A", "TCTGCTCCTGGCCCATGATCGTATAACTTTCAAATTT"),
    ("GCGCGGACT", "TAAATCCTGG"),
]

def test_edit_distance():
    assert ed("", "") == 0
    assert ed("", "A") == 1
    assert ed("A", "B") == 1
    assert ed("A", "A") == 0
    assert ed("A", "AB") == 1
    assert ed("BA", "AB") == 2
    for s, t in STRING_PAIRS:
        assert ed(s, "") == len(s)
        assert ed("", s) == len(s)
        assert ed(s, t) == ed(t, s)


def test_edit_distance_bytes():
    assert ed(b"", b"") == 0
    assert ed(b"", b"A") == 1
    assert ed(b"A", b"B") == 1
    assert ed(b"A", b"A") == 0
    assert ed(b"A", b"AB") == 1
    assert ed(b"BA", b"AB") == 2
    for s, t in STRING_PAIRS:
        s = s.encode("ascii")
        t = t.encode("ascii")
        assert ed(s, b"") == len(s)
        assert ed(b"", s) == len(s)
        assert ed(s, t) == ed(t, s)


def test_edit_distance_banded():

    def assert_banded(s, t, maxdiff):
        banded_dist = ed(s, t, maxdiff=maxdiff)
        true_dist = ed(s, t)
        if true_dist > maxdiff:
            assert banded_dist > maxdiff
        else:
            assert banded_dist == true_dist
    
    for maxdiff in range(5):
        assert_banded("ABC", "", maxdiff)

        for s, t in STRING_PAIRS:
            assert_banded(s, "", maxdiff)
            assert_banded("", s, maxdiff)
            assert_banded(s, t, maxdiff)
            assert_banded(t, s, maxdiff)

def test_wfa():

    def expected_error(allele, query, state):
        try:
            _ = aligner.align(allele=allele, query=query, state=state)
            assert False
        except RuntimeError:
            # expected
            pass

    
    aligner = WFAWrapper()

    '''
    Most of the alignments give distance 0 because the bandwidth is 0.
    '''
    # Testing alignment type 0
    assert(aligner.align(allele="AATGC", query="AATGC", state=0) == 0)
    expected_error(allele="AATGC", query="ATGC", state=0)
    assert(aligner.align(allele="AATGC", query="ATTGC", state=0) == 0)
    expected_error(allele="AATGC", query="ATC", state=0)
    assert(aligner.align(allele="AATGC", query="AGTCC", state=0) == 0)
    expected_error(allele="AAAAA", query="T", state=0)
    assert(aligner.align(allele="AAAAA", query="TTTTT", state=0) == 0)
    expected_error(allele="AAAAAA", query="TTTT", state=0)

    # Testing alignment type 2
    assert(aligner.align(allele="AAAAA", query="T", state=2) == 0)
    assert(aligner.align(allele="TAAAA", query="T", state=2) == 0)
    assert(aligner.align(allele="GCGAA", query="GCA", state=2) == 0)
    assert(aligner.align(allele="AAAAGCG", query="GCG", state=2) == 0)
    assert(aligner.align(allele="AATGC", query="AT", state=2) == 0)
    assert(aligner.align(allele="AATGC", query="TG", state=2) == 0)
    assert(aligner.align(allele="AATGC", query="ATG", state=2) == 0)

    # Testing alignment type 1
    assert(aligner.align(allele="AAAAA", query="T", state=1) == 0)
    assert(aligner.align(allele="AAAAT", query="T", state=1) == 0)
    assert(aligner.align(allele="AAGCG", query="GCA", state=1) == 0)
    assert(aligner.align(allele="GCGAAAA", query="GCG", state=1) == 0)
    assert(aligner.align(allele="AATGC", query="TG", state=1) == 0) # Because of the bandwidth, the alignment happens between GC and TG which has an ed of 2
    assert(aligner.align(allele="AATGC", query="AT", state=1) == 0) # Same bandwidth artefact
    assert(aligner.align(allele="AATGC", query="ATG", state=1) == 0)    # Same bandwidth artefact
    assert(aligner.align(allele="AAAATGC", query="AAATG", state=1) == 0)    # Same bandwidth artefact
    assert(aligner.align(allele="AAAAATGC", query="AAAATG", state=1) == 1)  # Should have bandwidth of at least 1 which allows AAAATG to align with AAAATGC
    
    # Testing alignment type 3
    assert(aligner.align(allele="AAAAA", query="T", state=3) == 0)
    assert(aligner.align(allele="AATAA", query="T", state=3) == 0)
    assert(aligner.align(allele="TTGCATT", query="GCA", state=3) == 0)
    assert(aligner.align(allele="TTGACAGTT", query="GCG", state=3) == 0)


    # Testing alignment type 0
    assert(aligner.align(allele="AAATGC", query="AAATGC", state=0) == 0)
    expected_error(allele="AAATGC", query="AATGC", state=0)
    assert(aligner.align(allele="AAATGC", query="AATTGC", state=0) == 1)
    expected_error(allele="AAATGC", query="AATC", state=0)
    assert(aligner.align(allele="AAATGC", query="AAGTCC", state=0) == 1)
    expected_error(allele="AAAAAA", query="T", state=0)
    assert(aligner.align(allele="AAAAAA", query="TTTTTT", state=0) == 1)
    expected_error(allele="AAAAAAA", query="TTTTT", state=0)
 
    # Testing alignment type 3
    assert(aligner.align(allele="TTGACAGTT", query="TGACAGT", state=3) == 0)
    assert(aligner.align(allele="TTGACAGTT", query="TGATAGT", state=3) == 1)
    assert(aligner.align(allele="TTGACAGTT", query="TGATTGT", state=3) == 1)
    
    
    '''
    # Testing alignment type 0
    print(f"Type 0 (AATGC vs AATGC) -> Expected: 0, Got: {aligner.align(allele='AATGC', query='AATGC', state=0)}")
    print(f"Type 0 (AATGC vs ATGC)  -> Expected: 1, Got: {aligner.align(allele='AATGC', query='ATGC', state=0)}")
    print(f"Type 0 (AATGC vs ATTGC) -> Expected: 1, Got: {aligner.align(allele='AATGC', query='ATTGC', state=0)}")
    print(f"Type 0 (AATGC vs ATC)   -> Expected: 2, Got: {aligner.align(allele='AATGC', query='ATC', state=0)}")
    print(f"Type 0 (AATGC vs AGTCC) -> Expected: 2, Got: {aligner.align(allele='AATGC', query='AGTCC', state=0)}")
    print(f"Type 0 (AAAAA vs T)     -> Expected: 5, Got: {aligner.align(allele='AAAAA', query='T', state=0)}")

    # Testing alignment type 2
    print(f"Type 2 (AAAAA vs T)     -> Expected: 1, Got: {aligner.align(allele='AAAAA', query='T', state=2)}")
    print(f"Type 2 (TAAAA vs T)     -> Expected: 0, Got: {aligner.align(allele='TAAAA', query='T', state=2)}")
    print(f"Type 2 (GCGAA vs GCA)   -> Expected: 1, Got: {aligner.align(allele='GCGAA', query='GCA', state=2)}")
    print(f"Type 2 (AAAAGCG vs GCG) -> Expected: 3, Got: {aligner.align(allele='AAAAGCG', query='GCG', state=2)}")

    # Testing alignment type 1
    print(f"Type 1 (AAAAA vs T)     -> Expected: 1, Got: {aligner.align(allele='AAAAA', query='T', state=1)}")
    print(f"Type 1 (AAAAT vs T)     -> Expected: 0, Got: {aligner.align(allele='AAAAT', query='T', state=1)}")
    print(f"Type 1 (AAGCG vs GCA)   -> Expected: 1, Got: {aligner.align(allele='AAGCG', query='GCA', state=1)}")
    print(f"Type 1 (GCGAAAA vs GCG) -> Expected: 3, Got: {aligner.align(allele='GCGAAAA', query='GCG', state=1)}")

    # Testing alignment type 0 (additional)
    print(f"Type 0 (AATGC vs AATGC) -> Expected: 0, Got: {aligner.align(allele='AATGC', query='AATGC', state=0)}")
    print(f"Type 0 (AAAAA vs TTTTT) -> Expected: 5, Got: {aligner.align(allele='AAAAA', query='TTTTT', state=0)}")
    print(f"Type 0 (AAAAAA vs TTTT) -> Expected: 6, Got: {aligner.align(allele='AAAAAA', query='TTTT', state=0)}")

    # Testing alignment type 2 (additional)
    print(f"Type 2 (AATGC vs AT)    -> Expected: 1, Got: {aligner.align(allele='AATGC', query='AT', state=2)}")
    print(f"Type 2 (AATGC vs TG)    -> Expected: 2, Got: {aligner.align(allele='AATGC', query='TG', state=2)}")
    print(f"Type 2 (AATGC vs ATG)   -> Expected: 1, Got: {aligner.align(allele='AATGC', query='ATG', state=2)}")

    # Testing alignment type 1 (additional)
    print(f"Type 1 (AATGC vs TG)    -> Expected: 1, Got: {aligner.align(allele='AATGC', query='TG', state=1)}")
    print(f"Type 1 (AATGC vs AT)    -> Expected: 2, Got: {aligner.align(allele='AATGC', query='AT', state=1)}")
    print(f"Type 1 (AATGC vs ATG)   -> Expected: 1, Got: {aligner.align(allele='AATGC', query='ATG', state=1)}")

    # Testing alignment type 3
    print(f"Type 3 (AAAAA vs T)       -> Expected: 1, Got: {aligner.align(allele='AAAAA', query='T', state=3)}")
    print(f"Type 3 (AATAA vs T)       -> Expected: 0, Got: {aligner.align(allele='AATAA', query='T', state=3)}")
    print(f"Type 3 (TTGCATT vs GCA)   -> Expected: 0, Got: {aligner.align(allele='TTGCATT', query='GCA', state=3)}")
    print(f"Type 3 (TTGACAGTT vs GCG) -> Expected: 2, Got: {aligner.align(allele='TTGACAGTT', query='GCG', state=3)}")
    '''
    allele='ATATTAAGCATGAATAAACATTAGATACTATTAAAATCCTATATATTAACAAAGCCAAAAGTTTCAAACTTTACTTTTTCCCAACATTCTTGTGAAATATGACACATCCCAATCTTAACAGATGCTCATTTGGGATACTGTACTTGTGAGTGGAAGTGTGTATATTTGTGTGCAAGTGTGTACTCATATACTTCCACCTTACCACCCTAGAAAGGCATGATGAAAATTTAAGATAGAAGGAAAATATAAATTGAAAAAAAAAAACCTTAACAAATGATTCTGACAAATATCTTCTCTTCCAGGGAGAGTCACTGAGCCAGAATAAAATTGAACACTAAATATTCTAAGAAAAAAAGGAATCTAGTTTGTCAAAATGTGACTTGAATTAATAGATAAGGAGAGTCAGATGATAAGAGGGTCAAAATTATGTTTATCTTAGGAAAAGTAGAATAGAAAATTTATAAGCAGATTAAAAACACATAATAAAAGTAGTAAATAATAATGACAGTATCTCAAATCAGTGCAG'
    query='ATATTAAGCACTTTACTTTTTCCCAACATTCTTGTGAAATATGACACATCCCAATCTTAACAGATGCTCATTTGGGATACTGTACTTGTGAGTGGAAGTGTGTATATTTGTGTGCAAGTGTGTACTCATATACTTCCACCTTACCACCCTAGAAAGGCATGATGAAAATTTAAGATAGAAGGAAAATATAAATTGAAAAAAAAAAACCTTAACAAATGATTCTGACAAATATCTTCTCTTTCCAGGGAGAATCACTGAGCCAGAATAAAATTGAACACTAAATATTCTAAGAAAAAAGGAATCTAGTTTGTCAAAATGTGACTTGAATTAATAGATAAGGAGAGTCAGATGATAAGAGGGTCAAAATTATGTTTATCTTAGGAAAAGTAGAATAGAAAATTTATAAGCAGATTAAAAACACATAATAAAAGTAGTAAATAATAATGACAGTATCTCAAATCAGTGCAG'
    for i in range(10):
        assert aligner.align(allele, query, state=0) == 61
    