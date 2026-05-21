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

    def expected_error(text, pattern, state):
        try:
            _ = aligner.align(text=text, pattern=pattern, state=state)
            assert False
        except RuntimeError:
            # expected
            pass

    
    ### Testing with large bandwidth.
    aligner = WFAWrapper(bandwidth=10)
    
    # Testing alignment type 0
    assert(aligner.align(text="AATGC", pattern="AATGC", state=0) == 0)
    assert(aligner.align(text="AATGC", pattern="ATGC", state=0) == 1)
    assert(aligner.align(text="AATGC", pattern="ATTGC", state=0) == 1)
    assert(aligner.align(text="AATGC", pattern="ATC", state=0) == 2)
    assert(aligner.align(text="AATGC", pattern="AGTCC", state=0) == 2)
    assert(aligner.align(text="AAAAA", pattern="T", state=0) == 5)

    # Testing alignment type 1
    assert(aligner.align(text="AAAAA", pattern="T", state=1) == 1)
    assert(aligner.align(text="TAAAA", pattern="T", state=1) == 0)
    assert(aligner.align(text="GCGAA", pattern="GCA", state=1) == 1)
    assert(aligner.align(text="AAAAGCG", pattern="GCG", state=1) == 3)

    # Testing alignment type 2
    assert(aligner.align(text="AAAAA", pattern="T", state=2) == 1)
    assert(aligner.align(text="AAAAT", pattern="T", state=2) == 0)
    assert(aligner.align(text="AAGCG", pattern="GCA", state=2) == 1)
    assert(aligner.align(text="GCGAAAA", pattern="GCG", state=2) == 3)

    ### Testing with smaller bandwidth
    aligner = WFAWrapper(bandwidth=2)
    
    # Testing alignment type 0
    assert(aligner.align(text="AATGC", pattern="AATGC", state=0) == 0)
    assert(aligner.align(text="AAAAA", pattern="TTTTT", state=0) == 5)
    expected_error(text="AAAAAA", pattern="T", state=0)
    expected_error(text="AAAAAA", pattern="TT", state=0)
    expected_error(text="AAAAAA", pattern="TTT", state=0)
    assert(aligner.align(text="AAAAAA", pattern="TTTT", state=0) == 6)

    # Testing alignment type 1
    assert(aligner.align(text="AATGC", pattern="AT", state=1) == 1)
    assert(aligner.align(text="AATGC", pattern="TG", state=1) == 2)
    assert(aligner.align(text="AATGC", pattern="ATG", state=1) == 1)
    
    # Testing alignment type 2
    assert(aligner.align(text="AATGC", pattern="TG", state=2) == 1)
    assert(aligner.align(text="AATGC", pattern="AT", state=2) == 2)
    assert(aligner.align(text="AATGC", pattern="ATG", state=2) == 1)
    

    # Testing alignment type 3 (does not care about bandwidth)
    assert(aligner.align(text="AAAAA", pattern="T", state=3) == 1)
    assert(aligner.align(text="AATAA", pattern="T", state=3) == 0)
    assert(aligner.align(text="TTGCATT", pattern="GCA", state=3) == 0)
    assert(aligner.align(text="TTGACAGTT", pattern="GCG", state=3) == 2)
