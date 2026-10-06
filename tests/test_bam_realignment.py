import pytest
from collections import namedtuple

# Import your classes here (adjust the import path based on your package structure)
from giggles.variants import BAMReader

# --- MOCKS ---

class MockAligner:
    """A mock aligner that records the state it was called with and returns a dummy distance."""
    def __init__(self):
        self.called_state = None
        self.called_query = None
        self.called_allele = None

    def get_distance(self, query: str, allele: str, state: int):
        self.called_state = state
        self.called_query = query
        self.called_allele = allele
        return 0  # Dummy perfect score

class MockVariant:
    def __init__(self, position_on_ref, ref_allele, alt_allele):
        self.position_on_ref = position_on_ref
        self.reference_allele = ref_allele
        self.alternative_allele = [alt_allele]

class MockBamRead:
    def __init__(self, query_sequence, reference_start):
        self.query_sequence = query_sequence
        self.reference_start = reference_start

# --- CONSTANTS FOR TESTS ---

# 100 bp dummy reference chromosome
CHROM_SEQ = "A" * 30 + "G" * 40 + "C" * 30 
VARIANT_POS = 40
REF_ALLELE = "G" * 20  # Variant is 20bp long (pos 40 to 60)
ALT_ALLELE = "T" * 20
OVERHANG = 10

# CIGAR Operators
M, I, D, N, S, H = 0, 1, 2, 3, 4, 5


@pytest.fixture
def mock_aligner():
    return MockAligner()

@pytest.fixture
def sv_variant():
    return MockVariant(VARIANT_POS, REF_ALLELE, ALT_ALLELE)


# --- TEST CASES ---

def test_state_0_global_alignment(mock_aligner, sv_variant):
    """
    Test State 0 (Full Coverage).
    Read1 spans 20 to 70 (Covers the 20bp variant with 10bp overhang on both sides).
    Read2 spans 10 to 80 (Covers the 20bp variant with 10bp overhang on both sides) but with an additional 10bp.
    Read3 spans 10 to 80 (Covers the 20bp variant with 10bp overhang on both sides) but with a 10bp deletion.
    Read4 spans 35 to 65 (Covers the 20bp variant but not enough for the overhang)
    Read5 spans 39 to 70 with a soft-clip sequence
    """
    # 40 Match (covers 30 to 70 on ref)
    cigartuples = [(M, 50)]
    read = MockBamRead("A" * 10 + "G" * 40, reference_start=20)
    
    i = 0
    consumed = 20
    query_pos = 20 
    
    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )

    assert mock_aligner.called_state == 0, f"[Read1] Expected State 0, got {mock_aligner.called_state}"

    # 30 Match, 20 Mismatch, 5 Insertion and 20 Match
    cigartuples = [(M, 50), (I, 5), (M, 20)]
    read = MockBamRead("A" * 20 + "G" * 10 + "T" * 25 + "G" * 10 + "C" * 10, reference_start=10)
    
    i = 0
    consumed = 30
    query_pos = 30 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )

    assert mock_aligner.called_state == 0, f"[Read2] Expected State 0, got {mock_aligner.called_state}"

    # 30 Match, 5 Deletion and 20 Match
    cigartuples = [(M, 30), (D, 5), (M, 20)]
    read = MockBamRead("A" * 20 + "G" * 35, reference_start=10)
    
    i = 0
    consumed = 30
    query_pos = 30 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )

    assert mock_aligner.called_state == 0, f"[Read3] Expected State 0, got {mock_aligner.called_state}"

    # testing for when overhangs are not availble
    cigartuples = [(M, 30)]
    read = MockBamRead("G" * 30, reference_start=35)
        
    i = 0
    consumed = 5
    query_pos = 5 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )

    assert mock_aligner.called_state == 0, f"[Read4] Expected State 0, got {mock_aligner.called_state}"

    # 4 Soft Clip, 31 Match (covers 39 to 70 on ref)
    cigartuples = [(S, 4), (M, 31)]
    read = MockBamRead("X" * 4 + "G" * 31, reference_start=39)
    
    # Read starts exactly at 39
    i = 1
    consumed = 1
    query_pos = 5

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )

    assert mock_aligner.called_state == 0, f"[Read5] Expected State 0, got {mock_aligner.called_state}"
    
        

def test_state_1_partial_start(mock_aligner, sv_variant):
    """
    Test State 1 (Partial Start).
    Read1 starts at 40 (Variant start) with a soft clip, and extends past the variant.
    Read2 starts at 40 (Variant start) perfectly with NO soft clips, and extends past the variant.
    Read3 starts at 40 (Variant start), has an insertion in the middle, and extends past.
    Read4 starts at 40 (Variant start), has a deletion in the middle, and extends past. 
    """
    
    # 15 Soft Clip, 30 Match (covers 40 to 70 on ref)
    cigartuples = [(S, 5), (M, 30)]
    read = MockBamRead("X" * 5 + "G" * 30, reference_start=40)
    
    i = 1
    consumed = 0
    query_pos = 15 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 1, f"[Read1] Expected State 1, got {mock_aligner.called_state}"

    # 30 Match (covers 40 to 70 on ref)
    cigartuples = [(M, 30)]
    read = MockBamRead("G" * 30, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 1, f"[Read2] Expected State 1, got {mock_aligner.called_state}"

    # ---------------------------------------------------------
    # Read 3: Insertion inside the variant
    # ---------------------------------------------------------
    # 10 Match, 5 Insertion, 20 Match (covers 40 to 70 on ref)
    cigartuples = [(M, 10), (I, 5), (M, 20)]
    read = MockBamRead("G" * 10 + "T" * 5 + "G" * 20, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 1, f"[Read3] Expected State 1, got {mock_aligner.called_state}"

    # 10 Match, 5 Deletion, 20 Match (covers 40 to 80 on ref)
    # Total read length = 10 (M) + 0 (D) + 20 (M) = 30 bases
    cigartuples = [(M, 10), (D, 5), (M, 20)]
    read = MockBamRead("G" * 20 + "C" * 10, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 1, f"[Read4] Expected State 1, got {mock_aligner.called_state}"


def test_state_2_partial_end(mock_aligner, sv_variant):
    """
    Test State 2 (Partial End).
    Read1 starts at 20 (before variant) and cuts off at 55 (midway through variant).
    Read2 starts at 20 (before variant) and ends exactly at 60 (exact boundary, 0 right overhang).
    Read3 starts at 20 (before variant), has an insertion inside, and cuts off at 55.
    Read4 starts at 20 (before variant), has a deletion inside, and cuts off at 55.
    """

    # 35 Match (covers 20 to 55 on ref)
    cigartuples = [(M, 35)]
    read = MockBamRead("A" * 10 + "G" * 25, reference_start=20)
    
    i = 0
    consumed = 20
    query_pos = 20 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 2, f"[Read1] Expected State 2, got {mock_aligner.called_state}"

    # 40 Match (covers 20 to 60 on ref)
    cigartuples = [(M, 40)]
    read = MockBamRead("A" * 10 + "G" * 30, reference_start=20)
    
    i = 0
    consumed = 20
    query_pos = 20 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 2, f"[Read2] Expected State 2, got {mock_aligner.called_state}"

    # 25 Match, 10 Insertion, 10 Match (covers 20 to 55 on ref)
    cigartuples = [(M, 25), (I, 5), (M, 10)]
    read = MockBamRead("A" * 10 + "G" * 15 + "T" * 5 + "G" * 10, reference_start=20)
    
    i = 0
    consumed = 20
    query_pos = 20 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 2, f"[Read3] Expected State 2, got {mock_aligner.called_state}"

    # 25 Match, 5 Deletion, 5 Match (covers 20 to 55 on ref)
    cigartuples = [(M, 25), (D, 5), (M, 5)]
    read = MockBamRead("A" * 10 + "G" * 20, reference_start=20)
    
    i = 0
    consumed = 20
    query_pos = 20 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 2, f"[Read4] Expected State 2, got {mock_aligner.called_state}"


def test_state_3_contained(mock_aligner, sv_variant):
    """
    Test State 3 (Contained).
    Read1 starts at 40 (Variant start) and cuts off at 55 (midway through variant).
    Read2 starts at 40 and ends exactly at 60, surrounded by soft clips on both sides.
    Read3 starts at 40, has an insertion inside, and ends at 55 (midway through).
    Read4 starts at 40, has a deletion inside, and ends exactly at 60.
    """

    # 15 Match (covers 40 to 55 on ref)
    cigartuples = [(M, 15)]
    read = MockBamRead("G" * 15, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 3, f"[Read1] Expected State 3, got {mock_aligner.called_state}"

    # 1 Soft Clip, 20 Match, 1 Soft Clip (covers 40 to 60 on ref)
    cigartuples = [(S, 1), (M, 20), (S, 1)]
    read = MockBamRead("X" * 1 + "G" * 20 + "Y" * 1, reference_start=40)
    
    i = 1
    consumed = 0
    query_pos = 1

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 3, f"[Read2] Expected State 3, got {mock_aligner.called_state}"

    # 5 Match, 5 Insertion, 10 Match (covers 40 to 55 on ref)
    cigartuples = [(M, 5), (I, 5), (M, 10)]
    read = MockBamRead("G" * 5 + "T" * 5 + "G" * 10, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 3, f"[Read3] Expected State 3, got {mock_aligner.called_state}"

    # 10 Match, 5 Deletion, 5 Match (covers 40 to 60 on ref)
    cigartuples = [(M, 10), (D, 5), (M, 5)]
    read = MockBamRead("G" * 15, reference_start=40)
    
    i = 0
    consumed = 0
    query_pos = 0 

    mock_aligner.called_state = None
    _ = BAMReader.realign(
        mock_aligner, sv_variant, read, cigartuples, i, consumed, query_pos, CHROM_SEQ, OVERHANG
    )
    assert mock_aligner.called_state == 3, f"[Read4] Expected State 3, got {mock_aligner.called_state}"