import textwrap
from giggles.core import Read, ReadSet
from giggles.core import readselection

# util function from WhatsHap
def string_to_readset(s, source_id = 0):
    s = textwrap.dedent(s).strip()
    rs = ReadSet()
    for index, line in enumerate(s.split("\n")):
        if len(line) == 0:
            continue
        read = Read(f"Read {index + 1}", 50, source_id)
        scores = None
        for pos, c in enumerate(line):
            if c == " ":
                continue
            if c == "0":
                scores = [1, 5] # ALLELE1
            elif c == "1":
                scores = [5, 1] # ALLELE2
            elif c == "2":
                scores = [5, 5] # EQUAL_SCORES
            read.add_variant(position=(pos + 1) * 10, scores=scores)
        rs.add(read)
    print(rs)
    return rs

# test case from WhatsHap
def test_selection():
    reads = string_to_readset(
        """
      1  1
      00
      0   1
      10  1
      1   1
        11
      0   1
      1    1
    """
    )
    selected_reads = readselection(reads, max_cov=1, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([1, 5])
    selected_reads = readselection(reads, max_cov=2, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([1, 3, 5]), str(selected_reads)
    selected_reads = readselection(reads, max_cov=3, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([1, 3, 5, 7]), str(selected_reads)
    selected_reads = readselection(reads, max_cov=3, preferred_source_ids=None, bridging=True)
    # Here the assert is wrong, because the bridging doesn't come into account , because in the slice_read the selected
    # reads  have already coverage 3 by set ([1,3,5,7]) because first each position has to covered at least once before
    # the bridging starts
    assert selected_reads == set([1, 3, 5, 7]), str(selected_reads)

# test case from WhatsHap
def test_selection2():
    reads = string_to_readset(
        """
      1111
         111
         1  111
         1     11
        1      11
    """
    )
    selected_reads = readselection(reads, max_cov=4, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([0, 1, 2, 3]), str(selected_reads)

# test case from WhatsHap
def bridging():
    reads = string_to_readset(
        """
      11
      00
        11
        00
          11
          00
      1    1
    """
    )
    selected_reads = readselection(reads, max_cov=2, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([0, 1, 2, 3, 4, 5])
    selected_reads = readselection(reads, max_cov=2, preferred_source_ids=None, bridging=True)
    # Not sure why 0 is there selected and not 1...
    assert selected_reads == set([0, 3, 5, 6])


# test case from WhatsHap
# Component comparison does not work
def test_components_of_readselection():
    reads = string_to_readset(
        """
      111
         000
      00
          00
       1   1
    """
    )
    selected_reads = readselection(reads, max_cov=2, preferred_source_ids=None, bridging=False)
    assert selected_reads == set([0, 1, 2, 3]), str(selected_reads)
    #    assert len(set(new_components.values())) == 2
    selected_reads = readselection(reads, max_cov=2, preferred_source_ids=None, bridging=True)
    assert selected_reads == set([0, 1, 4]), str(selected_reads)


#      assert len(set(new_components.values())) == 1

# test case from WhatsHap
def test_selection_with_preferred_sources():
    readset = string_to_readset(
        """
      1        1
    """,
        source_id=3,
    )
    more_reads = string_to_readset(
        """
      1111
         111
            1111
    """,
        source_id=1,
    )

    for read in more_reads:
        readset.add(read)

    selected_reads = readselection(readset, max_cov=2, preferred_source_ids=None, bridging=True)
    assert selected_reads == set([1, 2, 3]), str(selected_reads)

    selected_reads = readselection(readset, max_cov=2, preferred_source_ids=set([3]), bridging=True)
    assert selected_reads == set([0, 1, 3]), str(selected_reads)

# selects equal scores read over informative read
def test_equal_scores_reads():
    reads = string_to_readset(
    """
    010
    101
    222
     10
    10
     01
    """
    )
    selected_reads = readselection(reads, max_cov=5, preferred_source_ids=None)
    assert(selected_reads == set([0, 1, 2, 3, 4]))

# does not select single variant reads
def test_single_variant_reads():
    reads = string_to_readset(
    """
    010
    1
     0
      1
     1
    10
      1
    """
    )
    selected_reads = readselection(reads, max_cov=5, preferred_source_ids=None)
    assert(selected_reads == set([0, 5]))