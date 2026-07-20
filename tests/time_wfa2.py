import pytest
from giggles.ext import WFAWrapper
from giggles.align import edit_distance as ed
from giggles.timer import StageTimer

@pytest.mark.skip(reason="need to time and not test.")
def time_wfa2():

    def read_alleles():
        alleles = []
        with open('tests/data/other/time-testing-wfa2/large-alleles.tsv', 'r') as file:
            for line in file:
                line = line.strip().split('\t')
                assert len(line) == 2 
                alleles.append(line[0])
                for alt in line[1].split(','):
                    alleles.append(alt)
        return alleles

    def read_queries():
        queries = []
        with open('tests/data/other/time-testing-wfa2/queries.tsv', 'r') as file:
            for line in file:
                line = line.strip().split('\t')
                assert len(line) == 2
                state = int(line[0])
                query = line[1]
                queries.append((state, query))
        return queries
    
    timer = StageTimer()
    alleles = read_alleles()
    queries = read_queries()

    # timing WFA2-lib
    count = 0
    times = {0: [], 2: []}
    for state, query in queries:
        q_len = len(query)
        print(f"Aligning #{count} -> State {state}")
        aligner = WFAWrapper()
        with timer('align'):
            if state == 0:
                for idx, allele in enumerate(alleles):
                    a_len = len(allele)
                    # checking for heuristic 1
                    if a_len >= 1.2*q_len or a_len <= q_len/1.2:
                        #print(f'\tScore for allele {idx} (allele={len(allele)};query={len(query)}): NaN (Heuristic 1)')
                        continue
                    score = aligner.align(allele=allele, query=query, state=state)
                    #score = ed(allele, query)
                    print(f'\tScore for allele {idx} (allele={len(allele)};query={len(query)}): {score}')
            if state == 2:
                for idx, allele in enumerate(alleles):
                    a_len = len(allele)
                    # checking for heuristic 1
                    if a_len <= q_len/1.2:
                        #print(f'\tScore for allele {idx} (allele={len(allele)};query={len(query)}): NaN (Heuristic 2)')
                        continue
                    score = aligner.align(allele=allele, query=query, state=state)
                    #score = ed(allele, query)
                    print(f'\tScore for allele {idx} (allele={len(allele)};query={len(query)}): {score}')
        times[state].append(timer.elapsed('align'))
        timer._elapsed['align'] = 0
        count += 1
    
    print(times)


if __name__=='__main__':
    time_wfa2()