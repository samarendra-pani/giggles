"""
Pedigree-related functions
"""
from abc import ABC, abstractmethod

class RecombinationCostComputer(ABC):
    @abstractmethod
    def compute(self, positions):
        pass

class UniformRecombinationCostComputer(RecombinationCostComputer):
    def __init__(self, recombination_rate, eff_pop_size):
        self._recombination_rate = recombination_rate
        self._eff_pop_size = eff_pop_size

    @staticmethod
    def uniform_recombination_map(recombrate, eff_pop_size, positions):

        # For a list of positions and a constant recombination rate (in cM/Mb),
        # return a list "results" of the same length as "positions" such that
        # results[i] is the phred-scaled recombination probability between
        # positions[i-1] and positions[i].
        
        return [(positions[i] - positions[i - 1])*recombrate*eff_pop_size*(4/(pow(10,6))) for i in range(1, len(positions))]

    def compute(self, positions):
        return self.uniform_recombination_map(self._recombination_rate, self._eff_pop_size, positions)
