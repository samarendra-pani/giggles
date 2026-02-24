"""
Run C++ Tests

Runs the tests for the core classes written in C++.
"""

from giggles.core import run_cpp_tests

def run_tests():
    run_cpp_tests()

def add_arguments(parser):
    pass

def validate(args, parser):
    pass
    
def main(args):
    run_tests(**vars(args))