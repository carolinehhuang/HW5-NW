# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    test = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    score, seq1_align, seq2_align = test.align(seq1, seq2)

    check_alignment = np.array([[0., -11., -12., -13.],
                                [-11., 5., -6., -7.],
                                [-12., -6., 4., -7.],
                                [-13., -7., -1., 5.],
                                [-14., -8., -6., 4.]])

    check_bt = np.array([[-1., 2., 2., 2.],
                         [1., 0., 2., 2.],
                         [1., 1., 0., 2.],
                         [1., 1., 0., 0.],
                         [1., 1., 0., 0.]])

    assert seq1_align == "MYQR", 'incorrect alignment'
    assert seq2_align == "M-QR", 'incorrect alignment'
    assert np.array_equal(test._align_matrix, check_alignment), 'incorrect alignment matrix'
    assert np.array_equal(test._back, check_bt), 'incorrect backtrace matrix'

    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    test_bt = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    score, seq3_align, seq4_align = test_bt.align(seq3, seq4)

    check_seq3 = "MAVHQLIRRP"
    check_seq4 = "M---QLIRHP"

    #check that project's aligned seq matches actual aligned sequence
    assert seq3_align == check_seq3, "Aligned seq3 is incorrect."
    assert seq4_align == check_seq4, "Aligned seq4 is incorrect."

    assert score == 17, "Score is incorrect."




