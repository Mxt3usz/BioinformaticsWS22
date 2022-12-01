from exercise_sheet3 import *
import pytest
import pprint
from matrix_helpers import (
    nw_init_correct,
    given_matrix_csv_maker,
    nw_forward_correct,
    zero_init_correct,
    build_all_traceback_paths_correct,
    previous_cells_correct,
    build_alignment_correct,
)

SCORING = [
    {"match": -1, "mismatch": 0, "gap_introduction": 1},
    {"match": -1, "mismatch": 0, "gap_introduction": 2},
    {"match": -1, "mismatch": 1, "gap_introduction": 2},
    {"match": -1, "mismatch": 2, "gap_introduction": 1},
]


@pytest.mark.parametrize(
    "seq1,seq2",
    [
        ("TCCGA", "TACGCGC")
    ]
)
def test_exercise_4a(seq1, seq2):
    expected_matrix = zero_init_correct(seq1, seq2)
    actual_matrix = zero_init(seq1, seq2)
    if actual_matrix != expected_matrix:
        print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}")
        print(f"Your matrix is:")
        pprint.pprint(actual_matrix)
        print("It is supposed to look like:")
        pprint.pprint(expected_matrix)
    assert expected_matrix == actual_matrix


@pytest.mark.parametrize(
    "seq1,seq2,scoring",
    [
        ("TCCGA", "TACGCGC", SCORING[0]),
        ("TCCGA", "TACGCGC", SCORING[1]),
        ("TCCGA", "TACGCGC", SCORING[2]),
    ]
)
def test_exercise_4b(seq1, seq2, scoring):
    expected_matrix = nw_init_correct(seq1, seq2, scoring)
    actual_matrix = nw_init(seq1, seq2, scoring)
    if actual_matrix != expected_matrix:
        print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}\nscoring: {scoring}")
        print(f"Your matrix is:")
        pprint.pprint(actual_matrix)
        print("It is supposed to look like:")
        pprint.pprint(expected_matrix)
    assert expected_matrix == actual_matrix


@pytest.mark.parametrize(
    "seq1,seq2,scoring",
    [
        ("TCCGA", "TACGCGC", SCORING[0]),
        ("TCCGA", "TACGCGC", SCORING[1]),
        ("TCCGA", "TACGCGC", SCORING[2]),
    ]
)
def test_exercise_4c(seq1, seq2, scoring):
    expected_matrix = nw_forward_correct(seq1, seq2, scoring)
    actual_matrix = nw_forward(seq1, seq2, scoring)
    if actual_matrix != expected_matrix:
        print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}\nscoring: {scoring}")
        print(f"Your matrix is:")
        pprint.pprint(actual_matrix)
        print("It is supposed to look like:")
        pprint.pprint(expected_matrix)
    assert expected_matrix == actual_matrix

@pytest.mark.parametrize(
    "seq1,seq2,scoring,cell",
    [
        ("TCCGA", "TACGCGC", SCORING[0], (1, 1)),
        ("TCCGA", "TACGCGC", SCORING[1], (5, 7)),
        ("TCCGA", "TACGCGC", SCORING[2], (3, 4)),
    ]
)
def test_exercise_4d(seq1, seq2, scoring, cell):
    nw_matrix = nw_forward_correct(seq1, seq2, scoring)
    expected_cells = previous_cells_correct(
        seq1, seq2, scoring, nw_matrix, cell
    )
    actual_cells = previous_cells(seq1, seq2, scoring, nw_matrix, cell)
    expected_cells = set(expected_cells)
    actual_cells = set(actual_cells)
    if actual_cells != expected_cells:
        print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}\ncell: {cell}\n"
              f"scoring: {scoring}")
        print(f"Your previous cells are:\n{actual_cells}")
        print(f"It should be:\n{expected_cells}")
    assert expected_cells == actual_cells


@pytest.mark.parametrize(
    "seq1,seq2,scoring",
    [
        ("TCCGA", "TACGCGC", SCORING[0]),
        ("TCCCGG", "TCAAA", SCORING[0]),
        ("TCCGA", "TACGCGC", SCORING[1]),
        ("TCCGA", "TACGCGC", SCORING[2]),
        ("AAA", "TTT", SCORING[3]),
    ]
)


def test_exercise_4e(seq1, seq2, scoring):
    nw_matrix = nw_forward_correct(seq1, seq2, scoring)
    expected_paths = build_all_traceback_paths_correct(
        seq1, seq2, scoring, nw_matrix
    )
    actual_paths = build_all_traceback_paths(seq1, seq2, scoring, nw_matrix)
    expected_paths = [tuple(x) for x in expected_paths]
    expected_paths = set(expected_paths)
    actual_paths = [tuple(x) for x in actual_paths]
    actual_paths = set(actual_paths)
    if actual_paths != expected_paths:
        print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}\nscoring: {scoring}")
        print(f"Your paths are:")
        pprint.pprint(actual_paths)
        print(f"It should be:")
        pprint.pprint(expected_paths)
    assert expected_paths == actual_paths


@pytest.mark.parametrize(
    "seq1,seq2,scoring",
    [
        ("TCCGA", "TACGCGC", SCORING[0]),
        ("TCCCGG", "TCAAA", SCORING[0]),
        ("TCCGA", "TACGCGC", SCORING[1]),
        ("TCCGA", "TACGCGC", SCORING[2]),
        ("AAA", "TTT", SCORING[3]),
    ]
)
def test_exercise_4f(seq1, seq2, scoring):
    nw_matrix = nw_forward_correct(seq1, seq2, scoring)
    expected_paths = build_all_traceback_paths_correct(
        seq1, seq2, scoring, nw_matrix
    )
    for traceback_path in expected_paths:
        expected_alignment = build_alignment_correct(seq1, seq2, traceback_path)
        actual_alignment = build_alignment(seq1, seq2, traceback_path)
        if actual_alignment != expected_alignment:
            print(f"\nFor the Test case:\nS1: {seq1}\nS2: {seq2}\nscoring: {scoring}\npath: {traceback_path}")
            print("Your alignment is:")
            print(f"S1: {actual_alignment[0]}\nS2: {actual_alignment[1]}")
            print("It should be:")
            print(f"S1: {expected_alignment[0]}\nS2: {expected_alignment[1]}")
        assert expected_alignment == actual_alignment


