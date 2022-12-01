from matrix_helpers import nw_init_csv_maker
from typing import Dict, List, Tuple

########################################################
############## Programming tasks #######################
########################################################


def zero_init(seq1, seq2):
    """
    Exercise 4 a
    Implement the function zero_init() which takes two sequences S1 and S2 and
    creates the Needleman-Wunsch matrix and initiates all the matrix values
    with zeroes. Hereby S1 should be represented by the rows and S2 by
    the columns.
    """
    return [[0 for i in range(len(seq2)+1)] for i in range(len(seq1)+1)]

def string_matrix(matrix,seq1,seq2):
    """
    Represent the Needleman-Wunsch matrix in a more readable format
    """
    nw_matrix = "   -  "
    for char in seq2:
        nw_matrix += char + "  "
    nw_matrix += "\n"
    c = 0
    while c != len(matrix):
        nw_matrix += seq1[c-1] + " " if c != 0 else "-" + " "
        nw_matrix += str(matrix[c]) + "\n" if c != len(matrix) - 1 else str(matrix[c])
        c += 1
    return nw_matrix

def nw_init(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 4 b
    Implement the function nw_init() which takes two sequences S1 and S2 as
    well as the scoring function and fills in the values for the first row and
    first column of the matrix with the correct values. Utilize a) in your
    implementation.
    """
    matrix = zero_init(seq1, seq2)
    c = 0
    while c != len(matrix):
        if c == 0 :
            for i in range(1,len(matrix[c])):
                matrix[c][i] =  scoring["gap_introduction"] + matrix[c][i-1]
        else:
            matrix[c][0] =  scoring["gap_introduction"] + matrix[c-1][0]
        c += 1
    return matrix

def nw_forward(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 4 c
    Implement the function nw_forward() which takes the two sequences S1 and
    S2 and the scoring function and output the complete matrix filled with
    the Needleman-Wunsch approach.
    """
    nw_matrix = nw_init(seq1,seq2,scoring)
    i = 1
    j = 1
    while i != len(nw_matrix):
        min_gap =  nw_matrix[i-1][j] + scoring["gap_introduction"] if nw_matrix[i-1][j] < nw_matrix[i][j-1] else nw_matrix[i][j-1] + scoring["gap_introduction"]
        mis_or_match = nw_matrix[i-1][j-1] + scoring["match"] if seq1[i-1] == seq2[j-1] else nw_matrix[i-1][j-1] + scoring["mismatch"]
        nw_matrix[i][j] = mis_or_match if mis_or_match < min_gap else min_gap
        j += 1 if j != len(nw_matrix[i])-1 else -j+1
        i += 1 if j == 1 else 0
    return nw_matrix

def previous_cells(
    seq1, seq2, scoring, nw_matrix, cell: Tuple[int, int]
) -> List[Tuple[int, int]]:
    """
    Exercise 4 d
    Implement the function previous_cells() which takes two sequences S1 and
    S2, scoring function, the filled in recursion matrix from the step c) and
    the cell coordinates (row, column). The function should output the list
    of all possible previous cells.
    """
    if cell == (0,0):
        return []
    if cell == ((0,1) or (1,0)):
        return [(0,0)]
    cells = []
    i = cell[0]
    j = cell[1]
    left_gap =  nw_matrix[i-1][j] + scoring["gap_introduction"]
    right_gap = nw_matrix[i][j-1] + scoring["gap_introduction"]
    mis_or_match = nw_matrix[i-1][j-1] + scoring["match"] if seq1[i-1] == seq2[j-1] else nw_matrix[i-1][j-1] + scoring["mismatch"]
    if left_gap == nw_matrix[i][j]:
        cells += [(i-1,j)]
    if right_gap == nw_matrix[i][j]:
        cells += [(i,j-1)]
    if mis_or_match == nw_matrix[i][j]:
        cells += [(i-1,j-1)]
    return cells
    
def build_all_traceback_paths(
    seq1, seq2, scoring, nw_matrix
) -> List[List[Tuple[int, int]]]:
    """
    Exercise 4 e
    Implement the function which builds all possible traceback paths.
    """
    
    global pa
    pa = [(len(seq1),len(seq2))]
    def recursive_traceback(cell,trace_back_paths = []):
        global pa
        previous_cell = previous_cells(seq1,seq2,scoring,nw_matrix,cell)
        for cells in previous_cell:
            pa += [cells]
            recursive_traceback(cells)
        trace_back_paths += [pa] if pa != [] and pa[-1] == (0,0) else []
        pa = pa[:-1]
        return trace_back_paths
    return recursive_traceback((len(seq1),len(seq2)))

    
def build_alignment(seq1, seq2, traceback_path) -> Tuple[str, str]:
    """
    Exercise 4 f
    Implement the function build_alignment() which takes two sequences and
    outputs the alignment.
    """
    allignement1 = ""
    allignement2 = ""
    c = 0
    letters_left1 = len(seq1)-1
    letters_left2 = len(seq2)-1
    while c != len(traceback_path)-1:
        if traceback_path[c][0] - traceback_path[c+1][0] == 1 and traceback_path[c][1] - traceback_path[c+1][1] == 0:
             allignement1 += seq1[letters_left1] if letters_left1 >= 0 else ""
             allignement2 += "-"
             letters_left1 -= 1
        elif traceback_path[c][0] - traceback_path[c+1][0] == 0 and traceback_path[c][1] - traceback_path[c+1][1] == 1:
            allignement1 += "-"
            allignement2 += seq2[letters_left2] if letters_left2 >= 0 else ""
            letters_left2 -= 1
        else:
            allignement1 += seq1[letters_left1] if letters_left1 >= 0 else ""
            allignement2 += seq2[letters_left2] if letters_left2 >= 0 else ""
            letters_left1 -= 1
            letters_left2 -= 1
        c += 1
    return (allignement1[::-1],allignement2[::-1])


if __name__ == "__main__":
    """
    You can run this to create csv files from two sequences. Further you can 
    import this file in excel or some similar program, where you can fill in
    the forward values yourself.    
    """
    seq1 = "TCCGA"
    seq2 = "TACGCGC"
    scoring = {"match": -1, "mismatch": 0, "gap_introduction": 1}

    print(string_matrix(nw_forward(seq1,seq2,scoring),seq1,seq2))
    print(build_all_traceback_paths(seq1,seq2,scoring,nw_forward(seq1,seq2,scoring)))
    print(build_alignment(seq1,seq2,build_all_traceback_paths(seq1,seq2,scoring,nw_forward(seq1,seq2,scoring))[0]))