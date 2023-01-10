from typing import List, Tuple, Dict


########################################################
############## Programming tasks #######################
########################################################


def sw_init(seq1, seq2):
    """
    Exercise 3 a
    Implement the function sw_init() which takes two sequences S1 and S2 and
    creates the Smith-Waterman matrix and initiates all the matrix values
    with zeroes. Hereby S1 should be represented by the rows and S2 by
    the columns.
    """
    return [[0 for i in range(len(seq2)+1)] for i in range(len(seq1)+1)]


def sw_forward(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 3 b
    Implement the function sw_forward() which takes the two sequences S1 and
    S2 and the scoring function and output the complete matrix filled with
    the Smith-Waterman approach.
    """
    matrix = sw_init(seq1,seq2)
    i,j = 1,1
    while i != len(matrix):
        matrix[i][j] = max(matrix[i-1][j-1] + (scoring["match"] if seq1[i-1] == seq2[j-1] else scoring["mismatch"]), matrix[i-1][j] + scoring["gap_introduction"], matrix[i][j-1] + scoring["gap_introduction"],0)
        j += 1
        if j == len(matrix[0]):
            i += 1
            j = 1
    return matrix


def previous_cells(
    seq1, seq2, scoring, sw_matrix, cell: Tuple[int, int]
) -> List[Tuple[int, int]]:
    """
    Exercise 3 c
    Implement the function previous_cells() which takes two sequences S1 and
    S2, scoring function, the filled in recursion matrix from the step c) and
    the cell coordinates (row, column). The function should output the list
    of all possible previous cells.
    """
    cells = []
    i,j = cell[0],cell[1]
    if j > 0:
        left = sw_matrix[i][j-1] + scoring["gap_introduction"]
        if left == sw_matrix[i][j]:
            cells += [(i,j-1)]
    if i > 0:
        top = sw_matrix[i-1][j] + scoring["gap_introduction"]
        if top == sw_matrix[i][j]:
            cells += [(i-1,j)]
    if i > 0 and j > 0:
        diagonal = sw_matrix[i-1][j-1] + (scoring["match"] if seq1[i-1] == seq2[j-1] else scoring["mismatch"])
        if diagonal == sw_matrix[i][j]:
            cells += [(i-1,j-1)]
    return cells

def find_all_max(matrix):
    maxi = 0
    for row in matrix:
        for col in row:
            if col > maxi:
                maxi = col
    if maxi == 0:
        return []
    max_lst = []
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            if matrix[row][col] == maxi:
                max_lst += [(row,col,0)]
    return max_lst

        
def build_all_traceback_paths(
    seq1, seq2, scoring, sw_matrix
) -> List[List[Tuple[int, int]]]:
    """
    Exercise 3 d
    Implement the function which builds all possible traceback paths.
    """
    queue =  find_all_max(sw_matrix)
    curr_trace = []
    trace = []
    c = 0
    while queue != []:
        curr = queue[-1]
        curr_trace += [(curr[0],curr[1])]
        queue = queue[:-1]
        c += 1
        if sw_matrix[curr[0]][curr[1]] == 0:
            trace += [curr_trace]
            if queue != [] and queue[-1][2] == 0:
                curr_trace = []
                c = 0
            else:
                if queue != []:
                    c = queue[-1][2]
                    curr_trace = curr_trace[:c]
            continue
        for prev in previous_cells(seq1,seq2,scoring,sw_matrix,curr):
            queue += [(prev[0],prev[1],c)]

    return trace

def build_alignment(seq1, seq2, traceback_path) -> Tuple[str, str]:
    """
    Exercise 3 e
    Implement the function build_alignment() which takes two sequences and
    outputs the alignment.
    """
    seq1_aligned = ""
    seq2_aligned = ""
    c = 0
    if traceback_path == []:
        return ("","")
    while c != len(traceback_path)-1:
        if traceback_path[c][0] > traceback_path[c+1][0] and traceback_path[c][1] > traceback_path[c+1][1]:
            seq1_aligned += seq1[traceback_path[c][0]-1]
            seq2_aligned += seq2[traceback_path[c][1]-1]
        elif traceback_path[c][0] > traceback_path[c+1][0]:
            seq1_aligned += seq1[traceback_path[c][0]-1]
            seq2_aligned += "-"
        else:
            seq1_aligned += "-"
            seq2_aligned += seq2[traceback_path[c][1]-1]
        c += 1
    return seq1_aligned[::-1],seq2_aligned[::-1]

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



if __name__ == "__main__":

    scoring = {"match": 2, "mismatch": -1, "gap_introduction": -2}
    seq1 = 'AAA'
    seq2 = 'TTT'
    print(string_matrix(sw_forward(seq1,seq2,scoring),seq1,seq2))
    print(build_all_traceback_paths(seq1,seq2,scoring,sw_forward(seq1,seq2,scoring)))
    print(build_alignment(seq1,seq2,build_all_traceback_paths(seq1,seq2,scoring,sw_forward(seq1,seq2,scoring))))
