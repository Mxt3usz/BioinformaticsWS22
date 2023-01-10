from typing import List, Tuple, Dict

########################################################
############## Programming tasks #######################
########################################################


def zero_init(seq1, seq2):
    """
    Exercise 4 a
    Implement the function zero_init() which takes two sequences S1 and S2 and
    creates the helper matrix and initiates all the matrix values
    with zeroes. You can then use this matrix for D, P and Q matrices.
    Hereby S1 should be represented by the rows and S2 by the columns.
    """
    return [[0 for i in range(len(seq2)+1)] for i in range(len(seq1)+1)]

import math

def construct_matrix(seq1, seq2,value,scoring = {}):
    i, j = 1 , 1
    matrix = zero_init(seq1,seq2)
    place_i = "-" if value == "p" else math.inf
    place_j = math.inf if value == "p" else "-"

    while j != len(matrix[0]) or i != len(matrix):
        if j < len(matrix[0]) :
            matrix[0][j] = scoring["gap_introduction"] + j * scoring["gap_extension"] if value == "d" else place_j
            j += 1
        if i < len(matrix) :
            matrix[i][0] = scoring["gap_introduction"] + i * scoring["gap_extension"] if value == "d" else place_i
            i += 1
    return matrix



def d_matrix_init(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 4 b
    Implement the function d_matrix_init which takes two sequences S1 and S2 and
    the scoring schema and initializes the D matrix of the Gotoh algorithm.
    Hereby S1 should be represented by the rows and S2 by the columns.
    """
    return construct_matrix(seq1,seq2,"d",scoring)



def p_matrix_init(seq1, seq2):
    """
    Exercise 4 c
    Implement the function p_matrix_init which takes two sequences S1 and S2 and
    initializes the P matrix of the Gotoh algorithm.
    Hereby S1 should be represented by the rows and S2 by the columns.
    Use inf for the infinity and '-' to indicate the empty values to complete this matrix.
    """
    return construct_matrix(seq1,seq2,"p")



def q_matrix_init(seq1, seq2):
    """
    Exercise 4 d
    Implement the function q_matrix_init which takes two sequences S1 and S2 and
    initializes the Q matrix of the Gotoh algorithm.
    Hereby S1 should be represented by the rows and S2 by the columns.
    Use inf for the infinity and '-' to indicate the empty values to complete this matrix.
    """
    return construct_matrix(seq1,seq2,"q")



def gotoh_init(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 4 e
    Implement the function gotoh_init() to complete Gotoh initialization step
    Return the D, P and Q matrices with the corresponding initial values.
    """
    d_matrix, p_matrix, q_matrix = d_matrix_init(seq1, seq2, scoring), p_matrix_init(seq1, seq2), q_matrix_init(seq1, seq2)
    return d_matrix, p_matrix, q_matrix


def gotoh_forward(seq1, seq2, scoring: Dict[str, int]):
    """
    Exercise 4 f
    Implement the function gotoh_forward() which takes two sequences S1 and S2 as
    well as the scoring function and fills in all the values in D, P and Q matrices
    """
    d_matrix, p_matrix, q_matrix = gotoh_init(seq1, seq2, scoring)
    i,j = 1, 1
    while i != len(d_matrix):
        # calculate matrix entries
        p_matrix[i][j] = min(d_matrix[i-1][j] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],p_matrix[i-1][j] + scoring["gap_extension"])
        q_matrix[i][j] = min(d_matrix[i][j-1] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],q_matrix[i][j-1] + scoring["gap_extension"])
        d_matrix[i][j] = min(d_matrix[i-1][j-1] + (scoring["match"] if seq1[i-1] == seq2[j-1] else scoring["mismatch"]),p_matrix[i][j],q_matrix[i][j])
        j += 1
        if j == len(d_matrix[0]):
            j = 1
            i += 1

    return d_matrix, p_matrix, q_matrix


def previous_cells(
    seq1, seq2, scoring, d_matrix, p_matrix, q_matrix, cell: Tuple[str, Tuple[int, int]]
) -> List[Tuple[str, Tuple[int, int]]]:

    """
    Exercise 4 g
    Implement the function previous_cells() which takes two sequences S1 and
    S2, scoring function, the filled in recursion matrices from the step f) and
    the cell coordinates (matrix, (row, column)) i.e. ("D", (1, 3)). The function should output the list
    of all possible previous cells.
    """
    trace = []
    i = cell[1][0]
    j = cell[1][1]
    m = cell[0]
    if m == "D":
        if (i == 0 and j != 0):
                trace += [("D",(i,j-1))]
                return trace
        if (j == 0 and i != 0):
                trace += [("D",(i-1,j))]
                return trace
        if (i == 0 and j == 0) or (i == 1 and j == 0) or (i == 0 and j == 1):
                return trace
        d,p,q = d_matrix[i-1][j-1] + (scoring["match"] if seq1[i-1] == seq2[j-1] else scoring["mismatch"]),p_matrix[i][j],q_matrix[i][j] 
        if d == d_matrix[i][j]:
            trace += [("D",(i-1,j-1))]
        if p == d_matrix[i][j]:
            trace += [("P",(i,j))]
        if q == d_matrix[i][j]:
            trace += [("Q",(i,j))]
    if m == "P":
        d,p = d_matrix[i-1][j] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],p_matrix[i-1][j] + scoring["gap_extension"]
        if d == p_matrix[i][j]:
            trace += [("D",(i-1,j))]
        if p == p_matrix[i][j]:
            trace += [("P",(i-1,j))]
    if m == "Q":
        d,q = d_matrix[i][j-1] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],q_matrix[i][j-1] + scoring["gap_extension"]
        if d == q_matrix[i][j]:
            trace += [("D",(i,j-1))]
        if q == q_matrix[i][j]:
            trace += [("Q",(i,j-1))]
    return trace

        
    
def build_all_traceback_paths(seq1, seq2, scoring, d_matrix, p_matrix, q_matrix) -> \
        List[List[Tuple[str, Tuple[int, int]]]]:
    """
    Exercise 4 h
    Implement the function which builds all possible traceback paths.
    """
    global curr_trace
    curr_trace = [("D",(len(d_matrix)-1,len(d_matrix[0])-1))]
    def rec_traceback(i,j,m,tracebacks = []):
        global curr_trace
        # with d matrix there are 3 cases the cell (int,int) could orginate from
        if m == "D":
            # Handle i gaps
            if (i == 0 and j != 0):
                curr_trace += [("D",(i,j-1))]
                rec_traceback(i,j-1,"D",tracebacks)
                curr_trace = curr_trace[:-1]
                return
            # Handle j gaps
            if (j == 0 and i != 0):
                curr_trace += [("D",(i-1,j))]
                rec_traceback(i-1,j,"D",tracebacks)
                curr_trace = curr_trace[:-1]
                return
            # We reached one of the starting entries (terminate traceback)
            if (i == 0 and j == 0) or (i == 1 and j == 0) or (i == 0 and j == 1):
                    tracebacks += [curr_trace]
                    return
            # to determine where it came from we can look at the min values of each matrix, it should equal the score of the current cell
            d,p,q = d_matrix[i-1][j-1] + (scoring["match"] if seq1[i-1] == seq2[j-1] else scoring["mismatch"]),p_matrix[i][j],q_matrix[i][j] 
            # depending on where the current score came from we trace back to there
            if d == d_matrix[i][j]:
                curr_trace += [("D",(i-1,j-1))]
                rec_traceback(i-1,j-1,"D",tracebacks)
                curr_trace = curr_trace[:-1]
            if p == d_matrix[i][j]:
                curr_trace += [("P",(i,j))]
                rec_traceback(i,j,"P",tracebacks)
                curr_trace = curr_trace[:-1]
            if q == d_matrix[i][j]:
                curr_trace += [("Q",(i,j))]
                rec_traceback(i,j,"Q",tracebacks)
                curr_trace = curr_trace[:-1]
        if m == "P":
            d,p = d_matrix[i-1][j] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],p_matrix[i-1][j] + scoring["gap_extension"]
            if d == p_matrix[i][j]:
                curr_trace += [("D",(i-1,j))]
                rec_traceback(i-1,j,"D",tracebacks)
                curr_trace = curr_trace[:-1]
            if p == p_matrix[i][j]:
                curr_trace += [("P",(i-1,j))]
                rec_traceback(i-1,j,"P",tracebacks)
                curr_trace = curr_trace[:-1]
        if m == "Q":
            d,q = d_matrix[i][j-1] + scoring["gap_introduction"] + 1 * scoring["gap_extension"],q_matrix[i][j-1] + scoring["gap_extension"]
            if d == q_matrix[i][j]:
                curr_trace += [("D",(i,j-1))]
                rec_traceback(i,j-1,"D",tracebacks)
                curr_trace = curr_trace[:-1]
            if q == q_matrix[i][j]:
                curr_trace += [("Q",(i,j-1))]
                rec_traceback(i,j-1,"Q",tracebacks)
                curr_trace = curr_trace[:-1]
        return tracebacks
    return rec_traceback(len(d_matrix)-1,len(d_matrix[0])-1,"D")


def build_alignment(seq1, seq2, traceback_path) -> Tuple[str, str]:
    """
    Exercise 4 i
    Implement the function build_alignment() which takes two sequences and
    outputs the alignment.
    """
    seq1_al = ""
    seq2_al = ""
    c = 0
    while c != len(traceback_path)-1:
        # diagonal
        if traceback_path[c][1][0] > traceback_path[c+1][1][0] and traceback_path[c][1][1] > traceback_path[c+1][1][1]:
            seq1_al += seq1[traceback_path[c][1][0]-1]
            seq2_al += seq2[traceback_path[c][1][1]-1]
        # left
        elif traceback_path[c][1][1] > traceback_path[c+1][1][1]:
            seq1_al += "-"
            seq2_al += seq2[traceback_path[c][1][1]-1]
        # up
        elif traceback_path[c][1][0] > traceback_path[c+1][1][0]:
            seq1_al += seq1[traceback_path[c][1][0]-1]
            seq2_al += "-"
        c += 1
    return seq1_al[::-1],seq2_al[::-1]


def display_matrices(seq1,seq2,matrices):
    string = ""
    string += "     "
    for base in seq2:
        string += base + "  " 
    seq2 = string
    string = ""
    c = 0
    for matrix in matrices:
        s = -1
        if c == 0:
            string += "(Di,j) :" + "\n" + seq2 + "\n"
        if c == 1:
            string += "(Pi,j) :" + "\n" + seq2 + "\n"
        if c == 2:
            string += "(Qi,j) :" + "\n" + seq2 + "\n"
        c += 1
        for row in matrix:
            if s != -1:
                string += seq1[s]
            else:
                string += " "
            string += str(row) + "\n"
            s += 1
    return string




if __name__ == "__main__":

    seq1 = "CC"
    seq2 = "ACCT"
    scoring = {"match": 0, "mismatch": 1, "gap_introduction": 4, "gap_extension": 1}
    print(display_matrices(seq1,seq2,gotoh_forward(seq1,seq2,scoring)))
    b = build_all_traceback_paths(seq1,seq2,scoring,gotoh_forward(seq1,seq2,scoring)[0],gotoh_forward(seq1,seq2,scoring)[1],gotoh_forward(seq1,seq2,scoring)[2])
    print(b)
    print(build_alignment(seq1,seq2,b[0]))
