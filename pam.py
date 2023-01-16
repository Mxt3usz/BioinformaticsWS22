from typing import Tuple, List
import numpy as np
import math

def nucleotide_freq_calculation(sym_align: Tuple[str, str]) -> List[float]:
    """
    Implement the calculation of the nucleotide frequencies for a symmetric
    alignment. Make sure that the indexing works the followin way:
    [0] = A
    [1] = C
    [2] = T
    [3] = G
    sym align: Tuple of strings representing the symmetric alignment
    """
    freq_list = [0,0,0,0]  # Frequencies should be in order A,C,T,G
    base_to_num = {"A" : 0, "C" : 1, "T" : 2, "G" : 3}
    for base in sym_align[0]:
        freq_list[base_to_num[base]] += 1/len(sym_align[0])
    return freq_list


def mutation_calculation(sym_align: Tuple[str, str]) -> List[List[float]]:
    """
    Implement the calculation of the mutations remember indexing in any
    dimension should be done in order:
    [0] = A
    [1] = C
    [2] = T
    [3] = G
    Hint: You can use your already implemented functions here
    """
    base_to_num = {"A" : 0, "C" : 1, "T" : 2, "G" : 3}
    first_seq, second_seq = sym_align
    freq_matrix = [[0] * 4 for _ in range(4)]
    for i in range(len(first_seq)):
        freq_matrix[base_to_num[first_seq[i]]][base_to_num[second_seq[i]]] += 1/len(first_seq)
    return freq_matrix
 

def scores_calculation(sym_align: Tuple[str, str]) -> List[List[int]]:
    """
    Implement the calculation of the scores remember indexing in any
    dimension should be done in order:
    [0] = A
    [1] = C
    [2] = T
    [3] = G
    Hint: You can use your already implemented functions here
    """
    mutation_matrix = mutation_calculation(sym_align)
    freq = nucleotide_freq_calculation(sym_align)
    i,j= 0,0
    while i != 4:
        mutation_matrix[i][j] = round(10 * math.log10(mutation_matrix[i][j]/(freq[i] * freq[j])))
        j += 1
        if j == 4:
            i += 1
            j = 0
    return mutation_matrix


def gamma_calculation(sym_align: Tuple[str, str]) -> float:
    """
    Implement the calculation of gamma.
    Hint: You can use your already implemented functions here
    """
    mutation = mutation_calculation(sym_align)
    return 0.01/sum([sum(mutation[row]) - mutation[row][row] for row in range(len(mutation))]) # subtract exx to only get exy


def probabilities_calculation(sym_align: Tuple[str, str]) -> List[List[float]]:
    """
    Implement the calculation of probabilities matrix.
    Hint: You can use your already implemented functions here
    """
    mutation_matrix = mutation_calculation(sym_align)
    freq = nucleotide_freq_calculation(sym_align)
    i,j= 0,0
    while i != 4:
        mutation_matrix[i][j] = mutation_matrix[i][j]/freq[i]
        j += 1
        if j == 4:
            i += 1
            j = 0
    return mutation_matrix


def norm_probabilities_calculation(sym_align: Tuple[str, str]) -> List[List[float]]:
    """
    Implement the calculation of normalized probabilities matrix.
    Hint: You can use your already implemented functions here
    """
    prob_matrix = probabilities_calculation(sym_align)
    gamma = gamma_calculation(sym_align)
    i,j= 0,0
    while i != 4:
        if i != j:
            prob_matrix[i][j] = gamma * prob_matrix[i][j] # exy
        j += 1
        if j == 4:
            prob_matrix[i][i] = 1 - (sum(prob_matrix[i]) - prob_matrix[i][i]) # exx
            i += 1 
            j = 0
    return prob_matrix



def pam_calculation(sym_align: Tuple[str, str], power: int) -> List[List[int]]:
    """
    Implement the calculation of pam matrix Make sure to round values to
    integers. ( Notice that casting to int does not do the correct
    rounding)
    power: the power of the PAM matrix e.g. PAM1 or PAM250
    Hint: You can use your already implemented functions here
    """
    norm_prob = np.linalg.matrix_power(norm_probabilities_calculation(sym_align), power)
    norm_prob = [[col for col in row] for row in norm_prob] # convert from numpy
    freq = nucleotide_freq_calculation(sym_align)
    i,j= 0,0
    while i != 4:
        norm_prob[i][j] = round(10 * math.log10(norm_prob[i][j]/freq[j]))
        j += 1
        if j == 4:
            i += 1
            j = 0
    return norm_prob


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

def main():
    seq1 = "AAGTACTTTAGGTAACACGTTTAGTCAAAATTCCTAAGTTTACCGGGTTAATCA"
    seq2 = "AAATTCCTAAGTTTACCGGGTTAATCAAAGTACTTTAGGTAACACGTTTAGTCA"

    sim_align = seq1, seq2

    print(string_matrix(mutation_calculation(sim_align),"ACTG","ACTG"))
    print(string_matrix(scores_calculation(sim_align),"ACTG","ACTG"))
    print(gamma_calculation(sim_align))
    print(string_matrix(probabilities_calculation(sim_align),"ACTG","ACTG"))
    print(string_matrix(norm_probabilities_calculation(sim_align),"ACTG","ACTG"))
    print(string_matrix(pam_calculation(sim_align,1),"ACTG","ACTG"))
 

if __name__ == "__main__":
    main()