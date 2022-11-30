


def levenshtein_substitution(sequence1, sequence2):
    """
    Implement the function levenshtein_substitution() which takes two sequences
    of the same length and computes the minimum number of substitutions to
    transform one into another.
    """
    number_substitutions = 0
    i = 0
    while i != len(sequence1):
        if sequence1[i] != sequence2[i]:
            number_substitutions += 1
        i += 1
    return number_substitutions


def levenshtein_deletions(sequence1, sequence2):
    """
    Implement the function levenshtein_deletion() which takes two sequences
    and returns the positions of characters from the longest sequences which
    should be deleted to transform the sequence into the other one.
    This should be returned as a list of indices (int). 
    If such deletion can not be done the function should return None. 
    Also, if there are no editing operations needed the function should return an empty list.
    """
    deletions_indexes = []
    i = 0
    i_static = 0
    if len(sequence1) > len(sequence2):
        shorter = sequence2
        longer = sequence1
    else:
        shorter = sequence1
        longer = sequence2
    
    while i != len(longer):
        if i >= len(longer):
            return None
        if i >= len(shorter):
            longer = longer[:i] + longer[i+1:]
            deletions_indexes += [i_static]
            i -= 1
        else:
            if longer[i] != shorter[i]:
                longer = longer[:i] + longer[i+1:]
                deletions_indexes += [i_static]
                i -= 1
        i += 1
        i_static += 1
    return deletions_indexes if longer == shorter else None

print(levenshtein_deletions("CGCGCTGCT","CGCTGT"))