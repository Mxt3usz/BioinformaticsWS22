def zero_init_correct(seq1, seq2):
    return [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]


def nw_init_correct(seq1, seq2, scoring):
    match, mismatch, gap = (
        scoring["match"],
        scoring["mismatch"],
        scoring["gap_introduction"],
    )
    matrix = zero_init_correct(seq1, seq2)
    first_row = [i * gap for i in range(len(seq2) + 1)]
    matrix[0] = first_row
    for index, column in enumerate(matrix):
        column[0] = index * gap
    return matrix


def nw_init_csv_maker(seq1, seq2, scoring, csv_file_name="matrix.csv"):
    matrix = nw_init_correct(seq1, seq2, scoring)
    with open(csv_file_name, "w") as f:
        f.write(",".join("  " + seq2) + "\n")

        left_column = " " + seq1
        for index, char in enumerate(left_column):
            f.write(
                char + "," + ",".join([str(x) for x in matrix[index]]) + "\n"
            )


def given_matrix_csv_maker(seq1, seq2, matrix, csv_file_name):
    with open(csv_file_name, "w") as f:
        f.write(",".join("  " + seq2) + "\n")

        left_column = " " + seq1
        for index, char in enumerate(left_column):
            f.write(
                char + "," + ",".join([str(x) for x in matrix[index]]) + "\n"
            )


def nw_forward_correct(seq1, seq2, scoring):
    matrix = nw_init_correct(seq1, seq2, scoring)
    match_score, mismatch_score, gap_score = (
        scoring["match"],
        scoring["mismatch"],
        scoring["gap_introduction"],
    )

    for row_index, row in enumerate(matrix[1:], 1):
        for column_index, column in enumerate(row[1:], 1):
            char_seq_1, char_seq_2 = (
                seq1[row_index - 1],
                seq2[column_index - 1],
            )
            no_gap_score = (
                match_score if char_seq_1 == char_seq_2 else mismatch_score
            )

            diagonal = matrix[row_index - 1][column_index - 1] + no_gap_score
            left = matrix[row_index][column_index - 1] + gap_score
            top = matrix[row_index - 1][column_index] + gap_score
            min_val = min(top, left, diagonal)
            matrix[row_index][column_index] = min_val
    return matrix


def previous_cells_correct(seq1, seq2, scoring, nw_matrix, cell):
    prev_cells = []
    row, column = cell

    top = (row - 1, column) if row > 0 else None
    left = (row, column - 1) if column > 0 else None
    diagonal = (row - 1, column - 1) if (row > 0 and column > 0) else None

    cur_val = nw_matrix[row][column]
    char_first, char_second = seq1[row - 1], seq2[column - 1]
    match_score = (
        scoring["match"] if char_first == char_second else scoring["mismatch"]
    )
    gap_score = scoring["gap_introduction"]

    if diagonal:
        diagonal_val = nw_matrix[diagonal[0]][diagonal[1]]
        if (diagonal_val + match_score) == cur_val:
            prev_cells.append(diagonal)

    if top:
        top_val = nw_matrix[top[0]][top[1]]
        if (top_val + gap_score) == cur_val:
            prev_cells.append(top)

    if left:
        left_val = nw_matrix[left[0]][left[1]]
        if (left_val + gap_score) == cur_val:
            prev_cells.append(left)

    return prev_cells


def build_all_traceback_paths_correct(seq1, seq2, scoring, nw_matrix):
    list_traceback_paths = []

    cell = len(nw_matrix) - 1, len(nw_matrix[0]) - 1
    frontier = [[cell]]
    while frontier:
        partial_path = frontier.pop()
        last_cell_partial = partial_path[-1]
        next_steps = previous_cells_correct(
            seq1, seq2, scoring, nw_matrix, last_cell_partial
        )
        for next_step in next_steps:
            new_traceback_path = partial_path + [next_step]
            if next_step == (0, 0):
                list_traceback_paths.append(new_traceback_path)
            else:
                frontier.append(new_traceback_path)

    return list_traceback_paths


def build_alignment_correct(seq1, seq2, alignment_path):
    align_seq1 = ""
    align_seq2 = ""

    alignment_path = alignment_path[::-1]

    prev_cell = 0, 0
    for cell in alignment_path[1:]:
        prev_row, prev_column = prev_cell
        row, column = cell

        if (row > prev_row) and (column > prev_column):
            align_seq1 += seq1[row - 1]
            align_seq2 += seq2[column - 1]

        elif row > prev_row:
            align_seq1 += seq1[row - 1]
            align_seq2 += "-"

        else:
            align_seq1 += "-"
            align_seq2 += seq2[column - 1]

        prev_cell = cell

    return align_seq1, align_seq2