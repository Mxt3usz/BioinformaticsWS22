import argparse
from matrix_helpers import nw_init_csv_maker, given_matrix_csv_maker, nw_forward_correct


parser = argparse.ArgumentParser(description='Parameters for the NW')
parser.add_argument('--first', type=str,
                    help='first sequence')
parser.add_argument('--second', type=str,
                    help='second sequence')

parser.add_argument('--match', type=int,
                    help='match score')

parser.add_argument('--mismatch', type=int,
                    help='mismatch score')

parser.add_argument('--gap_introduction', type=int,
                    help='gap introduction score')


parser.add_argument('--file_name', type=str, default="matrix.csv",
                    help='file name')

parser.add_argument('--check', type=bool, default=False,
                    help='if set to True, the output will be the '
                         'correctly filled matrix (default: False)')

args = parser.parse_args()

seq_first = args.first
seq_second = args.second

match = args.match
mismatch = args.mismatch
gap = args.gap_introduction

file_name = args.file_name

scoring = {"match": match, "mismatch": mismatch, "gap_introduction": gap}

if not args.check:
    nw_init_csv_maker(seq_first, seq_second, scoring, file_name)
else:
    matrix = nw_forward_correct(seq_first, seq_second, scoring)
    given_matrix_csv_maker(seq_first, seq_second, matrix, file_name)