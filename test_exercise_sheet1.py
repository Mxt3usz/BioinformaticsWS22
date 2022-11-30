from random import randint
from exercise_sheet1 import *
from helpers import translation_table


def check_none(*args):
    if None in args:
        print("you have not filled in all the fields")
        raise ValueError


def check_zero_in_dict(d):
    val = d.values()
    if 0 in val:
        print("you have not filled in all the fields")
        raise ValueError


def correct_transcription(dna_seq):
    return dna_seq.replace("T", "U")


def correct_translation(rna_seq):
    return "".join([translation_table[rna_seq[3*i:(3*i+3)]]
                    for i in range(len(rna_seq)//3)])


def correct_central_dogma(dna_seq):
    return correct_translation(correct_transcription(dna_seq))


def dna_generator():
    return "".join(["AGCT"[randint(0, 3)] for _ in range(randint(42, 57))])


def rna_generator():
    rna = "AUG"
    for i in range(randint(30, 40)):
        codon = "".join(["AGCU"[randint(0, 3)] for _ in range(3)])
        if translation_table[codon] != "*":
            rna += codon
    rna += ["UAA", "UAG", "UGA"][randint(0, 2)]
    return rna



def test_exercise_4a():
    for _ in range(10):
        dna = dna_generator()
        correct_rna = correct_transcription(dna)
        provided_dna = transcription(dna)
        if correct_rna != provided_dna:
            print(f"Your solution has failed\n"
                  f"Input DNA sequence {dna}\n"
                  f"Your RNA sequence {provided_dna}\n"
                  f"Correct RNA sequence {correct_rna}")

            assert correct_rna == provided_dna


def test_exercise_4b():
    for _ in range(10):
        rna = rna_generator()
        correct_protein = correct_translation(rna)
        provided_protein = translation(rna)
        if correct_protein != provided_protein:
            print(f"Your solution has failed\n"
                  f"Input RNA sequence {rna}\n"
                  f"Your protein sequence {provided_protein}\n"
                  f"Correct protein sequence {correct_protein}")

            assert correct_protein == provided_protein


def test_exercise_4c():
    for _ in range(10):
        rna = rna_generator()
        dna = rna.replace("U", "T")
        correct_protein = correct_central_dogma(dna)
        provided_protein = central_dogma(dna)
        if correct_protein != provided_protein:
            print(f"Your solution has failed\n"
                  f"Input DNA sequence {dna}\n"
                  f"Your protein sequence {provided_protein}\n"
                  f"Correct protein sequence {correct_protein}")

            assert correct_protein == provided_protein
