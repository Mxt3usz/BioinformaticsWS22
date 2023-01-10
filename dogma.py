

def transcription(dna_string):
    """
    Implement the transcription function which takes a DNA sequence as an input and returns the corresponding
    RNA sequence.
    Assume that the provided DNA sequence will be an uppercase string containing only the correct characters
    (from the domain "AGCT"). The function output RNA sequence also has to be in the upper case.
    """
    rna_string = ""
    for base in dna_string:
        if base == "T":
            rna_string += "U"
        else:
            rna_string += base
    return rna_string


def translation(rna_string):
    """
     Implement the translation function which takes an mRNA sequence and produces the corresponding protein sequence.
     Assume that the provided RNA sequence will be an uppercase string containing only the correct characters
     (from the domain "AGCU"). Assume that the sequence starts with the start codon (AUG) and ends with one of the
     end codons. It does not have any end codons in the middle.The number of characters in the input RNA sequence is
     a multiple of three.
     The function output sequence also has to be in the upper case protein sequence. To make it easier for you,
     we provided you with the translation table which is represented as a python dictionary. You can translate a
     single codon by simply using it:

     translation_table["AUG"]
    """
    translation_table = dict(AUA='I', AUC='I', AUU='I', AUG='M', ACA='T', ACC='T', ACG='T', ACU='T', AAC='N', AAU='N',
                         AAA='K', AAG='K', AGC='S', AGU='S', AGA='R', AGG='R', CUA='L', CUC='L', CUG='L', CUU='L',
                         CCA='P', CCC='P', CCG='P', CCU='P', CAC='H', CAU='H', CAA='Q', CAG='Q', CGA='R', CGC='R',
                         CGG='R', CGU='R', GUA='V', GUC='V', GUG='V', GUU='V', GCA='A', GCC='A', GCG='A', GCU='A',
                         GAC='D', GAU='D', GAA='E', GAG='E', GGA='G', GGC='G', GGG='G', GGU='G', UCA='S', UCC='S',
                         UCG='S', UCU='S', UUC='F', UUU='F', UUA='L', UUG='L', UAC='Y', UAU='Y', UAA='*', UAG='*',
                         UGC='C', UGU='C', UGA='*', UGG='W')
    protein = ""
    triplet = ""
    for base in rna_string:
        triplet += base
        if len(triplet) == 3:
            protein += translation_table[triplet]
            triplet = ""
    return protein


def central_dogma(dna_seq):
    """
    Combine the functions you implemented for the tasks a) and b) and obtain the central dogma function which
    gives you the protein sequence with a given DNA.
    """
    return translation(transcription(dna_seq))
