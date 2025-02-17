# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    gap_open = -10
    gap_extend = -1

    #initialize NW
    test = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=gap_open, gap_extend=gap_extend,)

    #calculate alignment score and store aligned sequences for each speicies
    score_hs_gg, hs_align, gg_align = test.align(hs_seq, gg_seq)
    score_hs_mm, hs, mm_align = test.align(hs_seq, mm_seq)
    score_hs_br, hs_align, br_align = test.align(hs_seq, br_seq)
    score_hs_tt, hs_align, tt_align = test.align(hs_seq, tt_seq)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sequence_list = [(score_hs_gg, 'Gallus_gallus'), (score_hs_mm, 'Mus_musculus'), (score_hs_br, 'Balaeniceps_rex'), (score_hs_tt,'Tursiops-truncatus')]

    #sort the list of species by their alignment scores
    sorted_sequence = sorted(sequence_list, reverse = True)

    print("Original nameList:", [x[1] for x in sequence_list])
    print("Sorted nameList:", [x[1] for x in sorted_sequence])

    

if __name__ == "__main__":
    main()
