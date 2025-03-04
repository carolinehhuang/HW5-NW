# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        #check that the parameters have the correct format and are provided
        if not isinstance(seqA, str) or not isinstance(seqB, str):
            raise TypeError("seqA and seqB must be strings")

        if not seqA and not seqB:
            raise ValueError("Sequences are not provided")


        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores and backtracing
        rows, cols = len(seqA) + 1, len(seqB) + 1
        self._align_matrix = np.zeros([rows,cols]) #to store alignment scores
        self._back = np.full((rows,cols), fill_value = -1) #to store directions

        #define the alignment scores and backtracing values for the first row and column
        # alignment scores are just the gap penalties for the first row and column
        for i in range(1, rows):
            self._align_matrix[i][0] = self.gap_open + i * self.gap_extend
            self._back[i][0] = 1 #for the backtracing matrix, define 1 as "up", meaning that the current alignment score was derived from the value stored in the row above it

        for j in range(1, cols):
            self._align_matrix[0][j] = self.gap_open + j * self.gap_extend
            self._back[0][j] = 2 #for the bt matrix, define 2 as "left", meaning that the current alignment score was derived from the value stored in the column to its left

        for i in range(1, rows): #for every other row and column in the matrix...
            for j in range(1, cols):
                match = self._align_matrix[i-1][j-1] + self.sub_dict.get((seqA[i-1], seqB[j-1])) #match score is the substitution score plus the previous alignment score before
                empty = self._align_matrix[i-1][j] + (self.gap_extend if self._back[i - 1][j] == 1 else self.gap_extend + self.gap_open) #score in circumstance where seqB has a deletion and we need to insert a space in seqB for it to align with seqA; if seqB is just continuing an existing deletion, add gap_extension, but if  creating a new gap, add gap extension+gap open score
                insert = self._align_matrix[i][j-1] + (self.gap_extend if self._back[i][j - 1] == 2 else self.gap_extend + self.gap_open) #score in circumstance where seqB has an insertion and seqA needs an insertion to align with seqB

                best = max(match, empty, insert) #find the maximum value of 3 scenarios
                self._align_matrix[i][j] = best #assign best score to alignment matrix

                #fill backtrace matrix with direction depending on the beset alignment value
                if best == match:
                    self._back[i][j] = 0 #represents diagonal align[i-1,j-1]
                elif best == empty:
                    self._back[i][j] = 1 #represents up align[i-1, j]
                else:
                    self._back[i][j] = 2 #represents left align[i, j-1]

        self.alignment_score = self._align_matrix[-1,-1] #alignment score is the bottom right corner
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        rows, cols = len(self._seqA), len(self._seqB)

        while rows > 0 or cols > 0 and self._back[rows][cols]!= -1:

            if rows > 0 and cols > 0 and self._back[rows][cols] == 0: #if the best match was diagonal
                self.seqA_align = self._seqA[rows-1] + self.seqA_align #add the current letter of seqA to the aligned sequence
                self.seqB_align = self._seqB[cols-1] + self.seqB_align  #add current letter of seq B to the aligned sequence

                rows -= 1
                cols -= 1
            elif rows > 0 and self._back[rows][cols] == 1:
                self.seqA_align = self._seqA[rows-1] + self.seqA_align
                self.seqB_align = "-" + self.seqB_align #add a "-" in seqB to represent the deletion
                # move backward to the score it was derived from, the up value
                rows -= 1
            elif cols > 0:
                self.seqB_align = self._seqB[rows-1] + self.seqB_align
                self.seqA_align = "-" + self.seqA_align #add a "-" in seqA to represent the deletion
                # move backward to the score it was derived from, the left value
                cols -= 1

        return self.alignment_score, self.seqA_align, self.seqB_align


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
