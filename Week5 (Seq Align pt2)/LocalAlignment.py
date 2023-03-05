import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\localalignment")
    filepaths = glob(os.path.join(dirpath, "data*.txt"))

    for path in filepaths:
        seq1, seq2= ReadTest_LocalAlignment(path)
        score, final_seq1, final_seq2 = LocalAlignment(seq1, seq2)
    
    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n")
        f.write(final_seq2)

def ReadTest_LocalAlignment(filepath) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    return seq1, seq2

def LocalAlignment(seq1:str, seq2: str) -> tuple[int, str, str]:
    """
    LocalAlignment: Finds the highest scoring local alignment between the two strings.

    Input:
        seq1: The first peptide sequence.
        seq2: The second peptide sequence

    Output:
        The score of the longest alignment, and the strings that represent the alignments.
    """
    indel = 5 # Penalty for indel.
    score_matrix = [[0 for b in range(len(seq2)+1)] for a in range(len(seq1)+1)]
    path_matrix = [[
        {"src": False, "dia": False, "hor": False, "ver": False} for b in range(len(seq2)+1)
        ] for a in range(len(seq1)+1)]
    
    top_score = 0
    ts_position = [0,0] #index of the top score
    for r in range(1, len(seq1)+1):
        for c in range(1, len(seq2)+1):

            seq1_score = score_matrix[r-1][c] - indel              # Score for taking letter from seq1.
            seq2_score = score_matrix[r][c-1] - indel             # Score for taking letter from seq2.
            both_score = score_matrix[r-1][c-1] + pam250[seq1[r-1]][seq2[c-1]]   # Score for taking letters from both seq1 and seq2.
            
            score_matrix[r][c] = max(0, seq1_score,seq2_score, both_score)
            if score_matrix[r][c] > top_score:
                top_score = score_matrix[r][c]
                ts_position = [r, c]

            if score_matrix[r][c] == 0:
                path_matrix[r][c]["src"] = True
            if score_matrix[r][c] == both_score:
                path_matrix[r][c]["dia"] = True
            if score_matrix[r][c] == seq2_score:
                path_matrix[r][c]["hor"] = True
            if score_matrix[r][c] == seq1_score:
                path_matrix[r][c]["ver"] = True

    score_matrix[-1][-1] = top_score
    
    final_seq1, final_seq2 = LocalBacktrack(seq1, seq2, path_matrix, ts_position)

    return score_matrix[-1][-1], final_seq1, final_seq2

def LocalBacktrack(seq1: str, seq2: str, path_matrix: list[list[dict[str, bool]]], idx: list[int, int]) -> tuple[str, str]:
    """
    LocalBacktrack: Returns the substrings with highest local score based on a path_matrix generated when looking for
    largest local alignment.

    input:
        seq1, seq2: The two sequences used to generate the path matrix.
        pathmatrix: The path matrix generated for local alignment.
        idx: The last position of the longest alignment in pathmatrix.
    
    output:
        The longest local alignment sequence from seq1 and seq2.
    """
    r = idx[0]
    c = idx[1]

    final_seq1 = ""
    final_seq2 = ""

    while (r != 0) or (c != 0):
        if path_matrix[r][c]["dia"] == True:
            r -= 1
            c -= 1
            final_seq1 = "".join([seq1[r], final_seq1])
            final_seq2 = "".join([seq2[c], final_seq2])
        elif path_matrix[r][c]["hor"] == True:
            c -= 1
            final_seq1 = "".join(["-", final_seq1])
            final_seq2 = "".join([seq2[c], final_seq2])
        elif path_matrix[r][c]["ver"] == True:
            r -= 1
            final_seq1 = "".join([seq1[r], final_seq1])
            final_seq2 = "".join(["-", final_seq2])
        elif path_matrix[r][c]["src"] == True:
            r = 0
            c = 0

    return final_seq1, final_seq2

pam250 = {
        "A":{"A":  2, "C": -2, "D":  0, "E":  0, "F": -3, 
            "G":  1, "H": -1, "I": -1, "K": -1, "L": -2, 
            "M": -1, "N":  0, "P":  1, "Q":  0, "R": -2, 
            "S":  1, "T":  1, "V":  0, "W": -6, "Y": -3
            },
        "C":{"A": -2, "C":  12, "D": -5, "E": -5, "F": -4, 
            "G": -3, "H": -3, "I": -2, "K": -5, "L": -6, 
            "M": -5, "N": -4, "P": -3, "Q": -5, "R": -4, 
            "S":  0, "T": -2, "V": -2, "W": -8, "Y":  0
            },
        "D":{"A":  0, "C": -5, "D":  4, "E":  3, "F": -6, 
            "G":  1, "H":  1, "I": -2, "K":  0, "L": -4, 
            "M": -3, "N":  2, "P": -1, "Q":  2, "R": -1, 
            "S":  0, "T":  0, "V": -2, "W": -7, "Y": -4
            },
        "E":{"A":  0, "C": -5, "D":  3, "E":  4, "F": -5, 
            "G":  0, "H":  1, "I": -2, "K":  0, "L": -3, 
            "M": -2, "N":  1, "P": -1, "Q":  2, "R": -1, 
            "S":  0, "T":  0, "V": -2, "W": -7, "Y": -4
            },
        "F":{"A": -3, "C": -4, "D": -6, "E": -5, "F":  9, 
            "G": -5, "H": -2, "I":  1, "K": -5, "L":  2, 
            "M":  0, "N": -3, "P": -5, "Q": -5, "R": -4, 
            "S": -3, "T": -3, "V": -1, "W":  0, "Y":  7
            },
        "G":{"A":  1, "C": -3, "D":  1, "E":  0, "F": -5, 
            "G":  5, "H": -2, "I": -3, "K": -2, "L": -4, 
            "M": -3, "N":  0, "P":  0, "Q": -1, "R": -3, 
            "S":  1, "T":  0, "V": -1, "W": -7, "Y": -5
            },
        "H":{"A": -1, "C": -3, "D":  1, "E":  1, "F": -2, 
            "G": -2, "H":  6, "I": -2, "K":  0, "L": -2, 
            "M": -2, "N":  2, "P":  0, "Q":  3, "R":  2, 
            "S": -1, "T": -1, "V": -2, "W": -3, "Y":  0
            },
        "I":{"A": -1, "C": -2, "D": -2, "E": -2, "F":  1, 
            "G": -3, "H": -2, "I":  5, "K": -2, "L":  2, 
            "M":  2, "N": -2, "P": -2, "Q": -2, "R": -2, 
            "S": -1, "T":  0, "V":  4, "W": -5, "Y": -1
            },
        "K":{"A": -1, "C": -5, "D":  0, "E":  0, "F": -5, 
            "G": -2, "H":  0, "I": -2, "K":  5, "L": -3, 
            "M":  0, "N":  1, "P": -1, "Q":  1, "R":  3, 
            "S":  0, "T":  0, "V": -2, "W": -3, "Y": -4
            },
        "L":{"A": -2, "C": -6, "D": -4, "E": -3, "F":  2, 
            "G": -4, "H": -2, "I":  2, "K": -3, "L":  6, 
            "M":  4, "N": -3, "P": -3, "Q": -2, "R": -3, 
            "S": -3, "T": -2, "V":  2, "W": -2, "Y": -1
            },
        "M":{"A": -1, "C": -5, "D": -3, "E": -2, "F":  0, 
            "G": -3, "H": -2, "I":  2, "K":  0, "L":  4, 
            "M":  6, "N": -2, "P": -2, "Q": -1, "R":  0, 
            "S": -2, "T": -1, "V":  2, "W": -4, "Y": -2
            },
        "N":{"A":  0, "C": -4, "D":  2, "E":  1, "F": -3, 
            "G":  0, "H":  2, "I": -2, "K":  1, "L": -3, 
            "M": -2, "N":  2, "P":  0, "Q":  1, "R":  0, 
            "S":  1, "T":  0, "V": -2, "W": -4, "Y": -2
            },
        "P":{"A":  1, "C": -3, "D": -1, "E": -1, "F": -5, 
            "G":  0, "H":  0, "I": -2, "K": -1, "L": -3, 
            "M": -2, "N":  0, "P":  6, "Q":  0, "R":  0, 
            "S":  1, "T":  0, "V": -1, "W": -6, "Y": -5
            },
        "Q":{"A":  0, "C": -5, "D":  2, "E":  2, "F": -5, 
            "G": -1, "H":  3, "I": -2, "K":  1, "L": -2, 
            "M": -1, "N":  1, "P":  0, "Q":  4, "R":  1, 
            "S": -1, "T": -1, "V": -2, "W": -5, "Y": -4
            },
        "R":{"A": -2, "C": -4, "D": -1, "E": -1, "F": -4, 
            "G": -3, "H":  2, "I": -2, "K":  3, "L": -3, 
            "M":  0, "N":  0, "P":  0, "Q":  1, "R":  6, 
            "S":  0, "T": -1, "V": -2, "W":  2, "Y": -4
            },
        "S":{"A":  1, "C":  0, "D":  0, "E":  0, "F": -3, 
            "G":  1, "H": -1, "I": -1, "K":  0, "L": -3, 
            "M": -2, "N":  1, "P":  1, "Q": -1, "R":  0, 
            "S":  2, "T":  1, "V": -1, "W": -2, "Y": -3
            },
        "T":{"A":  1, "C": -2, "D":  0, "E":  0, "F": -3, 
            "G":  0, "H": -1, "I":  0, "K":  0, "L": -2, 
            "M": -1, "N":  0, "P":  0, "Q": -1, "R": -1, 
            "S":  1, "T":  3, "V":  0, "W": -5, "Y": -3
            },
        "V":{"A":  0, "C": -2, "D": -2, "E": -2, "F": -1, 
            "G": -1, "H": -2, "I":  4, "K": -2, "L":  2, 
            "M":  2, "N": -2, "P": -1, "Q": -2, "R": -2, 
            "S": -1, "T":  0, "V":  4, "W": -6, "Y": -2
            },
        "W":{"A": -6, "C": -8, "D": -7, "E": -7, "F":  0, 
            "G": -7, "H": -3, "I": -5, "K": -3, "L": -2, 
            "M": -4, "N": -4, "P": -6, "Q": -5, "R":  2, 
            "S": -2, "T": -5, "V": -6, "W":  17, "Y":  0
            },
        "Y":{"A": -3, "C":  0, "D": -4, "E": -4, "F":  7, 
            "G": -5, "H":  0, "I": -1, "K": -4, "L": -1, 
            "M": -2, "N": -2, "P": -5, "Q": -4, "R": -4, 
            "S": -3, "T": -3, "V": -2, "W":  0, "Y":  10
            },
        }

if __name__ == "__main__":
    main()

