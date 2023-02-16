import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\fitalignment")
    filepaths = glob(os.path.join(dirpath, "dat*.txt"))

    for path in filepaths:
        seq1, seq2= ReadTest_FitAlignment(path)
        score, final_seq1, final_seq2 = FitAlignment(seq1, seq2)
    
    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n" + final_seq2)


def ReadTest_FitAlignment(filepath) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    return seq1, seq2

def FitAlignment(s: str, t: str) -> int:
    if len(s) >= len(t):
        v = s
        w = t
    else:
        v = t
        w = s

    score_matrix = [[0 for b in range(len(w)+1)] for a in range(len(v)+1)]
    path_matrix = [[
    {"dia": False, "hor": False, "ver": False} for b in range(len(w)+1)
    ] for a in range(len(v)+1)]

    for c in range(1, len(w)+1):
        score_matrix[0][c] = score_matrix[0][c-1] - 1
        path_matrix[0][c]["hor"] = True

    gap = 1

    curr_max = 0
    end_idx = len(w)-1
    for r in range(1, len(v)+1):
        for c in range(1, len(w)+1):

            score = blossum62[v[r-1]][w[c-1]]

            diag = score_matrix[r-1][c-1] + score
            horz = score_matrix[r][c-1] - gap
            vert = score_matrix[r-1][c] - gap
            
            score_matrix[r][c] = max(diag, vert, horz)
            if score_matrix[r][c] == diag:
                path_matrix[r][c]["dia"] = True
            elif score_matrix[r][c] == horz:
                path_matrix[r][c]["hor"] = True
            elif score_matrix[r][c] == vert:
                path_matrix[r][c]["ver"] = True
        if score_matrix[r][len(w)] > curr_max:
            curr_max = score_matrix[r][len(w)]
            end_idx = r

    final_seq1, final_seq2 = FitBacktrack(v, w, path_matrix, end_idx)
    return score_matrix[end_idx][-1], final_seq1, final_seq2

def FitBacktrack(seq1: str, seq2: str, path_matrix: list[list[dict[str, bool]]], idx: int) -> tuple[str, str]:
    """
    FitBacktrack: Returns the substrings with highest fit alignment of two strings based on a generated path matrix.

    input:
        seq1, seq2: The two sequences used to generate the path matrix.
        pathmatrix: The path matrix generated for local alignment.
        idx: The last position of the longest alignment in pathmatrix.
    
    output:
        The best fit alignment sequence from seq1 and seq2.
    """
    r = idx
    c = len(seq2)

    final_seq1 = ""
    final_seq2 = ""

    while (c != 0):
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

    return final_seq1, final_seq2

blossum62 = {
"A":{"A":  4, "C":  0, "D": -2, "E": -1, "F": -2, 
	 "G":  0, "H": -2, "I": -1, "K": -1, "L": -1, 
	 "M": -1, "N": -2, "P": -1, "Q": -1, "R": -1, 
	 "S":  1, "T":  0, "V":  0, "W": -3, "Y": -2
	 },
"C":{"A":  0, "C":  9, "D": -3, "E": -4, "F": -2, 
	 "G": -3, "H": -3, "I": -1, "K": -3, "L": -1, 
	 "M": -1, "N": -3, "P": -3, "Q": -3, "R": -3, 
	 "S": -1, "T": -1, "V": -1, "W": -2, "Y": -2
	 },
"D":{"A": -2, "C": -3, "D":  6, "E":  2, "F": -3, 
	 "G": -1, "H": -1, "I": -3, "K": -1, "L": -4, 
	 "M": -3, "N":  1, "P": -1, "Q":  0, "R": -2, 
	 "S":  0, "T": -1, "V": -3, "W": -4, "Y": -3
	 },
"E":{"A": -1, "C": -4, "D":  2, "E":  5, "F": -3, 
	 "G": -2, "H":  0, "I": -3, "K":  1, "L": -3, 
	 "M": -2, "N":  0, "P": -1, "Q":  2, "R":  0, 
	 "S":  0, "T": -1, "V": -2, "W": -3, "Y": -2
	 },
"F":{"A": -2, "C": -2, "D": -3, "E": -3, "F":  6, 
	 "G": -3, "H": -1, "I":  0, "K": -3, "L":  0, 
	 "M":  0, "N": -3, "P": -4, "Q": -3, "R": -3, 
	 "S": -2, "T": -2, "V": -1, "W":  1, "Y":  3
	 },
"G":{"A":  0, "C": -3, "D": -1, "E": -2, "F": -3, 
	 "G":  6, "H": -2, "I": -4, "K": -2, "L": -4, 
	 "M": -3, "N":  0, "P": -2, "Q": -2, "R": -2, 
	 "S":  0, "T": -2, "V": -3, "W": -2, "Y": -3
	 },
"H":{"A": -2, "C": -3, "D": -1, "E":  0, "F": -1, 
	 "G": -2, "H":  8, "I": -3, "K": -1, "L": -3, 
	 "M": -2, "N":  1, "P": -2, "Q":  0, "R":  0, 
	 "S": -1, "T": -2, "V": -3, "W": -2, "Y":  2
	 },
"I":{"A": -1, "C": -1, "D": -3, "E": -3, "F":  0, 
	 "G": -4, "H": -3, "I":  4, "K": -3, "L":  2, 
	 "M":  1, "N": -3, "P": -3, "Q": -3, "R": -3, 
	 "S": -2, "T": -1, "V":  3, "W": -3, "Y": -1
	 },
"K":{"A": -1, "C": -3, "D": -1, "E":  1, "F": -3, 
	 "G": -2, "H": -1, "I": -3, "K":  5, "L": -2, 
	 "M": -1, "N":  0, "P": -1, "Q":  1, "R":  2, 
	 "S":  0, "T": -1, "V": -2, "W": -3, "Y": -2
	 },
"L":{"A": -1, "C": -1, "D": -4, "E": -3, "F":  0, 
	 "G": -4, "H": -3, "I":  2, "K": -2, "L":  4, 
	 "M":  2, "N": -3, "P": -3, "Q": -2, "R": -2, 
	 "S": -2, "T": -1, "V":  1, "W": -2, "Y": -1
	 },
"M":{"A": -1, "C": -1, "D": -3, "E": -2, "F":  0, 
	 "G": -3, "H": -2, "I":  1, "K": -1, "L":  2, 
	 "M":  5, "N": -2, "P": -2, "Q":  0, "R": -1, 
	 "S": -1, "T": -1, "V":  1, "W": -1, "Y": -1
	 },
"N":{"A": -2, "C": -3, "D":  1, "E":  0, "F": -3, 
	 "G":  0, "H":  1, "I": -3, "K":  0, "L": -3, 
	 "M": -2, "N":  6, "P": -2, "Q":  0, "R":  0, 
	 "S":  1, "T":  0, "V": -3, "W": -4, "Y": -2
	 },
"P":{"A": -1, "C": -3, "D": -1, "E": -1, "F": -4, 
	 "G": -2, "H": -2, "I": -3, "K": -1, "L": -3, 
	 "M": -2, "N": -2, "P":  7, "Q": -1, "R": -2, 
	 "S": -1, "T": -1, "V": -2, "W": -4, "Y": -3
	 },
"Q":{"A": -1, "C": -3, "D":  0, "E":  2, "F": -3, 
	 "G": -2, "H":  0, "I": -3, "K":  1, "L": -2, 
	 "M":  0, "N":  0, "P": -1, "Q":  5, "R":  1, 
	 "S":  0, "T": -1, "V": -2, "W": -2, "Y": -1
	 },
"R":{"A": -1, "C": -3, "D": -2, "E":  0, "F": -3, 
	 "G": -2, "H":  0, "I": -3, "K":  2, "L": -2, 
	 "M": -1, "N":  0, "P": -2, "Q":  1, "R":  5, 
	 "S": -1, "T": -1, "V": -3, "W": -3, "Y": -2
	 },
"S":{"A":  1, "C": -1, "D":  0, "E":  0, "F": -2, 
	 "G":  0, "H": -1, "I": -2, "K":  0, "L": -2, 
	 "M": -1, "N":  1, "P": -1, "Q":  0, "R": -1, 
	 "S":  4, "T":  1, "V": -2, "W": -3, "Y": -2
	 },
"T":{"A":  0, "C": -1, "D": -1, "E": -1, "F": -2, 
	 "G": -2, "H": -2, "I": -1, "K": -1, "L": -1, 
	 "M": -1, "N":  0, "P": -1, "Q": -1, "R": -1, 
	 "S":  1, "T":  5, "V":  0, "W": -2, "Y": -2
	 },
"V":{"A":  0, "C": -1, "D": -3, "E": -2, "F": -1, 
	 "G": -3, "H": -3, "I":  3, "K": -2, "L":  1, 
	 "M":  1, "N": -3, "P": -2, "Q": -2, "R": -3, 
	 "S": -2, "T":  0, "V":  4, "W": -3, "Y": -1
	 },
"W":{"A": -3, "C": -2, "D": -4, "E": -3, "F":  1, 
	 "G": -2, "H": -2, "I": -3, "K": -3, "L": -2, 
	 "M": -1, "N": -4, "P": -4, "Q": -2, "R": -3, 
	 "S": -3, "T": -2, "V": -3, "W":  11, "Y":  2
	 },
"Y":{"A": -2, "C": -2, "D": -3, "E": -2, "F":  3, 
	 "G": -3, "H":  2, "I": -1, "K": -2, "L": -1, 
	 "M": -1, "N": -2, "P": -3, "Q": -1, "R": -2, 
	 "S": -2, "T": -2, "V": -1, "W":  2, "Y":  7
	 },
}

if __name__ == "__main__":
    main()