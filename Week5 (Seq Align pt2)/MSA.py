import os, math
import numpy as np
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\msa")
    filepaths = glob(dirpath + "\\\data*.txt")
    for file in filepaths:
        seq1, seq2, seq3 = ReadTest_MSA(file)
        score, final_seq1, final_seq2, final_seq3 = MultipleAlignment(seq1, seq2, seq3)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, "w") as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n")
        f.write(final_seq2 + "\n")
        f.write(final_seq3)

def ReadTest_MSA(filepath:str) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
        seq3 = f.readline().strip()
    return seq1, seq2, seq3

def MultipleAlignment(s1: str, s2: str, s3: str) -> tuple[int, str, str, str]:
    """
    MultipleAlignment: Aligns 3 sequences using dynamic programming.
    Input:
        s1: sequence 1
        s2: sequence 2
        s3: sequence 3
    Tuple:
        The score and the alignment of the 3 strings that create the optimal score.
    """
    score_matrix = [[[0 for z in range(len(s3)+1)] for y in range(len(s2)+1)] for x in range(len(s1)+1)]
    path_matrix = [[[{"111": False, "110": False, "101": False, "011": False,
                      "100": False, "010": False, "001": False} for z in range(len(s3)+1)]
                       for y in range(len(s2)+1)] for x in range(len(s1)+1)]

    for x in range(1, len(s1)+1):
        path_matrix[x][0][0]["100"] = True
        for y in range(1, len(s2)+1):
            path_matrix[x][y][0]["110"] = True
        for z in range(1, len(s3)+1):
            path_matrix[x][0][z]["101"] = True
    
    for y in range(1, len(s2)+1):
        path_matrix[0][y][0]["010"] = True
        for z in range(1, len(s3)+1):
            path_matrix[0][y][z]["011"] = True

    for z in range(1, len(s3)+1):
        path_matrix[0][0][z]["001"] = True

    neighbors = ["111", "110", "101", "011", "100", "010", "001"]
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            for k in range(1, len(s3)+1):
                scores = {}
                for l in neighbors:
                    scores[l] = score_matrix[i-int(l[0])][j-int(l[1])][k-int(l[2])]
                    if (l == "111") and (s1[i - 1] == s2[j - 1]) and (s1[i - 1] == s3[k - 1]):
                        scores[l] += 1

                score_matrix[i][j][k] = max(scores.values())
                for key in scores:
                    if scores[key] == score_matrix[i][j][k]:
                        path_matrix[i][j][k][key] = True
    
    final_seq1, final_seq2, final_seq3 = Backtrack(s1, s2, s3, path_matrix)
    return score_matrix[i][j][k], final_seq1, final_seq2, final_seq3
                
def Backtrack(seq1: str, seq2: str, seq3: str, path_matrix: list[list[list[dict[str, bool]]]]) -> tuple[str, str, str]:
    final_seq1 = ""
    final_seq2 = ""
    final_seq3 = ""

    x = len(seq1)
    y = len(seq2)
    z = len(seq3)

    neighbors = ["111", "110", "101", "011", "100", "010", "001"]
    while (x > 0) or (y > 0) or (z > 0):
        for neighbor in neighbors:
            if path_matrix[x][y][z][neighbor] == True:
                x -= int(neighbor[0])
                y -= int(neighbor[1])
                z -= int(neighbor[2])
                nxt_char1 = "-"
                nxt_char2 = "-"
                nxt_char3 = "-"
                if int(neighbor[0]):
                    nxt_char1 = seq1[x]
                if int(neighbor[1]):
                    nxt_char2 = seq2[y]
                if int(neighbor[2]):
                    nxt_char3 = seq3[z]
                break

        final_seq1 = "".join([nxt_char1, final_seq1])
        final_seq2 = "".join([nxt_char2, final_seq2])
        final_seq3 = "".join([nxt_char3, final_seq3])

    return final_seq1, final_seq2, final_seq3


if __name__ == "__main__":
    main()
