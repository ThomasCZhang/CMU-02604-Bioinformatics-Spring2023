import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\editdistance")
    filepaths = glob(os.path.join(dirpath, "data*.txt"))

    for path in filepaths:
        seq1, seq2= ReadTest_EditDist(path)
        score = EditDistance(seq1, seq2)
    
    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(score))


def ReadTest_EditDist(filepath) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    return seq1, seq2

def EditDistance(s: str, t: str) -> int:
    score_matrix = [[0 for b in range(len(t)+1)] for a in range(len(s)+1)]

    for r in range(1, len(s)+1):
        score_matrix[r][0] = score_matrix[r-1][0] + 1     

    for c in range(1, len(t)+1):
        score_matrix[0][c] = score_matrix[0][c-1] + 1

    for r in range(1, len(s)+1):
        for c in range(1, len(t)+1):
            score = 1
            if s[r-1] == t[c-1]:
                score = 0
            diag = score_matrix[r-1][c-1] + score
            vert = score_matrix[r][c-1] + 1
            horz = score_matrix[r-1][c] + 1
            
            score_matrix[r][c] = min(diag, vert, horz)
    
    return score_matrix[-1][-1]


if __name__ == "__main__":
    main()