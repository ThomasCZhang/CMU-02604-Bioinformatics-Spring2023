import os
from glob import glob

def main():
    # path to directory containing test files
    dirpath = os.path.join(os.path.dirname(__file__), "files\\inputs\\lcs")
    filepaths = glob(dirpath + "\\input_0.txt")

    for path in filepaths:
        v, w = ReadLCSData(path)
        backtrack = LCSBacktrack(v, w)
        # print(len(backtrack), len(backtrack[0]))
        answer = OutputLCS(backtrack, v, len(v), len(w))

    answer_path = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answer_path, 'w') as f:
        f.write(answer)


def ReadLCSData(path: str) -> tuple[str, str]:
    """
    ReadLCSData: Reads in data for longest common string problem.

    Input:
        path: path to the data file.

    Output:
        v, w: the two strings used as input for longest common string.
    """
    with open(path) as f:
        v = f.readline()
        w = f.readline()
        v = v.strip()
        w = w.strip()
    
    return v, w

def LCSBacktrack(v: str, w: str):
    """
    LCSBacktrack: Finds the alignment of two strings that will create the longest common string between the two
    strings.

    Input:
        v: first string.

        w: second string.
    
    Output:
        backtrack: a matrix that provides info on how to backtrack from a alignment between v and w.
    """
    score_matrix = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    backtrack = [[[0, 0, 0] for j in range(len(w)+1)] for i in range(len(v)+1)]
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            score_matrix[i][j] = max(score_matrix[i-1][j],
                                     score_matrix[i][j-1],
                                     score_matrix[i-1][j-1]+match)
            if score_matrix[i][j] == score_matrix[i-1][j]:
                backtrack[i][j][0] = 1
            elif score_matrix[i][j] == score_matrix[i][j-1]:
                backtrack[i][j][1] = 1
            elif score_matrix[i][j] == score_matrix[i-1][j-1]+match:
                backtrack[i][j][2] = 1
    return backtrack

def OutputLCS(backtrack: list[list[list[int]]], v: str, i: int, j: int) -> str:
    """
    OutputLCS: Returns the path taken to reach element i, j in the matrix "backtrack".
    """
    traceback_path = ""
    while i != 0 and j!= 0:
        if backtrack[i][j][0] == 1:
            i -= 1
        elif backtrack[i][j][1] == 1:
            j -= 1
        elif backtrack[i][j][2] == 1:
            traceback_path = v[i-1] + traceback_path
            i -= 1
            j -= 1
    
    return traceback_path

    # if i== 0 or j == 0:
    #     return ""
    # if backtrack[i][j][0] == 1:
    #     return OutputLCS(backtrack, v, i-1, j)
    # elif backtrack[i][j][1] == 1:
    #     return OutputLCS(backtrack, v, i, j-1)
    # elif backtrack[i][j][2] == 1:
    #     return OutputLCS(backtrack, v, i-1, j-1) + v[i-1]

if __name__ == "__main__":
    main()