import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\inputs\\ManhattanTour")
    filepaths = glob(dirpath + "\\data*.txt")
    for filepath in filepaths:
        n, m, down, right = ReadTestInputs_ManhattanTour(filepath)
        answer = ManhattanTourist(n, m, down, right)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(answer))

def ReadTestInputs_ManhattanTour(FilePath: str) -> tuple[int, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        n: Number of rows - 1.

        m: Number of columns - 1.

        down: matrix of down edge weights.

        right: matrix of right edge weights.
    """
    edges = {"down": [], "right":[]}
    key = "down"
    with open(FilePath) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if idx == 0:
                arr = line.split(" ")
                n, m = int(arr[0]), int(arr[1])
            else:
                if line == "-":
                    key = "right"
                else:
                    int_line = [int(x) for x in line.split(" ")]
                    edges[key].append(int_line)

    return n, m, edges["down"], edges["right"]

def ManhattanTourist(n: int, m: int, down: list[list[int]], right: list[list[int]]) -> int:
    """
    ManhattanTourist: Finds the maximum weight path on a square grid graph.

    Input:
        n: the number of rows on the graph - 1.

        m: the number of cols in the graph - 1.

        down: a matrix of downward edge weights.

        right: a matrix of right edge weights.
    
    Output:
        The weight of the maximum weight path.
    """
    s = [[0 for j in range(m+1)] for i in range(n+1)]
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j]+down[i-1][j], s[i][j-1]+right[i][j-1])
    return s[n][m]

# ManhattanTourist(n, m, Down, Right)
#     s0, 0 ← 0
#     for i ← 1 to n
#         si, 0 ← si-1, 0 + downi-1, 0
#     for j ← 1 to m
#         s0, j ← s0, j−1 + right0, j-1
#     for i ← 1 to n
#         for j ← 1 to m
#             si, j ← max{si - 1, j + downi-1, j, si, j - 1 + righti, j-1}
#     return sn, m

if __name__ == "__main__":
    main()