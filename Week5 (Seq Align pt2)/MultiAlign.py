import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\MultiAlign")
    filepaths = glob(os.path.join(dirpath, "input*.txt"))

    for path in filepaths:
        seqs= ReadTest_MultiAlign(path)
        score, final_seq1, final_seq2 = MultiAlign(reward, mismatch, indel, seq1, seq2)
    
    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n")
        f.write(final_seq2)

    

def ReadTest_MultiAlign(filepath) -> list[str]:
    seqs = []
    with open(filepath) as f:
        for line in f:
            seqs.append(line.strip())
    return seqs

def Backtrack(seq1: str, seq2: str, path_matrix: list[list[int, int]]) -> tuple[str, str]:
    """
    Returns the "path" taken by Needleman Wunsch algorithm for aligning two strings.

    Input:
        seq1, seq2: The two strings being aligned

        path_matrix: A 3 dimensional numpy ndarray that contains info on the path taken by Needleman Wunsh
            The third dimension holds info on whether a letter was taken from seq1, seq2 or from both.
    
    Output:
        Two strings corresponding to the input sequences. "-" represents a gap.
    """
    final_seq1 = ""
    final_seq2 = ""
    
    row = len(seq1)
    col = len(seq2)
    
    while (row != 0) or (col != 0):
       
        if path_matrix[row][col][2] == 1:
            row -= 1
            col -= 1
            final_seq1 = seq1[row] + final_seq1
            final_seq2 = seq2[col] + final_seq2
        elif path_matrix[row][col][1] == 1:
            col -= 1
            final_seq1 = "-" + final_seq1
            final_seq2 = seq2[col] + final_seq2
        elif path_matrix[row][col][0] == 1:
            row -= 1
            final_seq1 = seq1[row] + final_seq1
            final_seq2 = "-" + final_seq2
    
    return final_seq1, final_seq2

def MultiAlign(char_match, char_mismatch, gap, seq1: str, seq2: str) -> tuple[int, str, str]:
    """
    GlobalAlignment: Finds the best global alignment of two strings.

    Input:
        char_match: the score when character's match.

        char_mismatch: the penalty when character's mismatch.

        gap: the penalty when there is a gap.

        seq1: the first string.

        seq2: the second string.

    Output:
        The score, and the alignment of the two strings. "-" represents a gap.
    """

    # Matrix that stores the scores of the subproblems that have been solved
    score_matrix = [[0 for i in range(len(seq2)+1)] for j in range(len(seq1) + 1)] 
    # Matrix that scores the "path" has been taken so far.
    path_matrix = [[[0, 0, 0] for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]
        # [i,j,0] take from str1,
        # [i,j,1] take from str2,
        # [i,j,2] take from both strings.

    for i in range(1, len(seq1)+1):
        score_matrix[i][0] = score_matrix[i-1][0] - gap
        path_matrix[i][0][0] = 1
    
    for i in range(1, len(seq2)+1):
        score_matrix[0][i] = score_matrix[0][i-1] - gap
        path_matrix[0][i][1] = 1

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            match_score = -char_mismatch
            if seq1[i-1] == seq2[j-1]:
                match_score = char_match

            seq1_score = score_matrix[i-1][j] - gap              # Score for taking letter from seq1.
            seq2_score = score_matrix[i][j-1] - gap             # Score for taking letter from seq2.
            both_score = score_matrix[i-1][j-1] + match_score   # Score for taking letters from both seq1 and seq2.

            score_matrix[i][j] = max(seq1_score, seq2_score, both_score)
            if score_matrix[i][j] == seq1_score:
                path_matrix[i][j][0] = 1
            if score_matrix[i][j] == seq2_score:
                path_matrix[i][j][1] = 1
            if score_matrix[i][j] == both_score:
                path_matrix[i][j][2] = 1

    final_seq1, final_seq2 = Backtrack(seq1, seq2, path_matrix)
    return (score_matrix[len(seq1)][len(seq2)], final_seq1, final_seq2)

if __name__ == "__main__":
    main()