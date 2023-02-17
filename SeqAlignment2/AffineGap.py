import os, math
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\affine")
    filepaths = glob(dirpath + "\\data*.txt")
    for file in filepaths:
        match_reward, mismatch_penalty, gap_open, gap_extend, seq1, seq2 = ReadTest_Affine(file)
        score, final_seq1, final_seq2 = AffineAlignment(match_reward, mismatch_penalty, gap_open, gap_extend, seq1, seq2)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, "w") as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n")
        f.write(final_seq2 + "\n")

def ReadTest_Affine(filepath:str) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        scores = [int(x) for x in f.readline().strip().split()]
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    return scores[0], scores[1] ,scores[2], scores[3], seq1, seq2

def AffineAlignment(match_reward: int, mismatch_penalty: int, gap_open: int, gap_extend: int, s: str, t: str) -> tuple[int, str, str]:
    """
    AffineGap: Finds the best alignment between two strings using the affine gap penalty algorithm.
    
    Input:
        match_reward: reward for matching chars
        mismatch_penalty: penalty for mismatching chars
        gap_open: penalty for starting a gap
        gap_extend: penalty for extending a gap.
    
    Output:
        The score of the alignment and the aligned sequences.
    """
    lower = [[-math.inf for j in range(len(t)+1)] for i in range(len(s)+1)]
    upper = [[-math.inf for j in range(len(t)+1)] for i in range(len(s)+1)]
    middle = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]

    lower_path = [[{"extend": False, "start": False} for j in range(len(t)+1)] for i in range(len(s)+1)]
    upper_path = [[{"extend": False, "start": False} for j in range(len(t)+1)] for i in range(len(s)+1)]
    middle_path = [[{"diag": False, "lower": False, "upper": False} for j in range(len(t)+1)] for i in range(len(s)+1)]

    for i in range(1, len(s)+1):
        if i == 1:
            lower_path[i][0]["start"] = True
            lower[i][0] = middle[i-1][0] - gap_open
        else:
            lower_path[i][0]["extend"] = True
            lower[i][0] = lower[i-1][0] - gap_extend
        middle_path[i][0]["lower"] = True
        middle[i][0] = lower[i][0]

    for j in range(1, len(t)+1):
        if j == 1:
            upper_path[0][j]["start"] = True
            upper[0][j] = middle[0][j-1] - gap_open
        else:
            upper_path[0][j]["extend"] = True
            upper[0][j] = upper[0][j-1] - gap_extend
        middle_path[0][j]["upper"] = True
        middle[0][j] = upper[0][j]

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            extend_gap_score = lower[i-1][j]-gap_extend
            start_gap_score = middle[i-1][j]-gap_open
            lower[i][j] = max(extend_gap_score, start_gap_score)

            if lower[i][j] == start_gap_score:
                lower_path[i][j]["start"] = True
            if lower[i][j] == extend_gap_score:
                lower_path[i][j]["extend"] = True

            extend_gap_score = upper[i][j-1]-gap_extend
            start_gap_score = middle[i][j-1]-gap_open
            upper[i][j] = max(extend_gap_score, start_gap_score)

            if upper[i][j] == start_gap_score:
                upper_path[i][j]["start"] = True
            if upper[i][j] == extend_gap_score:
                upper_path[i][j]["extend"] = True   

            char_score = -mismatch_penalty
            if s[i-1] == t[j-1]:
                char_score = match_reward
            close_lower = lower[i][j]
            close_upper = upper[i][j]
            diag_score = middle[i-1][j-1]+char_score

            middle[i][j] = max(close_lower, close_upper, diag_score)
            if middle[i][j] == diag_score:
                middle_path[i][j]["diag"] = True
            if middle[i][j] == close_lower:
                middle_path[i][j]["lower"] = True
            if middle[i][j] == close_upper:
                middle_path[i][j]["upper"] = True

    three_layer_path = [lower_path, middle_path, upper_path]
    final_seq1, final_seq2 = Backtrack(s, t, three_layer_path)

    return middle[i][j], final_seq1, final_seq2

def Backtrack(seq1:str, seq2: str, path_matrix: list[list[list[dict[str, bool]]]]) -> tuple[str, str]:
    """
    Backtrack: Backtracks through the path matrix of the affine gap problem. Returns the stings corresponding to the 
    path.

    Input:
        seq1, seq2: The two sequences.
    """
    i = len(seq1)
    j = len(seq2)
    layer = 1 

    final_seq1 = ""
    final_seq2 = ""    
    while (i > 0) or (j > 0):
        seq1_addition = ""
        seq2_addition = ""
        if layer == 1:
            if path_matrix[layer][i][j]["diag"] == True:
                i -= 1
                j -= 1
                seq1_addition = seq1[i]
                seq2_addition = seq2[j]
            elif path_matrix[layer][i][j]["lower"] == True:
                layer = 0
            elif path_matrix[layer][i][j]["upper"] == True:
                layer = 2    
        elif layer == 0:
            if path_matrix[layer][i][j]["start"] == True:
                i -= 1
                layer = 1
                seq1_addition = seq1[i]
                seq2_addition = "-"
            elif path_matrix[layer][i][j]["extend"] == True:
                i -= 1
                seq1_addition = seq1[i]
                seq2_addition = "-"
        elif layer == 2:
            if path_matrix[layer][i][j]["start"] == True:
                j -= 1
                layer = 1
                seq1_addition = "-"
                seq2_addition = seq2[j]
            elif path_matrix[layer][i][j]["extend"] == True:
                j -= 1
                seq1_addition = "-"
                seq2_addition = seq2[j]
        final_seq1 = "".join([seq1_addition, final_seq1])
        final_seq2 = "".join([seq2_addition, final_seq2])

    return final_seq1, final_seq2

if __name__ == "__main__":
    main()