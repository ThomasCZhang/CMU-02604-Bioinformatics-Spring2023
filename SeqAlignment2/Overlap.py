import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "file\\overlap")
    filepaths = glob(dirpath + "\\input*.txt")
    for file in filepaths:
        match_reward, mismatch_penalty, indel_penalty, seq1, seq2 = ReadTest_Overlap(file)
        score, final_seq1, final_seq2 = OverlapAlignment(match_reward, mismatch_penalty, indel_penalty, seq1, seq2)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, "w") as f:
        f.write(str(score) + "\n")
        f.write(final_seq1 + "\n")
        f.write(final_seq2 + "\n")

def ReadTest_Overlap(filepath:str) -> tuple[int, int, int, str, str]:
    with open(filepath) as f:
        scores = [int(x) for x in f.readline().strip().split()]
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    return scores[0], scores[1] ,scores[2], seq1, seq2


def OverlapAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                    s: str, t: str) -> tuple[int, str, str]:
    
    score_matrix = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    path_matrix = [[{"dia": False, "ver": False, "hor": False} for j in range(len(t)+1)] for i in range(len(s)+1)]
    
    for j in range(1, len(t)+1):
        score_matrix[0][j] = indel_penalty + score_matrix[0][j-1]
        path_matrix[0][j]["hor"] = True

    for i in range(1, len(s)+1):
        max_score = max(indel_penalty + score_matrix[i-1][0], 0)
        score_matrix[i][0] = max_score
        if max_score == score_matrix[i-1][0]+indel_penalty:
            path_matrix[i][0]["ver"] = True

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            score = mismatch_penalty
            if s[i-1] == t[j-1]:
                score = match_reward
            dia_score = score_matrix[i-1][j-1] + score
            ver_score = score_matrix[i-1][j] + indel_penalty
            hor_score = score_matrix[i][j-1] + indel_penalty

            score_matrix[i][j] = max(dia_score, ver_score, hor_score)
            if score_matrix[i][j] == dia_score:
                path_matrix[i][j]["dia"] = True
            if score_matrix[i][j] == hor_score:
                path_matrix[i][j]["hor"] = True
            if score_matrix[i][j] == ver_score:
                path_matrix[i][j]["ver"] = True 
    
    max_index = 0
    curr_max = score_matrix[len(s)][0]
    for j in range(len(t)):
        if score_matrix[len(s)][j] > curr_max:
            curr_max = score_matrix[len(s)][j]
            max_index = j

    f_seq1, f_seq2 = Backtrack(s, t, path_matrix, max_index)
    return curr_max, f_seq1, f_seq2          
    
def Backtrack(seq1: str, seq2: str, path_matrix: list[list[dict[str, int]]], max_idx: int) -> tuple[str, str]:
    """
    Backtrack takes two sequences and a path_matrix then returns an alignment of the two sequences based on the 
    path_matrix.    
    """
    i = len(seq1)
    j = max_idx
    final_seq1 = ""
    final_seq2 = ""

    while j != 0:
        if path_matrix[i][j]["dia"] == True:
            i -= 1
            j -= 1
            final_seq1 = "".join([seq1[i], final_seq1])
            final_seq2 = "".join([seq2[j], final_seq2])
        if path_matrix[i][j]["hor"] == True:
            j -= 1 
            final_seq1 = "".join(["-", final_seq1])
            final_seq2 = "".join([seq2[j], final_seq2])
        if path_matrix[i][j]["ver"] == True:
            i -= 1
            final_seq1 = "".join([seq1[i], final_seq1])
            final_seq2 = "".join(["-", final_seq2])
    return final_seq1, final_seq2
    
if __name__ == "__main__":
    main()