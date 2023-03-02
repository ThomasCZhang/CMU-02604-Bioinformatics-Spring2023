import os
import numpy as np
from glob import glob
from HMMPseudoCount import ProfileHMMPseudo

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "SeqAlign\\input_1.txt"))
    for path in paths:
        x, threshold, pseudocount, alphabet, alignment= ReadProfileHMM(path)
        answer = SeqAlignHMM(x, threshold, pseudocount, alphabet, alignment)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        for idx, state in enumerate(answer):
            if idx != 0:
                f.write(" ")
            f.write(state)

def ReadProfileHMM(path):
    items = [[] for _ in range(5)]
    item_num = 0
    with open(path) as f:
        for line in f:
            if line.strip() == "--------":
                item_num += 1
            elif item_num == 0:
                items[0] = line.strip()
            elif item_num == 1:
                line_list = line.strip().split()
                items[item_num] = float(line_list[0])
                item_num += 1
                items[item_num] = float(line_list[1])
            elif item_num == 3:
                items[item_num] = line.strip().split()
            elif item_num >= 4:
                items[item_num].append(line.strip())

    return tuple(items)

def SeqAlignHMM(x: str, threshold: float, pseudocount: float, alphabet:list[str], alignment: list[str]) -> list[str]:
    """
    SeqAlignHMM: Returns the most likely sequence of states to match x against an pre-existing sequence alignment.

    Input:
        x: The sequence being matched up against an alignment.
        threshold: The threshold gap frequency before a column in alignment becomes an insert column.
        pseudocount: The pseudocount values to add to the transition and emission matrix.
        alignment: The existing sequence alignment.
    
    Output:
        A list of strings denoting the sequence of the states.
    """
    transition, emission = ProfileHMMPseudo(threshold, pseudocount, alphabet, alignment)
    # Taking the log of both transition and emission probabilities so we can add instead of multiply.
    transition = np.log(transition) 
    emission = np.log(emission)

    # Alphabet dictionary to keep track of column index for each letter in emissions.
    alphabet_dict = {} 
    for idx, key in enumerate(alphabet):
        alphabet_dict[key] = idx

    num_states = transition.shape[0]
    score_matrix = np.zeros((num_states-2, len(x)+1))
    path_matrix = np.zeros((num_states-2, len(x)+1, 2))

    # score_matrix[0,0] = 1 # This is needed if we are multiplying. 

    # Fill in first column 
    for row in range(1, score_matrix.shape[0]):
        if row%3 == 2:
            prev_row = row-3
            if prev_row < 0:
                prev_row += 1
            score_matrix[row, 0] = transition[row-2, row+1]+score_matrix[prev_row, 0]
            path_matrix[row, 0, :] = [prev_row, 0]
        else:
            score_matrix[row, 0] = -np.inf
    
    # Fill in first three rows because they behave differently
    for col in range(1, score_matrix.shape[1]):
        prev_state = 1
        if col == 1:
            prev_state = 0
        
        score_matrix[0, col] = emission[1, alphabet_dict[x[col-1]]] + transition[prev_state, 1] + score_matrix[0, col-1]
        path_matrix[0, col, :] = [0, col-1]

        score_matrix[1, col] = emission[2, alphabet_dict[x[col-1]]] + transition[prev_state, 2] + score_matrix[0, col-1]
        path_matrix[1, col, :] = [0, col-1]

        score_matrix[2, col] = transition[1, 3] + score_matrix[0, col]
        path_matrix[2, col, :] = [0, col]


    for col in range(1, score_matrix.shape[1]): # Fill in column by column
        for row in range(3, score_matrix.shape[0]):
            if row%3 == 0: # Insert node
                if col > 3 and row > 14:
                    print("Hi")
                m_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row-1, row+1] + score_matrix[row-2, col-1]
                d_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row, row+1] + score_matrix[row-1, col-1]
                i_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row+1, row+1] + score_matrix[row, col-1]
                score_matrix[row, col] = max(m_score, d_score, i_score)
                if score_matrix[row, col] == m_score:
                    path_matrix[row, col, :] = [row-2, col-1]
                elif score_matrix[row, col] == d_score:
                    path_matrix[row, col, :] = [row-1, col-1]
                elif score_matrix[row, col] == i_score:
                    path_matrix[row, col, :] = [row, col-1]

            elif row%3 == 1: # Match Node
                m_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row-2, row+1] + score_matrix[row-3, col-1]
                d_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row-1, row+1] + score_matrix[row-2, col-1]
                i_score = emission[row+1, alphabet_dict[x[col-1]]] + transition[row, row+1] + score_matrix[row-1, col-1]
                score_matrix[row, col] = max(m_score, d_score, i_score)
                if score_matrix[row, col] == m_score:
                    path_matrix[row, col, :] = [row-3, col-1]
                elif score_matrix[row, col] == d_score:
                    path_matrix[row, col, :] = [row-2, col-1]
                elif score_matrix[row, col] == i_score:
                    path_matrix[row, col, :] = [row-1, col-1]
            elif row%3 == 2: # Delete node
                m_score = transition[row-3, row] + score_matrix[row-4, col]
                d_score = transition[row-3, row] + score_matrix[row-3, col]
                i_score = transition[row-3, row] + score_matrix[row-2, col]
                score_matrix[row, col] = max(m_score, d_score, i_score)
                if score_matrix[row, col] == m_score:
                    path_matrix[row, col, :] = [row-4, col]
                elif score_matrix[row, col] == d_score:
                    path_matrix[row, col, :] = [row-3, col]
                elif score_matrix[row, col] == i_score:
                    path_matrix[row, col, :] = [row-2, col]

    max_score = -np.inf
    start_pos = [0, 0]
    for idx in range(1,4):
        num_rows, num_cols = score_matrix.shape
        if score_matrix[num_rows-idx, -1] > max_score:
            max_score = score_matrix[num_rows-idx, num_cols-1]
            start_pos = [num_rows-idx, num_cols-1]

    backtrack_path = backtrack(path_matrix, start_pos)
    return backtrack_path

def backtrack(path_matrix: np.ndarray[int], start_position = list[int, int]) -> list[str]:
    """
    Backtracks through the viterbi HMM graph from a given start position using a path matrix.
    """
    row, col = start_position
    path = []
    while row>0 or col>0:
        state_num = int((row+2)//3)
        if row%3 == 0:
            state_letter = "I"
        if row%3 == 1:
            state_letter = "M"
        if row%3 == 2:
            state_letter = "D"
        state = state_letter + str(state_num)
        path = [state] + path
        row, col = path_matrix[int(row), int(col)]
    
    return path

    

if __name__ == "__main__":
    main()