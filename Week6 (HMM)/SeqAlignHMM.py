import os
import numpy as np
from typing import Iterable
from glob import glob
from HMMPseudoCount import ProfileHMMPseudo

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "SeqAlign\\data*.txt"))
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
    transition = np.log(transition, out = -np.inf*np.ones_like(transition), where= transition!= 0) 
    emission = np.log(emission, out = -np.inf*np.ones_like(emission), where = emission != 0)

    # Alphabet dictionary to keep track of column index for each letter in emissions.
    alphabet_dict = {} 
    for idx, key in enumerate(alphabet):
        alphabet_dict[key] = idx

    num_states = transition.shape[0]
    score_matrix = np.ones((num_states-2, len(x)+1))*(-np.inf)
    path_matrix = np.zeros((num_states-2, len(x)+1, 2))

    score_matrix[0,0] = 0 

    # Dictionary to hold functions so they can be called efficiently. (Instead of having a bunch of if statements)
    switch_dict = {0: InsertNode, 1: MatchNode, 2: DeleteNode} 

    for col in range(score_matrix.shape[1]): # Fill in column by column
        for row in range(score_matrix.shape[0]):
            if row == 0 and col == 0: # Must handle first position differently than the rest.
                # Checking Child Delete
                if TwoDimInBounds(row+2, col, score_matrix.shape) and \
                score_matrix[row, col]+transition[row, row+3] > score_matrix[row+2, col]:
                    score_matrix[row+2, col] = score_matrix[row, col]+transition[row, row+3]
                    path_matrix[row+2,col] = [row,col]
                # Checking Child Insert    
                if TwoDimInBounds(row, col+1, score_matrix.shape) and \
                score_matrix[row, col]+transition[row, row+1]+emission[row+1, alphabet_dict[x[col]]] > score_matrix[row, col+1]:
                    score_matrix[row, col+1] = score_matrix[row, col]+transition[row, row+1]+emission[row+1, alphabet_dict[x[col]]]                               
                    path_matrix[row,col+1] = [row,col]
                # Checking child match
                if TwoDimInBounds(row+1, col+1, score_matrix.shape) and\
                score_matrix[row, col] + transition[row, row+2] + emission[row+2, alphabet_dict[x[col]]] > score_matrix[row+1, col+1]:
                    score_matrix[row+1, col+1] = score_matrix[row, col] + transition[row, row+2] + emission[row+2, alphabet_dict[x[col]]]
                    path_matrix[row+1,col+1] = [row,col]
            else:
                node_func = switch_dict[int(row%3)]
                node_func(row, col, x, score_matrix, path_matrix, transition, emission, alphabet_dict)

    # Find the ending node of the max score path.
    max_score = -np.inf
    start_pos = [0, 0]
    for idx in range(1,4):
        num_rows, num_cols = score_matrix.shape
        if score_matrix[num_rows-idx, -1] > max_score:
            max_score = score_matrix[num_rows-idx, num_cols-1]
            start_pos = [num_rows-idx, num_cols-1]

    backtrack_path = backtrack(path_matrix, start_pos)
    return backtrack_path

def MatchNode(row: int, col:int, x: str, score_matrix: np.ndarray[float], path_matrix: np.ndarray[int],
               transition: np.ndarray[float], emission: np.ndarray[float], alphabet_dict: dict[str, int]):
    """
    MatchNode: Updates the children of a match node.
    Input:
        row, col: row and col indicies of node being evaluated.
        x: the string being evaluated.
        score_matrix: a score_matrix used for DP.
        path_matrix: a path matrix used for backtracking.
        transition: Transition matrix.
        emission: emission matrix.
    Output:
        None
    """
    # Checking Child Delete
    if TwoDimInBounds(row+4, col, score_matrix.shape) and\
    score_matrix[row,col]+transition[row+1, row+5]>score_matrix[row+4,col]:
        score_matrix[row+4,col] = score_matrix[row,col]+transition[row+1, row+5]
        path_matrix[row+4,col] = [row,col]
    # Checking Child Insert
    if TwoDimInBounds(row+2, col+1, score_matrix.shape) and\
    score_matrix[row,col]+transition[row+1,row+3]+emission[row+3, alphabet_dict[x[col]]] > score_matrix[row+2, col+1]:
        score_matrix[row+2, col+1] = score_matrix[row,col]+transition[row+1,row+3]+emission[row+3, alphabet_dict[x[col]]]
        path_matrix[row+2,col+1] = [row,col]
    # Checking Child Match
    if TwoDimInBounds(row+3, col+1, score_matrix.shape) and\
    score_matrix[row,col]+transition[row+1, row+4]+emission[row+4, alphabet_dict[x[col]]] > score_matrix[row+3, col+1]:
        score_matrix[row+3, col+1] = score_matrix[row,col]+transition[row+1, row+4]+emission[row+4, alphabet_dict[x[col]]]
        path_matrix[row+3,col+1] = [row,col]
    
def InsertNode(row: int, col:int, x: str, score_matrix: np.ndarray[float], path_matrix: np.ndarray[int],
               transition: np.ndarray[float], emission: np.ndarray[float], alphabet_dict: dict[str, int]):
    """
    InsertNode: Updates the children of a match node.
    Input:
        row, col: row and col indicies of node being evaluated.
        x: the string being evaluated.
        score_matrix: a score_matrix used for DP.
        path_matrix: a path matrix used for backtracking.
        transition: Transition matrix.
        emission: emission matrix.
    Output:
        None
    """
    # Checking Child Delete
    if TwoDimInBounds(row+2, col, score_matrix.shape) and \
    score_matrix[row, col]+transition[row+1, row+3] > score_matrix[row+2, col]:
        score_matrix[row+2, col] = score_matrix[row, col]+transition[row+1, row+3]
        path_matrix[row+2,col] = [row,col]
    # Checking Child Insert    
    if TwoDimInBounds(row, col+1, score_matrix.shape) and \
    score_matrix[row, col]+transition[row+1, row+1]+emission[row+1, alphabet_dict[x[col]]] > score_matrix[row, col+1]:
        score_matrix[row, col+1] = score_matrix[row, col]+transition[row+1, row+1]+emission[row+1, alphabet_dict[x[col]]]                               
        path_matrix[row,col+1] = [row,col]
    # Checking child match
    if TwoDimInBounds(row+1, col+2, score_matrix.shape) and\
    score_matrix[row, col] + transition[row+1, row+2] + emission[row+2, alphabet_dict[x[col]]] > score_matrix[row+1, col+1]:
        score_matrix[row+1, col+1] = score_matrix[row, col] + transition[row+1, row+2] + emission[row+2, alphabet_dict[x[col]]]
        path_matrix[row+1,col+1] = [row,col]


def DeleteNode(row: int, col:int, x: str, score_matrix: np.ndarray[float], path_matrix: np.ndarray[int],
               transition: np.ndarray[float], emission: np.ndarray[float], alphabet_dict: dict[str, int]):
    """
    DeleteNode: Updates the children of a match node.
    Input:
        row, col: row and col indicies of node being evaluated.
        x: the string being evaluated.
        score_matrix: a score_matrix used for DP.
        path_matrix: a path matrix used for backtracking.
        transition: Transition matrix.
        emission: emission matrix.
    Output:
        None
    """
    # Checking Child Delete
    if TwoDimInBounds(row+3, col, score_matrix.shape) and \
    score_matrix[row,col] + transition[row+1, row+4]>score_matrix[row+3,col]:
        score_matrix[row+3, col] = score_matrix[row,col] + transition[row+1, row+4]
        path_matrix[row+3,col] = [row,col]
    # Checking Child Insert
    if TwoDimInBounds(row+1, col+1, score_matrix.shape) and\
    score_matrix[row,col] + transition[row+1, row+2] + emission[row+2, alphabet_dict[x[col]]] > score_matrix[row+1, col+1]:
        score_matrix[row+1, col+1] = score_matrix[row,col] + transition[row+1, row+2] + emission[row+2, alphabet_dict[x[col]]]
        path_matrix[row+1, col+1] = [row,col]
    # Checking Child Match
    if TwoDimInBounds(row+2, col+1, score_matrix.shape) and\
    score_matrix[row,col] +transition[row+1, row+3] + emission[row+3, alphabet_dict[x[col]]] > score_matrix[row+2, col+1]:
        score_matrix[row+2, col+1] = score_matrix[row,col] +transition[row+1, row+3] + emission[row+3, alphabet_dict[x[col]]]
        path_matrix[row+2,col+1] = [row,col]
    pass

def TwoDimInBounds(row: int, col: int, shape: Iterable[int]) -> bool:
    """
    TwoDimInBounds: Determines if a given row and column are in bounds a two dimensional matrix.
    Input:
        row: the input row
        col: the input col
        shape: The shape of the matrix.
    Output:
        True if row and col are in bounds, false otherwise.
    """
    if row < 0 or row >= shape[0] or col < 0 or col >= shape[1] :
        return False
    return True

def backtrack(path_matrix: np.ndarray[int], start_position: list[int]) -> list[str]:
    """
    Backtracks through the viterbi HMM graph from a given start position using a path matrix.
    Input:
        path_matrix: The path_matrix that will be used for backtracking.
        start_position: The starting position from which to backtrack.
    Output:
        A list of strings corresponding to the path.
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