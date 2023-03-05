import os
from glob import glob
import numpy as np

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "inputs\GreedySort")
    filepaths = glob(dirpath+"\input_0.txt")
    for path in filepaths:
        permutation = ReadPermutationData(path)
        answer = GreedySorting(permutation)
    
    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        for idx0, permutation in enumerate(answer):
            if idx0 != 0:
                f.write("\n")
            for idx1, x in enumerate(permutation):
                if idx1 != 0:
                    f.write(" ")
                if x > 0:
                    f.write("+")
                f.write(str(x))
    
def ReadPermutationData(path: str) -> np.ndarray:
    with open(path) as f:
        permutation_list = [int(x) for x in f.readline().strip().split()]
    return np.asarray(permutation_list)

def GreedySorting(p: np.ndarray) -> list[list[int]]:
    """
    GreedySorting: Sorts a list of numbers via "burnt pancake flipping"
    Algorithm is greedy.
    Input:
        p: The list of integers to be sorted. Assumed that the absolute value of all integers in p are unique.
    Output:
        A list of the the permutations required to sort p.
    """
    permutation_list = []
    for idx in range(p.shape[0]):
        if np.abs(p[idx] != idx+1):
            p = PlaceValueInPosition(idx, p)
            permutation_list.append(p)
        if p[idx] < 0: # if current value is negative we need to flip again
            p = p.copy()
            p[idx] *= -1
            permutation_list.append(p)
    return permutation_list

def PlaceValueInPosition(val: int, p: np.ndarray)->np.ndarray:
    """
    PlaceValueInPosition: Moves the the value of val+1 into the val index of p via flipping.
    Note: Sign of specified value irrelevant. Assumed that |p| contains unique integers starting from 1 and ending
    at the length of p. 
    Input:
        val: The absolute value of the value we are trying to move to the correct position via flipping.
        p: The array of values.
    Output:
        The new array.
    """

    idx_of_val = np.nonzero(np.abs(p)==val+1)[0][0]
    if val-1 < 0:
        p = np.concatenate([-1*p[idx_of_val::-1], p[idx_of_val+1:]])
    else:
        p = np.concatenate([p[:val],-1*p[idx_of_val:val-1:-1], p[idx_of_val+1:]])
    return p

if __name__ == "__main__":
    main()