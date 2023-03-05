import os
import numpy as np
from GreedySorting import ReadPermutationData

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "inputs\Breakpoint")
    filepath = os.path.join(dirpath, "input_0.txt")
    permutation = ReadPermutationData(filepath)
    answer = CountBreakPoints(permutation)
    
    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        f.write(str(answer))

def CountBreakPoints(p: np.ndarray[int]) -> int:
    """
    CountBreakPoints: Counts the number of breakpoints that occur in a list of integers.
    A pair of consecutive numbbers in p that don't follow counting sequence or inverted counting sequence contribute
    1 to a break point. 
    For instance: 3, 4 and -4, -3 are not break points, however 4, 3 and -3, 4 are. 1,3 also adds to break points
    because we skip over the number 2.
    
    Input:
        p: The permutation data stored as a numpy array.
    """
    num_breakpoints = 0 # Counter for number of break points.
    temp_p = np.concatenate([[0], p, [p.size+1]])
    for idx in range(temp_p.size-1):
        if temp_p[idx+1] != temp_p[idx]+1:
            num_breakpoints += 1
    return num_breakpoints


if __name__ == main():
    main()