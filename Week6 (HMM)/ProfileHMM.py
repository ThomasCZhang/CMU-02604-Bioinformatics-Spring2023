import os
import numpy as np
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "ProfileHMM\\data*.txt"))
    for path in paths:
        threshold, alphabet, alignment= ReadProfileHMM(path)
        transition, emission = ProfileHMM(threshold, alphabet, alignment)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    num_states = transition.shape[0]
    
    with open(anspath, "w") as f:
        # Write row headers
        for i in range(num_states):
            character = DetermineStateName(i, num_states)
            f.write('\t'+character)

        # Write rows
        for i in range(num_states):
            f.write('\n')
            character = DetermineStateName(i, num_states)
            f.write(character)         
            for j, val in enumerate(transition[i]):
                if val != 0:
                    val = str(val)
                else:
                    val = "0"
                f.write("\t"+val)

        f.write("\n"+ "-"*8+"\n")

        # Write Row headers
        for character in alphabet:
            f.write("\t"+character)

        # Write Rows
        for i in range(num_states):
            f.write('\n')
            character = DetermineStateName(i, num_states)
            f.write(character)
            for j, val in enumerate(emission[i]):
                if val != 0:
                    val = str(val)
                else:
                    val = "0"
                f.write("\t"+val)
         

def DetermineStateName(i, num_states):
    """
    Determines the state name for writing the hmm to file.
    Input:
        i: number of the current state
        num_states: total number of states
    """
    letter_dict = {0: "D", 1: "I", 2: "M"}
    if i == 0:
        character = 'S'
    elif i == 1:
        character = 'I0'
    elif i == num_states-1:
        character = 'E'
    else:
        num = (i+1)//3
        character = letter_dict[i%3]+str(num)
    return character

def ReadProfileHMM(path):
    items = [[] for _ in range(3)]
    item_num = 0
    with open(path) as f:
        for line in f:
            if line.strip() == "--------":
                item_num += 1
            elif item_num == 0:
                items[item_num] = float(line.strip())
            elif item_num == 1:
                items[item_num] = line.strip().split()
            elif item_num >= 2:
                items[item_num].append(line.strip())

    return tuple(items)

def ProfileHMM(threshold: float, alphabet: list[str], alignment: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    ProfileHMM: Returns the transition and emission matricies for a given sequence alignment.

    Input:
        threshold: The percentage of blanks a non-insert column can have.
        alphabet: the potential emission symbols
        alignment: the sequence alignment
    Output:
        transition matrix, emission matrix as a tuple of np.ndarrays.
    """

    alignment = np.asarray([[x for x in sequence] for sequence in alignment])
    cols = GetColsBelowThreshold(alignment, threshold)
    transition, emission = InitializeProfileHMM(len(cols), alphabet)

    alphabet_dict = {}
    for i, letter in enumerate(alphabet):
        alphabet_dict[letter] = i
    
    for row in alignment:
        last_match = 0
        curr_idx = 0 # The current index.
        last_idx = 0 
        for c_idx, character in enumerate(row):
            if c_idx in cols: # We are currently in a column that can be a match
                if character != "-": # We are matching
                    curr_idx = (last_match+1)*3-1
                    emission[curr_idx][alphabet_dict[character]] += 1
                else: # We are deleting.
                    curr_idx = (last_match+1)*3
                transition[last_idx,curr_idx] += 1
                last_match += 1
            else: # We are in a column that can be an insert.
                if character != "-": # We are inserting
                    curr_idx = (last_match)*3 + 1
                    transition[last_idx,curr_idx] += 1
                    emission[curr_idx][alphabet_dict[character]] += 1
            
            last_idx = curr_idx    
    
            if c_idx == len(row)-1:
                transition[last_idx,-1] += 1
    
    row_sum = np.sum(transition, axis = 1)[:, np.newaxis]
    transition = np.divide(transition, row_sum, out = np.zeros(transition.shape), where= row_sum!=0)
    row_sum = np.sum(emission, axis = 1)[:, np.newaxis]
    emission = np.divide(emission, row_sum, out = np.zeros(emission.shape), where = row_sum!= 0)

    return transition, emission

def InitializeProfileHMM(n: int, alphabet: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    Initializes the transition and emission matricies for a profile hmm.
    Input:
        n: the number of match states in the profile hmm
        alphabet: a list of the possible characters.
    Output:
        The initialized transition and emission matrices as np.ndarrays.
    """
    transition = np.zeros((3*n+3, 3*n+3))
    emission = np.zeros((3*n+3, len(alphabet)))
    return transition, emission

def GetColsBelowThreshold(alignment: np.ndarray[str], threshold: float) -> set[int]:
    """
    GetColsBelowThreshold: Gets the column index of columns in alignment that don't exceed the gap threshold.

    Input:
        alignment: The alignment
        threshold: The threshold percentage.
    Output:
        the number of columns above the threshold
    """
    n_rows = alignment.shape[0]
    gap_matrix = np.zeros(alignment.shape)
    gap_idx = np.asarray(alignment == "-").nonzero()
    gap_matrix[gap_idx] = 1
    num_gaps_for_each_col = np.sum(gap_matrix, axis = 0)
    percent_gaps_for_each_col = num_gaps_for_each_col/n_rows
    num_cols_below  = np.asarray(percent_gaps_for_each_col < threshold).nonzero()
    return set(num_cols_below[0])

if __name__ == "__main__":
    main()