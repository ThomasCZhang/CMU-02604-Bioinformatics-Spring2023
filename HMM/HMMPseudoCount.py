import os
import numpy as np
from glob import glob
from ProfileHMM import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "HMMPseudo\\trial*.txt"))
    for path in paths:
        threshold, pseudocount, alphabet, alignment= ReadProfileHMM(path)
        hmm, emission = ProfileHMMPseudo(threshold, pseudocount, alphabet, alignment)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    num_states = hmm.shape[0]
    
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
            for j, val in enumerate(hmm[i]):
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
    items = [[] for _ in range(4)]
    item_num = 0
    with open(path) as f:
        for line in f:
            if line.strip() == "--------":
                item_num += 1
            elif item_num == 0:
                line_list = line.strip().split()
                items[item_num] = float(line_list[0])
                item_num += 1
                items[item_num] = float(line_list[1])
            elif item_num == 2:
                items[item_num] = line.strip().split()
            elif item_num >= 3:
                items[item_num].append(line.strip())

    return tuple(items)

def ProfileHMMPseudo(threshold: float,pseudocount: float, alphabet: list[str], alignment: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    ProfileHMMPseudo: Creates a profile hmm with pseudocounts.
    Input:
        threshold: The percentage cutoff for how many gaps a position can have before it stops being considered a 
        "match" column.
        pseudocount: the pseudocount value to add to the transition and emission tables.
        alphabet: the possible emission characters.
        alignment: An alignment of strings consisting of characters from the alphabet. Used to generate transition
        and emission matricies.
    Output:
        The transition and emission tables with pseudocounts as np.ndarrays.
    """
    transition, emission = ProfileHMM(threshold, alphabet, alignment)
    alignment = np.asarray([[x for x in sequence] for sequence in alignment])
    cols = GetColsBelowThreshold(alignment, threshold)

    # Add pseudocount values.
    for idx in range(transition.shape[0]-1):
        match_idx = (idx+1)//3
        start_idx = (match_idx+1)*3-2
        stop_idx = start_idx + 3
        if match_idx == len(cols):
            stop_idx -= 1
        transition[idx][start_idx:stop_idx] += pseudocount
        if idx%3 != 0:
            emission[idx] += pseudocount

    row_sum = np.sum(transition, axis = 1)[:, np.newaxis]
    transition = np.divide(transition, row_sum, out = np.zeros_like(transition), where= row_sum!=0)
    row_sum = np.sum(emission, axis = 1)[:, np.newaxis]
    emission = np.divide(emission, row_sum, out = np.zeros_like(emission), where = row_sum != 0)

    return transition, emission


if __name__ == "__main__":
    main()