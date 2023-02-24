import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "HMMProb\\data*.txt"))
    for path in paths:
        hiddenPath, state_names, hmm = ReadFiles_HMMProb(path)
        prob = HMMProb(hiddenPath, hmm)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        f.write(str(prob))

def ReadFiles_HMMProb(path):
    items = ["", [], {}]
    keys = []
    idx = 0
    with open(path) as f:
        for line in f:
            if line[0:3] == "---":
                idx += 1
            elif idx == 0:
                items[idx] = line.strip()
            elif idx == 1:
                items[idx] = line.strip().split()
            elif idx == 2:
                if len(keys) == 0:
                    keys = line.strip().split()
                else:
                    line_0 = line.strip().split()
                    items[idx][line_0[0]] = {}
                    for idx0, key in enumerate(keys):
                        items[idx][line_0[0]][key] = float(line_0[idx0 + 1])

    return items[0], items[1], items[2] 

def HMMProb(states: str, transition: dict[str, dict[str, float]]) -> float:
    """
    HMMProb: Calculates the probability of a hidden path based on a transition matrix.
    Input:
        states: the provided hidden path.
        transition: The transition probabilities between states.
    Output:
        The probability that the given transition state occurs based on the provided transition matrix.
    """
    prob = 1/(len(transition))
    curr_char = states[0]
    for next_char in states[1:]:
        prob *= transition[curr_char][next_char]
        curr_char = next_char
    return prob

if __name__ == "__main__":
    main()