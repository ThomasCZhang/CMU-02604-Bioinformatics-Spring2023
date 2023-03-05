import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "HMMEmission\\data*.txt"))
    for path in paths:
        output, alphabet, hiddenPath, state_names, hmm = ReadFiles_HMMEmission(path)
        prob = EmissionProb(output, alphabet, hiddenPath, state_names, hmm)
        print(prob)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        f.write(str(prob))

def ReadFiles_HMMEmission(path):
    items = ["", [], "", [], {}]
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
                items[idx] = line.strip()
            elif idx == 3:
                items[idx] = line.strip().split()
            elif idx == 4:
                if len(keys) == 0:
                    keys = line.strip().split()
                else:
                    line_0 = line.strip().split()
                    items[idx][line_0[0]] = {}
                    for idx0, key in enumerate(keys):
                        items[idx][line_0[0]][key] = float(line_0[idx0 + 1])

    return items[0], items[1], items[2], items[3], items[4]

def EmissionProb(output: str, alphabet: list[str], hidden_path: str, states: str, emission: dict[str, dict[str, float]]) -> float:
    """
    EmmissionProb calculates the the probability of a given emission sequence based on a emission matrix.
    Input:
        output: The sequence of emission outputs.
        hidden_path: A given hidden path.
        emission: A emission matrix
    """
    prob = 1
    for idx, curr_char in enumerate(output):
        curr_state = hidden_path[idx]
        prob *= emission[curr_state][curr_char]
        curr_char = curr_char
    return prob

if __name__ == "__main__":
    main()