import os, math
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "files")
    paths = glob(os.path.join(dirpath , "Viterbi\\Week*.txt"))
    for path in paths:
        output, emission_vals, state_names, transition_mm, emission_mm = ReadFiles_Viterbi(path)
        state_path = Viterbi(output, emission_vals, state_names, transition_mm, emission_mm)

    anspath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(anspath, "w") as f:
        f.write(str(state_path))

def ReadFiles_Viterbi(path):
    items = ["", [], [], {}, {}]
    keys = []
    idx = 0
    with open(path) as f:
        for line in f:
            if line[0:3] == "---":
                idx += 1
                keys = []
            elif idx == 0:
                items[idx] = line.strip()
            elif idx == 1:
                items[idx] = line.strip().split()
            elif idx == 2:
                items[idx] = line.strip().split()
            elif idx == 3:
                if len(keys) == 0:
                    keys = line.strip().split()
                else:
                    line_0 = line.strip().split()
                    items[idx][line_0[0]] = {}
                    for idx0, key in enumerate(keys):
                        items[idx][line_0[0]][key] = float(line_0[idx0 + 1])
            elif idx == 4:
                if len(keys) == 0:
                    keys = line.strip().split()
                else:
                    line_0 = line.strip().split()
                    items[idx][line_0[0]] = {}
                    for idx0, key in enumerate(keys):
                        if len(line_0[idx0+1].split("/")) > 1:
                            num, denom = line_0[idx0+1].split("/")
                            items[idx][line_0[0]][key] = float(num)/float(denom)
                        else:
                            items[idx][line_0[0]][key] = float(line_0[idx0 + 1])
                        
    
    return tuple(items)

def Viterbi(output:str, emission_vals: list[str], state_vals: list[str], transition_mm: dict[str, dict[str, float]],
            emission_mm: dict[str, dict[str, float]]) -> str:
    """
    Viterbi: Finds the most likely sequence of states that create a specific outcome based on a state transition
    matrix and an emission matrix.

    Input:
        output: The given output.
        emission_vals: The possible emission values.
        state_vals: The possible states.
        transition_mm: The markov model for state transitions.
        emission_mm: The model for emissions based on state.
    
    Output:
        The most likely sequence of states.
    """

    score_matrix = [[0 for y in output] for x in state_vals]
    for x in range(len(state_vals)):
        score_matrix[x][0] = emission_mm[state_vals[x]][output[0]]*transition_mm["F"][state_vals[x]]

    prior_dict = {}
    for x in state_vals:
        prior_dict[x] = False
    path_matrix = [[prior_dict.copy() for y in output] for x in state_vals]

    for x in range(1, len(output)):
        for y in range(len(state_vals)):
            current_state = state_vals[y]
            scores = {}
            for z in range(len(state_vals)):
                prior_state = state_vals[z]
                scores[prior_state] = score_matrix[z][x-1]*\
                    transition_mm[prior_state][current_state]*\
                    emission_mm[current_state][output[x]]
            score_matrix[y][x] = max(scores.values())
            
            for key in scores:
                if score_matrix[y][x] == scores[key]:
                    path_matrix[y][x][key] = True
    
    final_state = 0
    final_state_score = score_matrix[final_state][-1]
    for x in range(1, len(state_vals)):
        if score_matrix[x][-1] > final_state_score:
            final_state_score = score_matrix[x][-1]
            final_state = x

    state_path = Backtrack(final_state, state_vals, path_matrix)
    return state_path

def Backtrack(final_state: int, state_vals: list[str], path_matrix: list[list[dict[str, bool]]]):
    """
    Backtrack: Backtrack through a viterbi graph.
    Input:
        final_state: The state to start backtracking from.
        state_vals: a list of the possible states.
        path_matrix: A matrix that holds the path used to solve a viterbi graph.
    Output:
        The path taken as a string.
    """
    state_path = state_vals[final_state]
    current_state = final_state
    i = len(path_matrix[0])-1
    while i > 0:
        for idx, key in enumerate(state_vals):
            if path_matrix[current_state][i][key] == True:
                state_path = "".join([key, state_path])
                current_state = idx
                break
        i -= 1
    return state_path


if __name__ == "__main__":
    main()