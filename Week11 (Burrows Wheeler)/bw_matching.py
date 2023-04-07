import os
from bwt_decode import *

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "bwt_matching", "dataset_876292_8.txt")
    with open(filepath) as f:
        word = f.readline().strip()
        patterns = f.readline().strip().split()
    answer = bwt_multi_matching(word, patterns)

    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, "w") as f:
        for i, match_indexes in enumerate(answer):
            if i != 0:
                f.write(" ")
            f.write(str(len(match_indexes)))

def bwt_multi_matching(word: str, patterns: list[str]) -> list[list[int]]:
    """
    Finds the location of a list of patterns in a bwt encoded string.
    Input:
        word: the bwt encoded string being matched against.
        patterns: the patterns we are looking for in word.
    Output:
        A list of indicies indicating the starting location of where each pattern can be found in word.
    """
    ranks = get_ranks(word)
    idx_dict = get_char_positions(word)
    matches = [[] for _ in patterns]
    for i, pattern in enumerate(patterns):
        matches[i] = bwt_matching(word, pattern, ranks, idx_dict)
    
    return matches

def bwt_matching(word: str, pattern: str, ranks: list[int] = None, idx_dict: dict[int, dict[int, int]] = None) -> list[int]:
    """
    Finds all the locations of a pattern in a bwt encoded string.
    Input:
        word: the bwt encoded string
        pattern: The pattern being searched.
    Output:
        the indicies of the starting location of pattern in word.
    """
    if ranks is None:
        ranks = get_ranks(word)
    if idx_dict is None:
        idx_dict = get_char_positions(word)
    
    locations = []
    top_idx = ranks[pattern[-1]]
    bottom_idx = top_idx + len(idx_dict[pattern[-1]])
    for i, letter in enumerate(reversed(pattern[:-1])):
        found_first = False
        for j in range(top_idx, bottom_idx):
            if word[j] == letter:
                if i == len(pattern)-2:
                    locations.append(j)
                elif not found_first:
                    top_idx = ranks[letter] + idx_dict[letter][j]
                    bottom_idx = top_idx + 1
                    found_first = True
                else:
                    bottom_idx += 1
        if not found_first:
            break
                
    return locations


if __name__ == "__main__":
    main()