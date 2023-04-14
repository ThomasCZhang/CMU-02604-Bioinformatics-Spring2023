import os
from PartialSuffixArray import *
from BurrowsWheeler import *
from BWTDecode import *

def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "ApproximateMatching", "dataset_876296_6.txt")
    word, patterns, num_mismatches = ReadData_ApproximateMatching(filepath)
    answer = MultiApproximateMatching(word, patterns, num_mismatches)

    answer_path = os.path.join(dirname, "answer.txt")
    with open(answer_path, 'w') as f:
        for key in answer:
            f.write(f'{key}:')
            for val in answer[key]:
                f.write(f' {val}')
            f.write('\n')

def ReadData_ApproximateMatching(filepath: str) -> tuple[str, list[str], int]:
    with open(filepath) as f:
        word = f.readline().strip()
        patterns = f.readline().strip().split()
        num_mismatches = int(f.readline().strip())
    return word, patterns, num_mismatches


def MultiApproximateMatching(word: str, patterns: list[str], num_mismatches: int):
    """
    Finds approximate pattern matches in a word for each pattern in a list of patterns. An approximate match is a match
    that has a limit on the number of mismatched characters.
    Input:
        word: The word being searched for patterns.
        patterns: A list of patterns.
        num_mismatches: The number of allowed mismatches.
    Output:
        A list of starting indexes for each pattern. The starting index indicates the position in the word that the
        pattern starts at.
    """
    word = word + "$"
    suffix_arr = CreateSuffixArray(word)
    bw_word = bwt(word)
    ranks = GetRanks(bw_word)
    count_matrix = GetCharPositions(bw_word)
    match_dict = {}
    for pattern in patterns:
        match_dict[pattern] = ApproximateMatching(pattern, bw_word, num_mismatches, suffix_arr, ranks, count_matrix)

    return match_dict


def ApproximateMatching(
    pattern: str,
    bw_word: str,
    num_mismatches: int,
    suffix_arr: list[int],
    ranks: dict[str, int],
    count_matrix: dict[str, list[int]],
):
    """
    Performs approximate pattern matching for one pattern against a given word.
    Input:
        bw_word: The Burrows-Wheeler transform of the word being searched for a pattern.
        pattern: The pattern being searched for.
        num_mismatches: The number of mismatches allowed.
        suffix_arr: The suffix array corresponding to the bw transformed word
        ranks: the number of characters that are smaller in bw_word.
        count_matrix: The position for each character in bw_word as well as the number of times that character has 
        appeared. 
    Output:
        The locations of the approximate matches in word.
    """
    alphabet = set(ranks.keys())
    alphabet.remove("$")
    top_bot_dict = {"": [0, len(bw_word), 0]}  # Top, bot, number of mismatches so far.
    locations = [] # Will hold all the locations where pattern matches with <= max mismatches allowed.
    for i, letter in enumerate(reversed(pattern)):
        keys = set(top_bot_dict.keys())
        for key in keys:
            top, bot, mm_count = top_bot_dict[key]
            del top_bot_dict[key]
            
            for a in alphabet:
                # Don't bother trying this letter if using it exceeds permissable mismatches.
                if a != letter and mm_count+1 > num_mismatches: 
                    continue
    
                if i != len(pattern) - 1:
                    new_top, new_bot = FindNewTopBot(bw_word, a, top, bot, ranks, count_matrix)
                    if new_top == None:
                        continue
                    elif a != letter:
                        top_bot_dict[key+a] = [new_top, new_bot, mm_count + 1]
                    else:
                        top_bot_dict[key+a] = [new_top, new_bot, mm_count]
                else:
                    locations.extend(FindNewTopBotFinalLetter(bw_word, a, top, bot, ranks, count_matrix))

    for idx, val in enumerate(locations):
        locations[idx] = suffix_arr[val]
    return locations

def FindNewTopBot(
    bw_word: str, letter: str, top_idx: int, bot_idx: int, ranks: dict[str, int], count_matrix: dict[str, list[int]]
):
    """
    Finds the top and bottom index for the next iteration of BW suffix pattern matching.
    Input:
        bw_word: The Burrows-Wheeler transform of the word being searched for a pattern.
        letter: The letter being used to determine the next top and bot index.
        top_idx: The current top_idx
        bot_idx: the current bot_idx
        ranks: the number of characters that are smaller in bw_word.
        count_matrix: The position for each character in bw_word as well as the number of times that character has 
        appeared. 
    Return:
        A tuple of two integers. (new top idx, new bot idx)
    """
    new_top = None
    new_bot = None
    for i in range(top_idx, bot_idx):
        if bw_word[i] == letter:
            if new_top == None:
                new_top = ranks[letter] + count_matrix[letter][i]
                new_bot = new_top + 1
            else:
                new_bot += 1
    return new_top, new_bot

def FindNewTopBotFinalLetter(
    bw_word: str, letter: str, top_idx: int, bot_idx: int, ranks: dict[str, int], count_matrix: dict[str, list[int]]
):
    """
    Finds the top and bottom index for the next iteration of BW suffix pattern matching.
    Input:
        bw_word: The Burrows-Wheeler transform of the word being searched for a pattern.
        letter: The letter being used to determine the next top and bot index.
        top_idx: The current top_idx
        bot_idx: the current bot_idx
        ranks: the number of characters that are smaller in bw_word.
        count_matrix: The position for each character in bw_word as well as the number of times that character has 
        appeared. 
    Return:
        A tuple of two integers. (new top idx, new bot idx)
    """
    final_rows = []
    for i in range(top_idx, bot_idx):
        if bw_word[i] == letter:
            final_rows.append(ranks[letter]+count_matrix[letter][i])
    return final_rows


def CreateSuffixArray(word: str) -> int:
    """
    Creates the suffix array of a word.
    Input:
        word: The word of a suffix array.
    """
    suffixes = [[] for _ in range(len(word))]
    for i in range(len(word)):
        suffixes[i] = [i, word[i:] + word[:i]]
    suffixes = sorted(suffixes, key=lambda x: x[1])
    suffixes = [val[0] for val in suffixes]
    return suffixes


if __name__ == "__main__":
    main()
