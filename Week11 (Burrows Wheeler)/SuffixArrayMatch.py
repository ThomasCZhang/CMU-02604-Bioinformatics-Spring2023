import os
from PartialSuffixArray import *
from BurrowsWheeler import *
from BWTDecode import *


def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "SuffixArrayMatch", "dataset_876295_4.txt")
    word, patterns = ReadData(filepath)
    match_locations = SuffixArrayMatch(word, patterns)

    answerpath = os.path.join(dirname, 'answer.txt')
    with open(answerpath, 'w') as f:
        for key in match_locations:
            f.write(f"{key}:")
            for val in match_locations[key]:
                f.write(f' {val}')
            f.write('\n')

def ReadData(filepath):
    with open(filepath) as f:
        word = f.readline().strip()
        patterns = f.readline().strip().split()
    return word, patterns


def SuffixArrayMatch(word: str, patterns: list[str]) -> dict[str, list[int]]:
    """
    Performs suffix array pattern matching.
    Input:
        word: The word we are trying to match
        patterns: the patterns we are trying to match.
    Output:
        A dictionary that holds the starting locations of each pattern.
    """
    k = 5
    word = word + "$"
    bw_word = bwt(word)
    ranks = GetRanks(bw_word)
    partial_suffix = PartialSuffixArray(word, k)
    partial_count = PartialCountMatrix(word, k)
    match_locations = {}
    for pattern in patterns:
        match_locations[pattern] = PartialSuffixArrayMatch(pattern, bw_word, ranks, partial_suffix, partial_count, k)
    return match_locations


def PartialCountMatrix(word: str, k: int) -> dict[str, list[int]]:
    """
    Creates and returns a partical count dictionary for a bwt word.
    Input:
        bwt_word: A Burrows-Wheeler transformed word.
        k: The increment used to create the partical count matrix.
    Output:
        The partical count matrix as a dictionary.
    """
    bw_word = bwt(word)
    alphabet = set(bw_word)
    counter_dict = {}
    count_matrix_dict = {}
    for letter in alphabet:
        counter_dict[letter] = 0
        count_matrix_dict[letter] = []

    for i, letter in enumerate(bw_word):
        if i % 5 == 0:
            for key in counter_dict:
                count_matrix_dict[key].append(counter_dict[key])
        counter_dict[letter] += 1

    return count_matrix_dict


def PartialSuffixArrayMatch(
    pattern: str,
    bw_word: str,
    ranks: dict[str, int],
    partial_suffix: dict[int, int],
    partial_count: dict[str, list[int]],
    k: int,
):
    """
    Finds the locations in a word that matches a given pattern.
    Input:
        pattern: The pattern being matched.
        bw_word: The Burrows-Wheeler transformed word.
        ranks: A dictionary that holds the "rank" of each character that appears in the word. The rank
            of a character is the number of "smaller" characters that appear in the word.
        partial_suffx: The partial suffix array of the word.
        partial_count: The partial count matrix of the word.
        k: The "increment size" of the partial suffix array and partial count matrix.
    """
    alphabet = sorted(ranks.keys())
    alphabet_ind = alphabet.index(pattern[-1])
    top_idx = ranks[alphabet[alphabet_ind]]
    if alphabet_ind != len(alphabet) - 1:
        bot_idx = ranks[alphabet[alphabet_ind + 1]]
    else:
        bot_idx = len(bw_word)

    locations = []
    # print(f"Word: {bw_word}\nRanks {ranks}\tPattern: {pattern}")
    for i, letter in enumerate(reversed(pattern[:-1])):
        # print(f"Top Index: {top_idx}\tBot Index: {bot_idx}\tLetter: {letter}")
        found_first = False
        for j in range(top_idx, bot_idx):
            if bw_word[j] == letter:
                if i == len(pattern) - 2:
                    locations.append(
                        ranks[letter] + GetCountFromPartialCountMatrix(bw_word, partial_count, letter, j, k)
                    )
                elif not found_first:
                    top_idx = ranks[letter] + GetCountFromPartialCountMatrix(bw_word, partial_count, letter, j, k)
                    bot_idx = top_idx + 1
                    found_first = True
                else:
                    bot_idx += 1
        if not found_first:
            break
    # print(f"locations: {locations}")

    suffix_indexes = []
    for i in locations:
        suffix_indexes.append(GetStartingPosition(bw_word, ranks, partial_suffix, partial_count, i, k))
    # print(suffix)
    return suffix_indexes


def GetStartingPosition(
    bw_word: str,
    ranks: dict[str, int],
    partial_suffix: dict[int, int],
    partial_count: dict[str, list[int]],
    i: int,
    k: int,
) -> int:
    """
    Gets the suffix number of a given position using a partial suffix array and a Burrows-Wheeler transformed word.
    Input:
        bw_word: A Burrows-Wheeler transformed word.
        ranks: A dictionary that holds the "rank" of each character that appears in the word. The rank
            of a character is the number of "smaller" characters that appear in the word.
        partial_suffix: The partial suffix array.
        partial_count: The partial count matrix
        i: The index we are trying to find the suffix value of.
        k: The increment size of the partial suffix array and partial count.
    Output:
        The suffix array index as an integer
    """
    num_steps = 0
    current_idx = i
    while True:
        if current_idx in partial_suffix:  # Stop condition.
            break
        pass
        letter = bw_word[current_idx]
        current_idx = ranks[letter] + GetCountFromPartialCountMatrix(bw_word, partial_count, letter, current_idx, k)
        num_steps += 1
    return partial_suffix[current_idx] + num_steps


def GetCountFromPartialCountMatrix(bw_word: str, partial_count: dict[str, list[int]], letter: str, i: int, k: int):
    """
    Determines the count of a character using a partial count matrix. Count being the number of times a specific
    letter appears in a word before and at a given index.
    Input:
        bw_word: the word we are using to count
        partial_count: The partial count matrix
        letter: The letter we are counting
        i: The index to count until
        k: the increment size of the partial count matrix
    Output:
        The number of times a letter appears in a word from the beginning of the word to the given index.
    """
    count = partial_count[letter][i // k]
    for x in range(i // k * k, i):
        if bw_word[x] == letter:
            count += 1
    return count


if __name__ == "__main__":
    main()
