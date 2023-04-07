import os

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "bwt_decode", "input_0.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    answer = bwt_decode(word)
    
    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, "w") as f:
        f.write(answer)


def bwt_decode(word: str) -> str:
    """
    Decodes a Burrows-Wheeler encoded word.
    Input:
        word: The Burrows-Wheeler encoded word.
    """
    ranks = get_ranks(word)
    idx_dict = get_char_positions(word)
    decoded = "$"
    next_idx = 0
    for _ in range(len(word)-1):
        next_char = word[next_idx]
        decoded = "".join([next_char, decoded])
        next_idx = ranks[next_char] + idx_dict[next_char][next_idx]
   
    return decoded
    

def get_char_positions(word: str) -> dict[str, dict[int, int]]:
    """
    gets the indicies that each unique character in a word appears at.
    Input:
        word: The word being analyzed
    """
    char_idx_dict = {}
    for i, letter in enumerate(word):
        if letter in char_idx_dict:
            order = len(char_idx_dict[letter])
            char_idx_dict[letter][i] = order
        else:
            char_idx_dict[letter] = {i: 0}
    return char_idx_dict

def get_ranks(word: str) -> dict[str, int]:
    """
    Determines the rank of each character in word. The rank of a character is the number of times characters that are 
    "smaller" than that character appear in the word. A "smaller" character, is a character that occurs earlier
    when ordered lexicographically.
    Input:
        word: the word being analyzed.
    Output:
        The rank of each unique character in word stored in a dictionary.
    """
    letter_counts = count_letters(word)
    sorted_chars = sorted(letter_counts.keys())
    rank = {}
    total = 0
    for char in sorted_chars:
        rank[char] = total
        total += letter_counts[char]
    return rank

def count_letters(word: str) -> dict[str, int]:
    """
    Counts the number of times a character occurs in word.
    Input:
        word: the word being analyzed.
    Output:
        A dictionary were the key is the character being analyzed, and the value is the number of times that letter
        occurs in word.
    """
    count_dict = {}
    for character in word:
        if character in count_dict:
            count_dict[character] += 1
        else:
            count_dict[character] = 1
    return count_dict


if __name__ == "__main__":
    main()
