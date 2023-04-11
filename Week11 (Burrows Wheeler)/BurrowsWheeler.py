import os

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "bwt", "dataset_876289_5.txt")
    word = ReadData(filepath)
    answer = bwt(word)

    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, "w") as f:
        f.write(answer)

def ReadData(filepath):
    with open(filepath) as f:
        return f.readline().strip()

def bwt(word: str) -> str:
    """
    Returns the Burrows-Wheeler Transform on a word.
    Input:
        word: the word we are performing Burrow's Wheeler on.
    Output:
        The Burrows-Wheeler Tranform.
    """
    word_rotations = [word[i:]+word[:i] for i in range(len(word))]
    word_rotations = sorted(word_rotations)
    bwt = ""
    for w in word_rotations:
        bwt = "".join([bwt, w[-1]])
    return bwt

if __name__ == "__main__":
    main()