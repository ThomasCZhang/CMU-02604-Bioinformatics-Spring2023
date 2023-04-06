import os

def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "suffixarray", "input_1.txt")
    s = ReadString(filepath)
    answer = SuffixArray(s)

    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, "w") as f:
        for i, val in enumerate(answer):
            if i != 0:
                f.write(' ')
            f.write(f"{val}")

def ReadString(filepath):
    with open(filepath) as f:
        return f.readline().strip()

def SuffixArray(s: str):
    """
    Gets the index order of the suffix array:
    """
    suffixes = [(i,s[i:]) for i in range(len(s))]
    suffixes = sorted(suffixes, key=lambda x: x[1])
    indexes = [x[0] for x in suffixes]
    return indexes

if __name__ == "__main__":
    main()