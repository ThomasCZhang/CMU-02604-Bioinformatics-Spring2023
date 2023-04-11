import os

def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "PartialSuffixArray", "dataset_876299_2.txt")
    word, k = ReadData(filepath)
    partial_suffix = PartialSuffixArray(word, k)
    
    answer_path = os.path.join(dirname, "answer.txt")
    with open(answer_path, 'w') as f:
        for key in partial_suffix:
            f.write(f'{key} {partial_suffix[key]}\n')
        
def ReadData(filepath):
    with open(filepath) as f:
        word = f.readline().strip()
        k = int(f.readline().strip())    
    return word, k

def PartialSuffixArray(s: str, k: int) -> dict[int, int]:
    """
    Generates a partial suffix array. Keeps every k'th element of the suffix array.
    Input:
        s: The word to create a partial suffix array of.
        k: The increment size for which elements of the suffix array to keep.
    Output:
        The partial suffix array as a list of integers.
    """
    suffixes = [(i,s[i:]) for i in range(len(s))]
    suffixes = sorted(suffixes, key=lambda x: x[1])
    indexes = [x[0] for x in suffixes]
    suffix_dict = {}
    for i, j in enumerate(indexes):
        if j%k == 0:
            suffix_dict[i] = j
    return suffix_dict

if __name__ == "__main__":
    main()