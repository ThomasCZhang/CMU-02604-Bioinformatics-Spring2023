import os

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "SharedKMers", "input_0.txt")
    k, sequence1, sequence2 = ReadSequences(filepath)
    answer = SharedKMers(sequence1, sequence2, k)
    
    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, "w") as f:
        f.write(str(answer))

def ReadSequences(filepath):
    with open(filepath) as f:
        n = int(f.readline().strip())
        sequence1 = f.readline().strip()
        sequence2 = f.readline().strip()
    return n, sequence1, sequence2

def SharedKMers(sequence1: str, sequence2: str, k: int) -> int:
    """
    Counts the number of shared k-mers between two strings.
    """
    reverse_complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    sequence2_reverse_complement = ""
    for letter in sequence2[::-1]:
        sequence2_reverse_complement = "".join([sequence2_reverse_complement, reverse_complement_dict[letter]])

    forward_count = CountSharedKmers(sequence1, sequence2, k)
    backward_count = CountSharedKmers(sequence1, sequence2_reverse_complement, k)

    return forward_count + backward_count

def CountSharedKmers(sequence1: str, sequence2: str, k: int):
    count = 0
    print(sequence1)
    for i in range(len(sequence1)-k+1):
        kmer = sequence1[i:i+k]
        print(kmer)
        for j in range(len(sequence2)-k+1):
            if kmer == sequence2[j:j+k]:
                count += 1

    return count

if __name__ == "__main__":
    main()