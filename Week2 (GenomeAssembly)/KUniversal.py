import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from GenomePath import *
from StringReconstruction import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "Files\\inputs\\KUniverse")
    FilePaths = glob(dirpath + "\\data*.txt")
    for Path in FilePaths:
        k= ReadTest_K_Universe(Path)
        Text = KUniversalString(k)
    
    answer_path = os.path.join(os.path.dirname(__file__), "KUniverse.txt")
    with open(answer_path, 'w') as f:
        f.write(Text)


def ReadTest_K_Universe(filepath: str) -> int:
    with open(filepath) as f:
        k = int(f.readline().strip())
    return k

def KUniversalString(k: int) -> str:
    """
    Generates the K-universal string.

    Input:
        k : The length of each kmer.
    
    Output:
        Text: the universal string.
    """
    patterns = GenerateBinaryKMers(k)
    Text = StringReconstruction(patterns)
    Text = Text[:len(Text)+1-k]
    return Text

def GenerateBinaryKMers(K: int) -> list[str]:
    """
    GenerateKMers generates all binary kmers of length k.

    Input:
        K : The length of each k-mer.
    
    Output:
        kmers : The list of k-mers.
    """
    kmers = []
    for i in range(2**K):
        kmers.append(f"{i:0{K}b}")
    return kmers

if __name__ == "__main__":
    main()