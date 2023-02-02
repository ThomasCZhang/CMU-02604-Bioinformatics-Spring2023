import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from GenomePath import *
from StringReconstruction import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "Files\\inputs\\KUniverse")
    FilePaths = glob(dirpath + "\\input*.txt")
    for Path in FilePaths:
        k= ReadTest_K_Universe(Path)
        patterns = GenerateKMers(k)
        Text = StringReconstruction(patterns)
    
    answer_path = os.path.join(os.path.dirname(__file__), "KUniverse.txt")
    with open(answer_path, 'w') as f:
        f.write(Text)


def ReadTest_K_Universe(filepath: str) -> int:
    with open(filepath) as f:
        k = f.readline().strip()
    return k

def GenerateKMers(K: int) -> list[str]:
    kmers = []
    for i in range(2**K):
        kmers.append(f"{i:0{K}b}")
    return kmers

if __name__ == "__main__":
    main()