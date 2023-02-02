import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from GenomePath import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "Files\\inputs\\Reconstruction")
    FilePaths = glob(dirpath + "\\input*.txt")
    for Path in FilePaths:
        Kmers = ReadTests_DeBruijnKmers(Path)
        AdjacencyDict = DeBruijnKmers(Kmers)

def StringReconstruction(Patterns: list[str]) -> str:
    dB = DeBruijnKmers(Patterns)
    path = EulerianPath(dB)
    Text = GenomePath(path)
    return Text

def ReadTest_StringReconstruction(filepath: str) -> tuple[int, list[str]]:
    with open(filepath) as f:
        for index, line in enumerate(f):
            if index == 0:
                k = line.strip()
            else:
                patterns = line.strip().split(" ")
    return k, patterns

if __name__ == "__main__":
    main()