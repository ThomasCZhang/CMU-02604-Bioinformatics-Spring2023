import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from GenomePath import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "Files\\inputs\\Reconstruction")
    FilePaths = glob(dirpath + "\\data*.txt")
    for Path in FilePaths:
        k, patterns = ReadTest_StringReconstruction(Path)
        Text = StringReconstruction(patterns)
    
    answer_path = os.path.join(os.path.dirname(__file__), "Reconstruction.txt")
    with open(answer_path, 'w') as f:
        f.write(Text)
    

def StringReconstruction(Patterns: list[str]) -> str:
    """
    StringReconstruction reconstructs a string from a list of k-mers.

    Input:
        Patterns: The list of k-mers.
    
    Output:
        Text: The reconstructed string.
    """
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