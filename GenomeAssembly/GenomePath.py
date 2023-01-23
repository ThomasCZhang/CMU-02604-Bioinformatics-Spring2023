import os
from glob import glob

def main():
    DirPath = os.path.join(os.path.dirname(__file__), "Files\\Inputs\\GenomePath")
    FilePaths = glob(DirPath + "\\*.txt")

    for path in FilePaths:
        KMer_Sequences = ReadTestFiles_GenomePath(path)
        FullSequence = GenomePath(KMer_Sequences)
        print(f"For file {os.path.basename(path)} the full DNA sequence is {FullSequence}.")

def GenomePath(path: list[str]) -> str:
    """
    Takes a sequence path of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Pattern_i are equal 
    to the first k-1 symbols of Pattern_i+1 for 1 ≤ i ≤ n-1. Returns a string of Text of length k+n-1 such that the 
    it-th k-mer in Text is equal to Pattern_i (for 1 ≤ i ≤ n).\n
    Input: \n
    path (list[str]): the list of patterns.\n
    Output:\n
    FinalString (str): the final string of text
    """
    k = len(path[0])
    FinalString = ""
    for i, string in enumerate(path):
        if i == 0:
            FinalString = FinalString + string
        else:
            FinalString = FinalString + string[k-1]
    return FinalString

def ReadTestFiles_GenomePath(FilePath: str) -> list[str]:
    """
    Reads the test files for GenomePath and returns the paths as a list of strings.\n
    Input:\n
    FilePath (str): The path to the test file.\n
    Output:\n
    Sequences (list[str]): The sequences as in a list of strings.
    """
    with open(FilePath) as f:
        Line = f.readline()
        Sequences = Line.split(" ")
    return Sequences

if __name__ == "__main__":
    main()