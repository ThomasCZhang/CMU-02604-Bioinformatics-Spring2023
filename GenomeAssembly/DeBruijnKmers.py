import os
from glob import glob


def main():
    DirPath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\DeBruijn")
    FilePaths = glob(DirPath + "\\data*.txt")
    for Path in FilePaths:
        Kmers = ReadTests_DeBruijnKmers(Path)
        AdjacencyDict = DeBruijnKmers(Kmers)
    
    AnswerPath = os.path.join(os.path.dirname(__file__), "DeBruijnAnswer.txt")
    with open(AnswerPath, "w") as f:
        for Key in AdjacencyDict:
            f.write(Key + " : ")
            for Text in AdjacencyDict[Key]:
                f.write(" "+Text)
            f.write("\n")

# Insert your DeBruijnKmers function here, along with any subroutines you need


def DeBruijnKmers(k_mers: list[str]) -> dict[str, list[str]]:
    """
    Creates an adjacency dict for the DebruijnKmers graph based on a list of k mers.
    """
    AdjacencyDict = {}
    for Kmer in k_mers:
        Key = Kmer[:-1]
        if Key in AdjacencyDict:
            AdjacencyDict[Kmer[:-1]].append(Kmer[1:])
        else:
            AdjacencyDict[Kmer[:-1]] = [Kmer[1:]]
    return AdjacencyDict


def ReadTests_DeBruijnKmers(FilePath: str) -> list[str]:
    with open(FilePath) as f:
        Line = f.readline()
        Line = Line.strip()
        Kmers = Line.split(" ")
    return Kmers


if __name__ == "__main__":
    main()
