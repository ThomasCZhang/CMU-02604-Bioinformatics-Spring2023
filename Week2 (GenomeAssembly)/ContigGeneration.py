import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from NonBranchingPath import *
from GenomePath import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "Files\\inputs\\ContigGeneration")
    FilePaths = glob(dirpath + "\\data*.txt")
    for Path in FilePaths:
        patterns = ReadTests_DeBruijnKmers(Path)
        answer = ContigGeneration(patterns)

    answer_path = os.path.join(os.path.dirname(__file__), "ContigGeneration.txt")
    with open(answer_path, 'w') as f:
        for idx, contig in enumerate(answer):
            if idx != 0:
                f.write(" ")
            f.write(contig)

def ContigGeneration(Patterns: list[str]) -> list[str]:
    """
    ContigGeneration finds the contigs from the DeBruijn graph generated from a collection of k-mers.

    Input:
        Patterns : The collection of K-mers.

    Output:
        contigs : The list of contigs.
    """
    G = DeBruijnKmers(Patterns)
    contig_paths = MaximalNonBranchingPaths(G)
    contigs = []
    for contig_path in contig_paths:
        contig = GenomePath(contig_path)
        contigs.append(contig)

    return contigs

if __name__ == "__main__":
    main()