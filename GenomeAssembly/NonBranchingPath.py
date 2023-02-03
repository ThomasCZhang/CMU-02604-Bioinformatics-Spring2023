import os
from glob import glob
from DeBruijnKmers import *
from Eulerian import *
from GenomePath import *
from StringReconstruction import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\inputs\\NonBranchingPath")
    filepaths = glob(dirpath + "\\data*.txt")
    G = ReadTests_Eulerian(filepaths[0])
    answer = MaximalNonBranchingPaths(G)
    
    answer_path = os.path.join(os.path.dirname(__file__), "NonBranchingPath.txt")
    with open(answer_path, 'w') as f:
        for idx1, traveled_path  in enumerate(answer):
            if idx1 != 0:
                f.write("\n")
            for idx2, vertex in enumerate(traveled_path):
                if idx2 == 0:
                    f.write(str(vertex))
                else:
                    f.write(" " + str(vertex))
            

def MaximalNonBranchingPaths(G: dict[int, list[int]]) -> list[list[int]]:
    traveled = {}
    in_degree, out_degree = InOutDegree(G)
    paths = []
    for vertex in G:
        if (in_degree[vertex] != 1) or (out_degree[vertex] != 1):
            if out_degree[vertex] > 0:
                traveled[vertex] = True
                for x in G[vertex]:
                    NonBranchingPath = [vertex, x]
                    traveled[x] = True
                    while (in_degree[x] == 1) and (out_degree[x] == 1):
                        NonBranchingPath.append(G[x][0])
                        traveled[G[x][0]] = True
                        x = G[x][0]
                    paths.append(NonBranchingPath)

    for vertex in G:
        if vertex not in traveled:
            if (in_degree[vertex] == 1) and (out_degree[vertex] == 1):
                NonBranchingPath = [vertex]
                traveled[vertex] = True
                
                x = G[vertex][0]
                while x not in traveled:
                    NonBranchingPath.append(x)
                    traveled[x] = True
                    x = G[x][0]
                NonBranchingPath.append(x)
                paths.append(NonBranchingPath)
    
    return paths

if __name__ == "__main__":
    main()