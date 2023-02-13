import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\longestpath")
    filepaths = glob(dirpath + "\\data*.txt")
    for filepath in filepaths:
        s, e, graph = ReadTestInputs_LongestPath(filepath)
        answer = LongestPath(s, e, graph)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        for idx0, ele in enumerate(answer):
            if idx0 == 0:
                f.write(str(ele))
            else:
                f.write("\n")
                for idx1, val in enumerate(ele):
                    if idx1 != 0:
                        f.write(" ")
                    f.write(str(val))



def ReadTestInputs_LongestPath(FilePath: str) -> tuple[int, int, dict[int, list[tuple[int, int]]]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        N: the number of leader peptides to keep.

        spectrum: the mass spectrum.
    """
    graph = {}
    with open(FilePath) as f:
        for idx, line in enumerate(f):
            line = [int(x) for x in line.strip().split(" ")]
            if idx == 0:
                s = line[0]
                e = line[1]
            else:
                vertex = line[0]
                edge = line[1:3]
                if line[0] in graph:
                    graph[vertex].append(tuple(edge))
                else:
                    graph[vertex] = [tuple(edge)]

    return s, e, graph

# s and t are the starting and the ending nodes of the path respectively
# E[u] is the list of neighbors of the vertex u, paired with corresponding edge weights
# LongestPath should return the length of the longest path betweeen s and t together with
# the list of nodes of the path
def LongestPath(s: int, t: int,
                E: dict[int, list[tuple[int, int]]]) -> tuple[int, list[int]]:
    """
    LongestPath returns the longest path in a DAG between a starting and ending vertex

    Input:
        s: starting vertex.
        t: ending vertex
        E: graph as a dictionary of directed edges. 

    Output: 
        length of longest path, and path as a list nodes in the path.
    """
    start_vertices = sorted(list(E.keys()))
    # end_vertices = set([x[0] for edges in E.values() for x in edges])
    # all_vertices = sorted(list(start_vertices.union(end_vertices)))

    start_idx = start_vertices.index(s)
    end_idx = start_idx
    while start_vertices[end_idx] < t:
        end_idx += 1
        if end_idx >= len(start_vertices):
            break
    
    linked_list = {}
    best_score = {s: 0}
    for i in range(start_idx,end_idx):
        curr_vertex = start_vertices[i]
        if curr_vertex in best_score:
            for edge in E[curr_vertex]:
                next_vertex = edge[0]
                score = best_score[curr_vertex] + edge[1] 
                if next_vertex not in best_score:
                    best_score[next_vertex] = score
                    linked_list[next_vertex] = curr_vertex
                elif score > best_score[next_vertex]:
                    best_score[next_vertex] = score
                    linked_list[next_vertex] = curr_vertex

    path = [t]
    while path[0] != s:
        path = [linked_list[path[0]], *path]

    return best_score[t], path

if __name__ == "__main__":
    main()