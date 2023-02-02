import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\inputs\\EulerianPath")
    filepaths = glob(dirpath + "\\input_1.txt")
    G = ReadTests_Eulerian(filepaths[0])
    # eulerian_cycle = EulerianCycle(G)
    # eulerian_path = EulerianPath(G)
    test = [1,2,3]
    print(test[2:-1]+test[:2])
    
    # answer_path = os.path.join(os.path.dirname(__file__), "eulerian_answer.txt")
    # with open(answer_path, 'w') as f:
    #     for ind, ele in enumerate(eulerian_cycle):
    #         if ind == 0:
    #             f.write(str(ele))
    #         else:
    #             f.write(" " + str(ele))

def ReadTests_Eulerian(FilePath: str) -> dict[int, int]:
    EdgeDict = {}
    with open(FilePath) as f:
        for _, line in enumerate(f):
            cur_line = line.strip().split(" ")
            Key = int(cur_line[0][:-1])
            Val = [int(x) for x in cur_line[1:]]
            EdgeDict[Key] = Val 
    return EdgeDict

def EulerianCycle(G: dict[int, list[int]]) -> list[int]:
    """
    Takes a graph containing an Eulerian cycle and returns the path of a eulerian cycle.
    G (dict[int, list[int]]): A dictionary of edges. Keys are a vertex. Values are list of vertices that edge connects 
    to.
    """
    num_edges = 0
    for key in G:
        num_edges += len(G[key])
    path = [list(G)[0]] # Arbitrary start vertex (since graph contains eulerian cycle).
    while num_edges > 0:
        last_vertex = path[len(path)-1] # last vertex visited.
        if not G[last_vertex]: # Check if last_vertex has any unvisited neighbors.
            # Need to choose another starting position
            for index, vertex in enumerate(path):
                if G[vertex]:
                    path = path[index:-1]+path[:index]+[vertex]
                    break
        else:
            # Keep traversing
            path.append(G[last_vertex][0])
            G[last_vertex] = G[last_vertex][1:] # Remove edge we've traversed on.
            num_edges -= 1
    return path

def EulerianPath(G: dict[int, list[int]]) -> list[int]:
    in_out_degree = InOutDegree(G)
    num_unbalanced = 0
    for key in in_out_degree:
        if in_out_degree[key] == -1:
            start_vertex = key
            num_unbalanced += 1
        elif in_out_degree[key] == 1:
            end_vertex = key
            num_unbalanced += 1
    if num_unbalanced > 2:
        print("Warning: Graph does not have Eulerian Path.")
    
    if end_vertex in G:
        G[end_vertex].append(start_vertex)
    else:
        G[end_vertex] = [start_vertex]

    eulerian_cycle = EulerianCycle(G)
    for index, val in enumerate(eulerian_cycle):
        next_val = eulerian_cycle[index+1]
        if (val == end_vertex) and (next_val == start_vertex):
            pass

        
def InOutDegree(G: dict[int, list[int]]) -> dict[int, int]:
    degrees = {}
    for key in G:
        for val in G[key]:
            if key in degrees:
                degrees[key] -= 1
            else:
                degrees[key] = -1
            if val in degrees:
                degrees[val] += 1
            else:
                degrees[val] = 1
    return degrees


if __name__ == "__main__":
    main()