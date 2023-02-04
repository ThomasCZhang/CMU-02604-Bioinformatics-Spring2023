import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\inputs\\EulerianPath")
    filepaths = glob(dirpath + "\\data*.txt")
    G = ReadTests_Eulerian(filepaths[0])
    # eulerian_cycle = EulerianCycle(G)
    answer = EulerianPath(G)
    
    answer_path = os.path.join(os.path.dirname(__file__), "eulerian_answer.txt")
    with open(answer_path, 'w') as f:
        for ind, ele in enumerate(answer):
            if ind == 0:
                f.write(str(ele))
            else:
                f.write(" " + str(ele))

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
    
    Input:
        G: A dictionary of edges. Keys are a vertex. Values are list of vertices that edge connects to.
    
    Output:
        path: A list of vertices corresponding to the eulerian cycle path.
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
    """
    EulerianPath finds the eulerian path in a graph.
    
    Input:
        G : A graph containing a eulerian path. Represented as a dictionary of edges.

    Output:
        path: The eulerian path, as a list of vertices.
    """
    in_degree, out_degree = InOutDegree(G)
    in_out_degree = {}
    for key in in_degree:
        in_out_degree[key] = in_degree[key]
    for key in out_degree:
        if key in in_out_degree:
            in_out_degree[key] -= out_degree[key]
        else:
            in_out_degree[key] = -out_degree[key]

    # Finding the starting and ending vertex
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
    
    # Adding an edge so we make the graph contain eulerian cycle.
    if num_unbalanced == 2:
        if end_vertex in G:
            G[end_vertex].append(start_vertex)
        else:
            G[end_vertex] = [start_vertex]

    cycle = EulerianCycle(G)
    if num_unbalanced == 2:
        for index, val in enumerate(cycle[:-1]):
            next_val = cycle[index+1]
            if (val == end_vertex) and (next_val == start_vertex):
                path = cycle[index+1:-1]+cycle[:index+1]
                return path
    else:
        return cycle

    print("Something broke if this prints.")
        
def InOutDegree(G: dict[int, list[int]]) -> tuple[dict[int, int], dict[int, int]]:
    """
    InOutDegree: Determines the in and out degrees of all vertices in a graph.

    Input:
        G: A graph represnted as an adjacency list. The adjacency list is stored as a dictionary.

    Output:
        in_degree, out_degree: Two dictionaries that store the in and out degrees of each vertex in G.
    """
    in_degree = {}
    out_degree = {}
    for key in G:
        for val in G[key]:
            if key in out_degree:
                out_degree[key] += 1
            else:
                out_degree[key] = 1
            
            if val in in_degree:
                in_degree[val] += 1
            else:
                in_degree[val] = 1

            if key not in in_degree:
                in_degree[key] = 0
            if val not in out_degree:
                out_degree[val] = 0
                
    return in_degree, out_degree


if __name__ == "__main__":
    main()