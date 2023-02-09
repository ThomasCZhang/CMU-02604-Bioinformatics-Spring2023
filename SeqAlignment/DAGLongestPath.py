# Insert your LongestPath function here, along with any subroutines you need
# s and t are the starting and the ending nodes of the path respectively
# E[u] is the list of neighbors of the vertex u, paired with corresponding edge weights
# LongestPath should return the length of the longest path betweeen s and t together with
# the list of nodes of the path
def LongestPath(s: int, t: int, E: dict[int, list[tuple[int, int]]]) -> tuple[int, list[int]]:\
    """
    LongestPath: Finds the longest path between a starting vertex and ending vertex in a DAG.
    
    Input:
        s: starting vertex.
        
        t: ending vertex.
        
        E: Assume that nodes are numbered by topological sort order in E.

    Output:
        The length of the longest path.
        
        The path between and s and t.
    """
    nodes = list(E)
    
    pass
