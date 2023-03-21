import os
import numpy as np

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath,"inputs", "LeafDistance", "dataset_876244_12.txt")
    num_leaves, adjacency_dict = ReadAdjacencyList(filepath)
    answer = LeafDistance(adjacency_dict)

    answerpath = os.path.join(dirpath, "answer.txt")
    np.savetxt(answerpath, answer,fmt="%d", delimiter="\t")

def ReadAdjacencyList(filepath):
    with open(filepath) as f:
        num_leaves = int(f.readline().strip())
        adjacency_dict = dict()
        for row in f:
            row = row.strip().split("->")
            key = int(row.pop(0))
            row = [int(x) for x in row[0].split(":")]
            if key in adjacency_dict:
                adjacency_dict[key].append(row)
            else:
                adjacency_dict[key] = [row]
    return num_leaves, adjacency_dict

def LeafDistance(adj_dict: dict[int, list[list[int]]]) -> np.ndarray[int]:
    """
    LeafDistance: Takes an tree represented as an adjacency dictionary and returns a matrix that contains the distance
    between the leaves of the tree.
    Input:
        adj_dict: The tree represented as an adjacency dictionary.
    Output:
        A matrix that contains the distances between leaves.
    """
    leaves = GetLeaves(adj_dict)
    visited = set()
    path = dict()
    current_path = list()
    
    current_vertex = leaves[0]
    weight = dict([(current_vertex, 0)])
    DFTraversal(current_vertex, adj_dict, visited, weight, path, current_path)
    
    distance_matrix = np.zeros((len(leaves), len(leaves)))
    for idx0 in range(len(leaves)):
        leaf0 = leaves[idx0]
        leaf0_path = set(path[leaf0])
        for idx1 in range(idx0+1, len(leaves)):
            leaf1 = leaves[idx1]
            leaf1_path = path[leaf1]
            for vertex in reversed(leaf1_path): # iterate through leaf1_path from the end.
                if vertex in leaf0_path:
                    lca = vertex
                    break
            leaf0_to_lca = abs(weight[lca]-weight[leaf0])
            leaf1_to_lca = abs(weight[lca]-weight[leaf1])
            leaf_to_leaf = leaf0_to_lca +leaf1_to_lca
            distance_matrix[idx0,idx1] = leaf_to_leaf
            distance_matrix[idx1,idx0] = leaf_to_leaf
    
    return distance_matrix

def DFTraversal(current_vertex: int, adj_dict: dict[int, list[list[int]]], visited: set[int], weight:dict[int, int],
                 path: dict[int, list[int]], current_path: list[int]):
    """
    DFTraversal: Traverses a tree and keeps track of the weight-distance from the starting point of the traversal
    as well as the path taken to the current vertex.
    Input:
        current_vertex: the current vertex. 
        adj_dict: An adjacency dictionary that represents the tree.
        visited: the visited vertices. 
        weight: a dictionary that keeps track of the weight distance from the start vertex of the traversal
        path: a dictionary for holding path of verticies to each leaf.
        current_path: a list of verticies indicating the path taken from the start vertex to the current vertex.
    """
    current_path = current_path.copy()
    current_path.append(current_vertex)
    if len(adj_dict[current_vertex]) == 1 and current_vertex not in visited: # At a leaf vertex.
        path[current_vertex] = current_path

    visited.add(current_vertex)
    for edge in adj_dict[current_vertex]:
        v, w = edge
        if v not in visited:
            weight[v] = weight[current_vertex] + w
            DFTraversal(v, adj_dict, visited, weight, path, current_path)

def GetLeaves(adj_dict: dict[int, list[list[int]]]) -> list[int]:
    """
    Finds the leaf nodes of a tree. Tree represented as an adjacency dict with edge weights.
    Input:
        adj_dict: The adjacency dictionary representation of a tree.
    Output:
        A list of keys of the adjacency dictionary that correspond to tree leaves.
    """
    leaves = []
    for key in adj_dict:
        if len(adj_dict[key]) == 1:
            leaves.append(key)
    return leaves


if __name__ == "__main__":
    main()