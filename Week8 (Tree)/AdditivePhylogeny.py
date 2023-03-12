import os
import numpy as np
from LimbLength import LimbLength
from LeafDistance import DFTraversal   

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "AdditivePhylogeny", "dataset_876246_6.txt")
    num_leaves, distance_matrix = ReadDistanceMatrix_AdditivePhylogeny(filepath)
    answer = AdditivePhylogeny(distance_matrix)
    
    answer_path = os.path.join(dirpath, "answer.txt")
    with open(answer_path, 'w') as f:
        for idx0, key in enumerate(answer):
            for idx1, edge in enumerate(answer[key]):
                if idx0 != 0 or idx1 != 0:
                    f.write('\n')
                f.write(str(key)+"->")
                f.write(str(edge[0])+":"+str(edge[1]))

def ReadDistanceMatrix_AdditivePhylogeny(filepath):
    with open(filepath) as f:
        num_leaves = int(f.readline().strip())
        distance_matrix = []
        for row in f:
            row = row.strip().split()
            row = [int(x) for x in row]
            distance_matrix.append(row)
        distance_matrix = np.array(distance_matrix)
    return num_leaves, distance_matrix

def AdditivePhylogeny(distance_matrix: np.ndarray[int], adjacency_dict: dict[int, list[int]] = None) -> dict[int, list[int]]:
# AdditivePhylogeny(D)
#     n ← number of rows in D
#     if n = 2
#         return the tree consisting of a single edge of length D1,2
#     limbLength ← Limb(D, n)
#     for j ← 1 to n - 1
#         Dj,n ← Dj,n - limbLength
#         Dn,j ← Dj,n
#     (i, k) ← two leaves such that Di,k = Di,n + Dn,k
#     x ← Di,n
#     D ← D﻿ with row n and column n removed
#     T ← AdditivePhylogeny(D)
#     v ← the (potentially new) node in T at distance x from i on the path between i and k
#     add leaf n back to T by creating a limb (v, n) of length limbLength
#     return T
    n = distance_matrix.shape[0]
    if adjacency_dict == None:
        adjacency_dict = dict()
        for i in range(n):
            adjacency_dict[i] = []

    if n == 2:
        weight = distance_matrix[0,1]
        AddEdge(adjacency_dict, 0, 1, weight)
        return adjacency_dict
    
    # Creating "bald" distance matrix
    limb = LimbLength(n-1, distance_matrix)
    for j in range(n-1):
        distance_matrix[j, n-1] = distance_matrix[j, n-1]-limb
        distance_matrix[n-1, j] = distance_matrix[j, n-1]

    found_leaves = False
    for i in range(n-1):
        for j in range(i+1, n-1):
            if distance_matrix[i, j] == distance_matrix[i, n-1] + distance_matrix[j, n-1]:
                leaf1 = i
                leaf2 = j
                found_leaves == True
                break
        if found_leaves:
            break
    
    leaf1_to_new_vertex = distance_matrix[leaf1, n-1]
    distance_matrix = distance_matrix[0:n-1, 0:n-1]
    adjacency_dict = AdditivePhylogeny(distance_matrix, adjacency_dict)
    # if len(adjacency_dict) == 5:
    #     print(f"Tree: {adjacency_dict}")
    #     print(f"Leaf 1: {leaf1}")
    #     print(f"Leaf to vertex: {leaf1_to_new_vertex}")
    path = GetPath(leaf1, leaf2, adjacency_dict)
    AddLeafToSimpleTree(adjacency_dict, n-1, len(adjacency_dict), leaf1, leaf1_to_new_vertex, limb, path)
    # if len(adjacency_dict) == 6:
    #     print(f"Tree: {adjacency_dict}")
    return adjacency_dict

def AddLeafToSimpleTree(adjacency_dict: dict[int, list[list[int]]], new_leaf: int, new_internal_node: int,
                                   vertex1: int, dist_to_new_vertex: int, limb_length: int, path: list[int]):
    """
    Adds a new leaf at a specified distance from a starting vertex along some path.
    """
    current_vertex = vertex1
    idx0 = 0
    current_dist = 0
    while current_dist < dist_to_new_vertex:
        idx0 += 1
        next_vertex = path[idx0]
        # Loop to find the edge corresponding to the next vertex.
        for edge in adjacency_dict[current_vertex]:
            if edge[0] == next_vertex:
                next_dist = current_dist + edge[1]
                break


        if next_dist == dist_to_new_vertex:  # New leaf is a child of existing vertex.
            AddEdge(adjacency_dict, next_vertex, new_leaf, limb_length)
            break
        elif next_dist > dist_to_new_vertex: # New leaf is not child of existing vertex
            # Create new internal node with leaf as one of its neighbors.
            AddEdge(adjacency_dict, new_internal_node, new_leaf, limb_length)
            # Add new internal node to adjacency list of current vertex.
            AddEdge(adjacency_dict, current_vertex, new_internal_node, dist_to_new_vertex-current_dist)
            # Add new internal node to adjacency list of next vertex.
            AddEdge(adjacency_dict, next_vertex, new_internal_node, next_dist-dist_to_new_vertex)
            # Delete the edge (current_vertex, next_vertex)
            DeleteEdge(adjacency_dict, current_vertex, next_vertex, next_dist-current_dist)
            break
        else:
            current_vertex = next_vertex
            current_dist = next_dist

def AddEdge(adjacency_dict: dict[int, list[int]], v1: int, v2: int, weight: int):
    """
    Adds an edge to an adjacency dictionary.
    Input:
        adjacency_dict: Adjacency dictionary to represent a tree. Has weighted edges.
        v1, v2: The vertices of the edge to add.
        weight: The weight of the edge to add
    """
    if v1 in adjacency_dict:
        adjacency_dict[v1].append([v2, weight])
    else:
        adjacency_dict[v1] = [[v2, weight]]
    if v2 in adjacency_dict:
        adjacency_dict[v2].append([v1, weight])
    else:
        adjacency_dict[v2] = [[v1, weight]]

def DeleteEdge(adjacency_dict: dict[int, list[int]], v1: int, v2: int, weight: int=None):
    """
    Deletes an edge from an adjacency dictionary.
    Input:
        adjacency_dict: Adjacency dictionary to represent a tree. Has weighted edges.
        v1, v2: The vertices of the edge to remove.
        weight: The weight of the edge being removed. If None, this will be determined by searching through the 
        adjacency dictionary for the existing edge.
    """
    if weight == None:
        for edge in adjacency_dict[v1]:
            if edge[0] == v2:
                weight = edge[1]
                break
    adjacency_dict[v1].remove([v2, weight])
    adjacency_dict[v2].remove([v1, weight])

def GetPath(vertex1: int, vertex2: int, adjacency_dict: dict[int, list[list[int]]]):
    """
    Finds the path between two verticies in a tree.
    Input:
        vertex1: The starting vertex.
        vertex2: The ending vertex.
        adjacency_dict: The tree, represented as an adjacency dictionary.
    Output:
        The path from vertex1 to vertex2 in the tree.
    """
    visited = set()
    path = dict()
    current_path = list()
    
    current_vertex = vertex1
    weight = dict([(current_vertex, 0)])
    DFTraversal(current_vertex, adjacency_dict, visited, weight, path, current_path)
    v1_v2_path = path[vertex2]
    return v1_v2_path

if __name__ == "__main__":
    main()