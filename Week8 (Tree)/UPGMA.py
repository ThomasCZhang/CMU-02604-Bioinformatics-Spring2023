import os
import numpy as np


def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "UPGMA", "dataset_876248_8.txt")
    n, distance_matrix = ReadData_UPGMA(filepath)
    answer = UPGMA(distance_matrix, n)
    
    answer_path = os.path.join(dirpath, "answer.txt")
    with open(answer_path, "w") as f:
        for idx0, key in enumerate(answer):
            for idx1, edge in enumerate(answer[key]):
                if idx0 != 0 or idx1 != 0:
                    f.write("\n")
                f.write(str(key) + "->")
                f.write(str(edge[0]) + ":" + str(edge[1]))


def ReadData_UPGMA(filepath):
    with open(filepath) as f:
        n = int(f.readline().strip())
        distance_matrix = []
        for row in f:
            row = row.strip().split()
            row = [float(x) for x in row]
            distance_matrix.append(row)
        distance_matrix = np.array(distance_matrix)

    return n, distance_matrix


def UPGMA(distance_matrix: np.ndarray, n: int) -> dict[int, list[int]]:
    """
    Creates a phylogentic tree using UPGMA (Unweighted Pair Group Method and Arithmetic Mean) based on a distance
    matrix.
    Input:
        distance_matrix: The given distance matrix.
        n: The number of leaves in the phylogentic tree.
    Output:
        The generated phylogenetic tree as a adjacency dictionary.
    """
    tree = {}  # Empty dictionary to initialize tree
    for i in range(n):
        tree[i] = []
    age_dict = {}
    for i in range(n):
        age_dict[i] = 0
    
    clusters = [{i} for i in range(n)]
    node_number = [i for i in range(n)]  # Tracker of what node a row/column in distance_matrix corresponds to.
    
    iterations = -1
    while len(clusters) > 1:
        iterations += 1
        c_i, c_j = FindClosestClusters(distance_matrix)
        
        UpdateTree(tree, age_dict, node_number[c_i], node_number[c_j], n+iterations, distance_matrix[c_i, c_j]/2)        

        distance_matrix = UpdateDistanceMatrix(distance_matrix, c_i, c_j, clusters)
        for index in sorted([c_i, c_j], reverse=True):
            del node_number[index]
        node_number.append(n + iterations)
        
        clusters = MergeClusters(clusters, c_i, c_j)
    
    return tree


def UpdateTree(
    adjacency_dict: dict[int, list[int]],
    age_dict: dict[int, float],
    vertex1: int,
    vertex2: int,
    new_vertex: int,
    new_vertex_age: float,
):
    """
    Updates a tree represented by an adjacency dictionary and an age dictionary.
    Input:
        adjacency_dict: An adjacency dictionary representing the edges and vertices of a tree.
        age_dict: The age of each node.
        vertex1, vertex2: An edge will be added between these two verticies and the new vertex.
        new_vertex: The new vertex.
        new_vertex_age: the age of the new vertex.
    """
    AddWeightedEdge(adjacency_dict, vertex1, new_vertex, new_vertex_age - age_dict[vertex1])
    AddWeightedEdge(adjacency_dict, vertex2, new_vertex, new_vertex_age - age_dict[vertex2])
    age_dict[new_vertex] = new_vertex_age


def FindClosestClusters(distance_matrix: np.ndarray) -> np.ndarray:
    """
    Finds the closest clusters based on a distance matrix.
    Input:
        distance_matrix: The distance matrix.
    Output:
        The indicies representing the two clusters which are closest.
    """
    np.fill_diagonal(distance_matrix, np.inf)
    min_index = np.unravel_index(np.argmin(distance_matrix, axis=None), distance_matrix.shape)
    np.fill_diagonal(distance_matrix, 0)
    return min_index


def AddWeightedEdge(adjacency_dict: dict[int, list[int]], vertex1: int, vertex2: int, edge_weight: int):
    """
    AddEdgeToTree: Adds an edge to a tree.
    Input:
        adjacency_dict: The tree represented as an adjacency dictionary.
        vertex1, vertex2: The two vertices of the edge being added.
        edge_weight: The weight of the edge being added.
    """
    if vertex1 in adjacency_dict:
        adjacency_dict[vertex1].append([vertex2, edge_weight])
    else:
        adjacency_dict[vertex1] = [[vertex2, edge_weight]]

    if vertex2 in adjacency_dict:
        adjacency_dict[vertex2].append([vertex1, edge_weight])
    else:
        adjacency_dict[vertex2] = [[vertex1, edge_weight]]


def UpdateDistanceMatrix(
    distance_matrix: np.ndarray, cluster1: int, cluster2: int, clusters: list[set[int]]
) -> tuple[np.ndarray, list[set[int]]]:
    """
    UpdateDistanceMatrix: Updates a distance matrix by combining two "clusters" into one and updating the distances
    between the clusters.
    Inputs:
        distance_matrix: A matrix that contains the distance between each cluster.
        cluster1: The first cluster being combined.
        cluster2: The second cluster being combined.
        cluster_weights: The weights of the clusters.
    Output:
        A new distance matrix and a list of clusters.
    """
    new_distance_matrix = np.zeros([distance_matrix.shape[0] - 1, distance_matrix.shape[1] - 1])
    temp_matrix = np.delete(distance_matrix, [cluster1, cluster2], axis=0)
    temp_matrix = np.delete(temp_matrix, [cluster1, cluster2], axis=1)
    new_distance_matrix[:-1, :-1] = temp_matrix

    CalculateDistancesToNewCluster(new_distance_matrix, distance_matrix, clusters, cluster1, cluster2)
    return new_distance_matrix


def MergeClusters(clusters: list[set[int]], cluster1: int, cluster2: int):
    """
    Deletes two clusters in a distance matrix into one.
    Inputs:
        distance_matrix: distance matrix between the clusters.
        r1: the first cluster
        r2: the second cluster
        clusters: A list of sets represeting the clusters.
    """
    new_clusters = clusters.copy()
    if cluster1 > cluster2:
        del new_clusters[cluster1]
        del new_clusters[cluster2]
    else:
        del new_clusters[cluster2]
        del new_clusters[cluster1]
    new_clusters.append(clusters[cluster1].union(clusters[cluster2]))

    return new_clusters


def CalculateDistancesToNewCluster(
    new_distance_matrix: np.ndarray,
    old_distance_matrix: np.ndarray,
    clusters: list[set[int]],
    cluster1: int,
    cluster2: int,
):
    """
    Calculates the distances to the new cluster formed from combining two clusters and updates the distance matrix.
    Input:
        new_distance_matrix: The new distance matrix, with non-updated distance values to the new cluster.
        old_distance_matrix: The old distance matrix.
        clusters: The clusters.
        cluster1: The first cluster that was merged to create the new distance matrix.
        cluster2: The second cluster that was merged to create the new distance matrix.
    """
    j = 0  # Counter for new matrix because cluster positions change after delete two rows/cols.
    for i in range(old_distance_matrix.shape[0]):
        if i != cluster1 and i != cluster2:
            new_distance_matrix[j, -1] = CacluateWeightedAverageDistance(
                old_distance_matrix[i, cluster1],
                old_distance_matrix[i, cluster2],
                len(clusters[cluster1]),
                len(clusters[cluster2]),
            )
            new_distance_matrix[-1, j] = new_distance_matrix[j, -1]
            j += 1


def CacluateWeightedAverageDistance(distance1: float, distance2: float, weight1: int, weight2: int) -> float:
    """
    Calculates the weighted average distance.
    Inputs:
        distance1, distance2: The distances being averaged.
        weight1, weight2: The weights to assign each distance.
    Output:
        weighted average distance.
    """
    weighted_average = (weight1 * distance1 + weight2 * distance2) / (weight1 + weight2)
    return weighted_average


if __name__ == "__main__":
    main()
