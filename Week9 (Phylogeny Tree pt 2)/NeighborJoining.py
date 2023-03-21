import os
import numpy as np

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "NeighborJoining", "dataset_876249_7.txt")
    num_leaves, distance_matrix = ReadData_NeighborJoining(filepath)
    answer = NeighborJoining(distance_matrix)

    answer_path = os.path.join(dirpath, "answer.txt")
    with open(answer_path, "w") as f:
        for idx0, key in enumerate(answer):
            for idx1, val in enumerate(answer[key]):
                if idx0 != 0 or idx1 != 0:
                    f.write("\n")
                f.write(str(key) + "->" + str(val[0]) + ":" + str(val[1]))


def ReadData_NeighborJoining(filepath):
    with open(filepath) as f:
        num_leaves = int(f.readline().strip())
        distance_matrix = []
        for row in f:
            row = row.strip().split()
            row = [float(x) for x in row]
            distance_matrix.append(row)
        distance_matrix = np.array(distance_matrix)
    return num_leaves, distance_matrix


def NeighborJoining(
    distance_matrix: np.ndarray, adjacency_dict: dict[int, list[int]] = None, vertex_number: list[int] = None
) -> dict[int, list[int]]:
    """
    Performs the neighbor joining algorithm.
    Input:
        distance_matrix: The distance between verticies that can be clustered
        adjacency_dict: The tree represented as an adjacency dictionary.
        vertex_number: The vertex that each row in distance matrix corresponds to.

    """
    n = distance_matrix.shape[0]

    if adjacency_dict == None:
        vertex_number = [i for i in range(n)]
        adjacency_dict = {}
        for i in range(n):
            adjacency_dict[i] = []
    if n == 2:
        AddEdge(adjacency_dict, vertex_number[0], vertex_number[1],distance_matrix[0, 1])
        return distance_matrix

    d_star = ConstructNeighborJoiningMatrix(distance_matrix)
    # Row, col indicies of the minimum value in the matrix.
    row_of_min, col_of_min = FindMatrixMinimumOffDiagonal(d_star)
    
    UpdateTree(adjacency_dict, distance_matrix, len(adjacency_dict), row_of_min, col_of_min, vertex_number)
    distance_matrix = UpdateDistanceMatrix(distance_matrix, row_of_min, col_of_min)
    UpdateVertexNumber(vertex_number, row_of_min, col_of_min)
    NeighborJoining(distance_matrix, adjacency_dict, vertex_number)
    return adjacency_dict


def UpdateVertexNumber(vertex_number: list[int], i: int, j: int):
    """
    Updates the list of vertex numbers by deleting two values and appending a new one. The new value will be
    max(list of vertex numbers) + 1.
    Input:
        i, j: The two positions to delete.
    """

    vertex_number.append(max(vertex_number) + 1)
    if i > j:
        del vertex_number[i]
        del vertex_number[j]
    else:
        del vertex_number[j]
        del vertex_number[i]


def FindMatrixMinimumOffDiagonal(m: np.ndarray) -> tuple[int]:
    """
    Finds the minimum element of a matrix that is not on the diagonal
    Input:
        m: the matrix
    """
    np.fill_diagonal(m, np.inf)
    indexes = np.unravel_index(np.argmin(m, axis=None), m.shape)
    np.fill_diagonal(m, 0)
    return indexes


def UpdateDistanceMatrix(distance_matrix: np.ndarray, i: int, j: int) -> np.ndarray:
    """
    Updates a leaf distance matrix by combining neighboring leaves.
    Input:
        distance_matrix: Old distance matrix
        i, j: The indicies of the two rows/cols that will be combined.
    Output:
        new distance matrix.
    """
    new_distance_matrix = np.zeros([distance_matrix.shape[0] + 1, distance_matrix.shape[1] + 1])
    new_distance_matrix[:-1, :-1] = distance_matrix
    for w in range(distance_matrix.shape[0]):
        if w != i and w != j:
            new_distance_matrix[w, -1] = 0.5 * (distance_matrix[i, w] + distance_matrix[j, w] - distance_matrix[i, j])
            new_distance_matrix[-1, w] = new_distance_matrix[w, -1]
    new_distance_matrix = np.delete(new_distance_matrix, [i, j], axis=0)
    new_distance_matrix = np.delete(new_distance_matrix, [i, j], axis=1)
    return new_distance_matrix


def ConstructNeighborJoiningMatrix(distance_matrix: np.ndarray) -> np.ndarray:
    """
    Creates the neighbor-joining matrix d*.
    Input:
        distance_matrix: The distance matrix from which to create d*
    Output:
        The neighbor-joining matrix d*
    """
    n = distance_matrix.shape[0]
    d_star = np.zeros(distance_matrix.shape)
    for i in range(n):
        for j in range(i + 1, n):
            d_star[i, j] = (
                (n - 2) * distance_matrix[i, j] - np.sum(distance_matrix[i, :]) - np.sum(distance_matrix[j, :])
            )
    return d_star


def UpdateTree(
    adjacency_dict: dict[int, list[int]],
    distance_matrix: np.ndarray,
    new_vertex: int,
    row1: int,
    row2: int,
    vertex_number: list[int],
):
    """
    Updates a tree using neighbor joining algorithm. Adds a vertex "new_vertex" to the tree. Then connects vertex1 and
    vertex2 to new_vertex. Weight is determined by limb length formula.
    Input:
        adjacency_dict: a graph represented as an adjacency dictionary.
        distance_matrix: a distance matrix between the vertices that can be connected via a parent vertex.
        new_vertex: the new vertex being added to the tree.
        row1, row2: the two rows in the distance matrix being joined.
        vertex_number: the vertex number corresponding to a specific row in the distance matrix.
    """

    delta_ij = CalculateDelta(distance_matrix, row1, row2)
    vertex1 = vertex_number[row1]
    vertex2 = vertex_number[row2]
    AddEdge(adjacency_dict, vertex1, new_vertex, 0.5 * (distance_matrix[row1, row2] + delta_ij))
    AddEdge(adjacency_dict, vertex2, new_vertex, 0.5 * (distance_matrix[row1, row2] - delta_ij))


def CalculateDelta(distance_matrix: np.ndarray, i: int, j: int):
    """
    Calculates (TotalDistance(i) - TotalDistance(j))/(n-2)
    Input:
        Distance_matrix: A distance matrix.
        i, j: The rows between which to calculate TotalDistance(i) - TotalDistance(j)
    Output:
        The value of delta_ij.
    """
    total_dist_i = np.sum(distance_matrix[i, :])
    total_dist_j = np.sum(distance_matrix[j, :])
    n = distance_matrix.shape[0]
    return (total_dist_i - total_dist_j) / (n - 2)


def AddEdge(adjacency_dict: dict[int, list[int]], vertex1: int, vertex2: int, weight: float):
    """
    Adds an undirected edge to a graph. The graph is represented by an adjacency dictionary.
    Input:
        adjacency_dict: the adjacency dictionary representation of a graph.
        vertex1: the first vertex of the edge being added.
        vertex2: the second vertex of the edge being added.
        weight: the weight of the edge.
    """
    if vertex1 in adjacency_dict:
        adjacency_dict[vertex1].append([vertex2, weight])
    else:
        adjacency_dict[vertex1] = [[vertex2, weight]]
    if vertex2 in adjacency_dict:
        adjacency_dict[vertex2].append([vertex1, weight])
    else:
        adjacency_dict[vertex2] = [[vertex1, weight]]


def DeleteEdge(adjacency_dict: dict[int, list[int]], v1: int, v2: int, weight: int = None):
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
    if weight != None:  # If weight still = none, then the edge v1-v2 is not in the graph.
        adjacency_dict[v1].remove([v2, weight])
        adjacency_dict[v2].remove([v1, weight])


if __name__ == "__main__":
    main()
