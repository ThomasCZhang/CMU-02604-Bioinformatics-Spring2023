import os
import numpy as np

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath,"inputs", "LimbLength", "input_1.txt")
    num_leaves,target_leaf,distance_matrix = ReadDistanceMatrix_LimbLength(filepath)
    # print(f"Number of Leaves: {num_leaves}\nTarget Leaf: {target_leaf}\nDistance Matrix:{distance_matrix.shape}\n {distance_matrix}")
    answer = LimbLength(target_leaf, distance_matrix)

    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(answer))

def ReadDistanceMatrix_LimbLength(filepath):
    with open(filepath) as f:
        num_leaves = int(f.readline().strip())
        target_leaf = int(f.readline().strip())
        distance_matrix = []
        for row in f:
            row = row.strip().split()
            row = [int(x) for x in row]
            distance_matrix.append(row)
        distance_matrix = np.array(distance_matrix)

    return num_leaves, target_leaf, distance_matrix

def LimbLength(target_leaf: int, distance_matrix: np.ndarray) -> int:
    """
    LimbLength: Finds the limb length of a given leaf in a simple tree.
    Input:
        target_leaf: The leaf that will have its limb length determined.
        distance_matrix: A distance matrix that contains the distance between each pair of leaves in the simple tree.
    Output:
        The limb length of target leaf.
    """
    min_dist = np.inf
    # for leaf1 in range(distance_matrix.shape[0]):
    #     if leaf1 != target_leaf:
    #         for leaf2 in range(leaf1+1, distance_matrix.shape[1]):
    #             if leaf2 != target_leaf:
    #                 dist1 = distance_matrix[leaf1, target_leaf]
    #                 dist2 = distance_matrix[leaf2, target_leaf]
    #                 dist3 = distance_matrix[leaf1, leaf2]
    #                 limb_length = int(0.5*(dist1+dist2-dist3))
    #                 if limb_length < min_dist:
    #                     min_dist = limb_length

    leaf1 = 0
    if target_leaf == 0:
        leaf1 = 1
    
    for leaf2 in range(distance_matrix.shape[0]):
        if leaf2 != target_leaf and leaf2 != leaf1:
            dist1 = distance_matrix[leaf1, target_leaf]
            dist2 = distance_matrix[leaf2, target_leaf]
            dist3 = distance_matrix[leaf1, leaf2]
            limb_length = int(0.5*(dist1+dist2-dist3))
            if limb_length < min_dist:
                min_dist = limb_length

    return min_dist

if __name__ == "__main__":
    main()