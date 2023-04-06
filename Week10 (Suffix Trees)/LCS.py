import os
from ConstructSuffixTree import *

def main():
    dirname = os.path.dirname(__name__)
    filepath = os.path.join(dirname, "inputs", "LCS", "dataset_876287_5.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    answer = LongestCommonSubstring(word)
    
    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, "w") as f:
        f.write(answer)


def LongestCommonSubstring(word: str):
    """
    Finds the longest common substring in a word.
    Input:
        word: The word for which to find the longest common substring
    """
    word = word + "$"
    tree = ConstructSuffixTree(word)
    path, _ = FindDeepestInternalNode(tree, [-1], 0)
    substring = ConstructSubstring(tree, path, word)
    return substring

def ConstructSubstring(tree: dict[int, list[int]], path: list[int], word: str) -> str:
    """
    Reconstructs a substring from a suffix tree using a path of nodes.
    Input:
        tree: the suffix tree represented as a directional adjacency dictionary.
        path: a list of node numbers that should be used to travel through the adjacency dictionary.
        word: word used to generate the suffix tree
    Output:
        The reconstructed substring
    """
    substring = ""
    for i, node in enumerate(path[:-1]):
        for edge in tree[node]:
            if edge[0] == path[i+1]:
                substring = "".join([substring, word[edge[1]:edge[1]+edge[2]]])
                break
    return substring

def FindDeepestInternalNode(
    tree: dict[int, list[int]], current_path: list[int], current_depth: int
) -> tuple[list[int], int]:
    """
    Recursively finds the deepest internal node in a suffix tree.
    Input:
        tree: The suffix tree.
        current_path: the path to the current node.
        best_depth: The depth of the deepest internal node found so far.
    Output:
        The path taken to the deepest internal node.
    """
    current_node = current_path[-1]
    best_path, best_depth = current_path, current_depth

    for edge in tree[current_node]:
        child_node = edge[0]
        if len(tree[child_node]) > 0: # If child node is not a leaf
            new_path, new_depth = FindDeepestInternalNode(
                tree, [*current_path, child_node], current_depth + edge[2]
            )
            if new_depth > best_depth:
                best_depth = new_depth
                best_path = new_path

    return best_path, best_depth


if __name__ == "__main__":
    main()
