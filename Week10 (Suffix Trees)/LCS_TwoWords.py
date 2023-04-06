import os
from ConstructSuffixTree import *
from LCS import ConstructSubstring


def main():
    dirname = os.path.dirname(__name__)
    filepath = os.path.join(dirname, "inputs", "LCS_TwoWords", "dataset_876287_6.txt")
    with open(filepath) as f:
        word1 = f.readline().strip()
        word2 = f.readline().strip()
    answer = LCS_TwoWords(word1, word2)

    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, "w") as f:
        f.write(answer)


def LCS_TwoWords(w1: str, w2: str) -> str:
    """
    Finds the longest common substring between two words using a suffix tree.
    w1: the first word.
    w2: the second word.
    """
    w1 = w1 + "$"
    w2 = w2 + "#"
    word = "".join([w1, w2])
    tree = ConstructSuffixTree(word)
    tracker_dict = InitializeTrackerDictionary(tree, len(w1))
    path, _ = FindDeepestInternalNode_TwoWords(tree, [-1], 0, tracker_dict)
    substring = ConstructSubstring(tree, path, word)

    if substring == "":
        return "nan"
    return substring

def FindDeepestInternalNode_TwoWords(
    tree: dict[int, list[int]],
    curr_path: list[int],
    curr_depth: int,
    tracker_dict: dict[int, list[bool]] = None,
) -> tuple[list[int], int]:
    """
    Finds the deepest internal node in a two word suffix tree that has children branches that are from different words.
    Input:
        tree: The suffix tree represented as an adjacency dictionary of directed edges.
              (Leaves have length 0 in the adjacnecy dictionary).
        curr_path: The path to the current node.
        curr_depth: The depth of the current node.
        tracker_dict: A dictionary that keeps track of whether a node is part of a suffix from word 1, word 2 or both.
    Output:
        The path to the deepest internal node that has children from each word.
        The depth of the deepest internal node.
    """

    best_path, best_depth = curr_path, curr_depth    
    curr_node = curr_path[-1]
    for edge in tree[curr_node]:
        child_node = edge[0]
        if len(tree[child_node]) == 0: # Child node is a leaf.
            UpdateTrackerDict(tracker_dict, curr_node, child_node)
        else:
            new_path, new_depth = FindDeepestInternalNode_TwoWords(
                tree, [*curr_path, child_node], curr_depth + edge[2], tracker_dict
            )
            UpdateTrackerDict(tracker_dict, curr_node, child_node)

            if all(tracker_dict[child_node]) and new_depth > best_depth:
                best_path, best_depth = new_path, new_depth
            
    return best_path, best_depth

def UpdateTrackerDict(tracker_dict: dict[int, list[int]], curr_node:int, child_node: int):
    """
    Updates the boolean values of the current node in the tracker dictionary using the boolean values of the child node
    in the tracker dictionary.
    Input:
        tracker_dict: the tracker dictionary
        curr_node: the current node.
        child_node: the child node.
    """
    if tracker_dict[child_node][0]:
        tracker_dict[curr_node][0] = True
    if tracker_dict[child_node][1]:
        tracker_dict[curr_node][1] = True


def InitializeTrackerDictionary(tree: dict[int, list[int]], w2_start: int) -> dict[int, list[bool]]:
    """
    Intializes the tracker dictionary for whether a node is part of a suffix from word1, word2 or both.
    Input:
        tree: the suffix tree 
        w2_start: the starting index of word 2.
    Output:
        A dictionary that represents whether a node in the suffix tree is part of word 1 or word2. Since this is an
        initialization, all of the internal nodes will be false for both word 1 and word2 except for the root, 
        which will be true for both words. The leaves will only be true for which ever word generated that leaf.
    """
    tracker_dict = {}
    for key in tree:
        if len(tree[key]) == 0:
            if key < w2_start:
                tracker_dict[key] = [True, False]
            else:
                tracker_dict[key] = [False, True]
        elif key == -1:
            tracker_dict[key] = [True, True]
        else:
            tracker_dict[key] = [False, False]
    
    return tracker_dict

if __name__ == "__main__":
    main()