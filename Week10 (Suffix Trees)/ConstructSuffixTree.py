from ConstructTrie import *


def main():
    dirname = os.path.dirname(__file__)
    # filepath = os.path.join(dirname, "inputs", "SuffixTree", "dataset_876287_4.txt")
    filepath = os.path.join(dirname, "inputs", "SuffixTree", "input_0.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    suffix_tree = ConstructSuffixTree(word)
    print(suffix_tree[-1])
    for k in suffix_tree:
        for edge in suffix_tree[k]:
            print(f"{k} -> node: {edge[0]} suffix: {word[edge[1]:edge[1]+edge[2]]}")
    word_list = GetAllEdgeLabels(suffix_tree, word)

    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, "w") as f:
        for i, w in enumerate(word_list):
            if i != 0:
                f.write(" ")
            f.write(w)


def ConstructSuffixTree(word: str) -> dict[int, list[list[int]]]:
    """
    Constructs a suffix tree by first creating the suffix trie then compressing the edges down.
    Input:
        word: The word being used to construct a suffix tree.
    """
    tree = {-1: []}
    for i in range(len(word)):
        tree[i] = []
    for i in range(len(word)):
        GrowTree(tree, word, i)
    return tree


def GetAllEdgeLabels(t: dict[int, list[list[int]]], word: str) -> list[str]:
    """
    Converts all the start position, weight labels of a suffix tree into the actual string they represent.
    Input:
        t: The suffix tree as an adjacency dictionary (directional).
        word: the word used to build the suffix tree.
    """
    word_list = []
    for key in t:
        for edge in t[key]:
            word_list.append(word[edge[1] : edge[1] + edge[2]])
    return word_list


def GrowTree(t: dict[int, list[list[int]]], word: str, start_idx: int):
    """
    Grows a suffix tree using a suffix of a word.
    Input:
        t: the suffix tree
        word: the word being used to grow the suffix tree
        start_idx: the index corresponding to start position of the word
    """
    new_node = len(t) - 1

    suffix_length = len(word) - start_idx
    first_mismatch_of_suffix, current_node, edge_index, num_matches_on_edge = FindBranchPoint(t, word, start_idx)
    if num_matches_on_edge == 0:  # Create branch from node
        t[current_node].append(
            [start_idx, start_idx + first_mismatch_of_suffix, suffix_length - first_mismatch_of_suffix]
        )
    else:  # split an edge into two paths
        edge = t[current_node][edge_index]
        t[new_node] = [  # Add two new children of internal node
            [edge[0], edge[1] + num_matches_on_edge, edge[2] - num_matches_on_edge],
            [start_idx, start_idx + first_mismatch_of_suffix, suffix_length - first_mismatch_of_suffix],
        ]
        t[current_node][edge_index] = [
            new_node,
            edge[1],
            num_matches_on_edge,
        ]  # replace edge with edge to new internal node


def FindBranchPoint(t: dict[int, list[list[int]]], word: str, start_index: int) -> tuple[int, int, int | None, int]:
    """
    Finds the branching point in a suffix tree for adding a new suffix.
    Input:
        t: The suffix tree
        word: the word the suffix tree is being constructed from
        start_index: the starting position in word of the current suffix being added.
    Output:
        The index of the first mismatched character on the suffix
        the deepest node in the suffix tree that had all matches with the suffix
        the index of the edge that the mismatch occured along
        the distance along the edge that the suffix matched
    """
    first_mismatch_of_suffix = 0  # Index of the first mismatched character in the suffix being added.
    current_node = -1

    while True:
        num_matches_on_edge = 0
        edge_index = None
        for edge_index, edge in enumerate(t[current_node]):
            edge_label = word[edge[1] : edge[1] + edge[2]]
            num_matches_on_edge = CountMatches(edge_label, word[start_index + first_mismatch_of_suffix :])
            if num_matches_on_edge != 0:
                break
        first_mismatch_of_suffix += num_matches_on_edge
        if len(t[current_node]) == 0 or num_matches_on_edge != edge[2]:  # We found the splitting point
            return first_mismatch_of_suffix, current_node, edge_index, num_matches_on_edge
        else:  # We reached a node that branches
            current_node = edge[0]


def CountMatches(s1: str, s2: str) -> int:
    """
    Returns the index of the first mismatch between two strings.
    s1, s2: the two strings being compared.
    """
    i = 0
    max_length = min(len(s1), len(s2))
    while i < max_length and s1[i] == s2[i]:
        i += 1
    return i


if __name__ == "__main__":
    main()
