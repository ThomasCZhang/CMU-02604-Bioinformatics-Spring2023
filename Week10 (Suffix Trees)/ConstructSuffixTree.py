from ConstructTrie import *


def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "SuffixTree", "dataset_876287_4.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    suffix_tree = ConstructSuffixTree(word)

    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, 'w') as f:
        for idx, key in enumerate(suffix_tree):
            for idx2, edge in enumerate(suffix_tree[key]):
                if idx != 0 or idx2 != 0:
                    f.write(' ')
                f.write(edge[1])


def ConstructSuffixTree(word: str) -> dict[int, list[list[int]]]:
    """
    Constructs a suffix tree by first creating the suffix trie then compressing the edges down.
    Input:
        word: The word being used to construct a suffix tree.
    """
    suffixes = [word[i:] for i in range(len(word))]
    trie = ConstructTrie(suffixes)
    CompressTrie(trie, 0)
    return trie


def CompressTrie(
    t: dict[int, list[list[int]]], current_node: int, current_edge: list[int | str] = None
):
    """
    Traverse down a suffix trie and get a list of "words" that represent non-branching chains of nodes.
    Input:
        t: the suffix trie
        current_node: the current node
        current_string: the current growing word.
    """
    if len(t[current_node]) == 1:
        edge = t[current_node][0]
        current_edge[0] = edge[0]
        current_edge[1] += edge[1]
        CompressTrie(t, edge[0], current_edge)
        del t[current_node]
    else:
        for edge in t[current_node]:
            CompressTrie(t, edge[0], edge)



if __name__ == "__main__":
    main()
