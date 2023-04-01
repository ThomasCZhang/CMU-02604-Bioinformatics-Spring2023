from ConstructTrie import *


def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "SuffixTree", "dataset_876287_4.txt")
    # filepath = os.path.join(dirname, "inputs", "SuffixTree", "input_4.txt")
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
    tree = {0: []}
    for s in suffixes:
        GrowTree(tree, s)
    return tree

def GrowTree(t: dict[int, list[list[int]]], s: str):
    """
    Grows a suffix tree using a string.
    Input:
        t: the suffix tree
        s: the string being used to grow the suffix tree
    """
    s_index = 0
    current_node = 0
    new_node = len(t)

    while True:
        mismatch_idx = 0
        for i, child in enumerate(t[current_node]):
            mismatch_idx = FindFirstMismatch(child[1], s[s_index:])
            if mismatch_idx != 0:
                break

        if mismatch_idx == 0: # Create branch from node
            t[current_node].append([new_node, s[s_index:]])
            t[new_node] = []
            break
        elif mismatch_idx == len(child[1]):
            current_node = child[0]
            s_index += mismatch_idx
        else: # create branch from inside an edge
            t[current_node].append([new_node, child[1][:mismatch_idx]]) # Add new internal node
            t[new_node] = [# Add two new children of internal node
                [child[0], child[1][mismatch_idx:]],
                [new_node+1, s[s_index+mismatch_idx:]]
            ]
            t[new_node+1] = [] 
            del t[current_node][i]
            break

def FindFirstMismatch(s1: str, s2: str) -> int:
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
