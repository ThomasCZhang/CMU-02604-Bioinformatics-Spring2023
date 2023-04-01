from ConstructTrie import *


def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "SuffixTree", "dataset_876287_4.txt")
    # filepath = os.path.join(dirname, "inputs", "SuffixTree", "input_2.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    suffix_tree = ConstructSuffixTree(word)
    word_list = GetAllEdgeLabels(suffix_tree, word)

    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, 'w') as f:
        for i, w in enumerate(word_list):
            if i != 0:
                f.write(' ')
            f.write(w)
            
def ConstructSuffixTree(word: str) -> dict[int, list[list[int]]]:
    """
    Constructs a suffix tree by first creating the suffix trie then compressing the edges down.
    Input:
        word: The word being used to construct a suffix tree.
    """
    tree = {-1:[]}
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
            word_list.append(word[edge[1]:edge[1]+edge[2]])
    return word_list

def GrowTree(t: dict[int, list[list[int]]], word: str, start_idx: int):
    """
    Grows a suffix tree using a suffix of a word.
    Input:
        t: the suffix tree
        word: the word being used to grow the suffix tree
        start_idx: the index corresponding to start position of the word
    """
    s = word[start_idx:]
    s_index = 0 # Index along s
    current_node = -1
    new_node = len(t)-1

    while True:
        num_matches = 0
        for i, edge in enumerate(t[current_node]):
            edge_label = word[edge[1]:edge[1]+edge[2]]
            num_matches = CountMatches(edge_label, s[s_index:])
            if num_matches != 0:
                break

        if num_matches == 0: # Create branch from node
            t[current_node].append([start_idx, start_idx+s_index, len(s)-s_index])
            # t[current_node].append([new_node, s[s_index:]])
            t[new_node] = []
            break
        elif num_matches == edge[2]: # We reached a node that branches
            current_node = edge[0]
            s_index += num_matches
        else: # create branch from inside an edge
            t[current_node][i] =[new_node, edge[1], num_matches] # replace edge with edge to new internal node
            t[new_node] = [# Add two new children of internal node
                [edge[0], edge[1]+num_matches, edge[2]-num_matches],
                [start_idx, start_idx+s_index+num_matches, len(s)-s_index-num_matches]
            ]
            break

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
