import os

def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "ConstructTrie", "dataset_876285_4.txt")
    words = ReadFile_ConstructTrie(filepath)
    trie = ConstructTrie(words)
    
    answerpath = os.path.join(dirname, "answer.txt")
    with open(answerpath, 'w') as f:
        for key in trie:
            for idx, edge in enumerate(trie[key]):
                if key != 0 or idx != 0:
                    f.write("\n")
                f.write(f"{key} {edge[0]} {edge[1]}")

def ReadFile_ConstructTrie(filepath) -> list[str]:
    with open(filepath) as f:
        line = f.readline().upper()
        return line.strip().split()

def ConstructTrie(words) -> dict[int, list[list[int, str]]]:
    """
    Constructs a trie from a list of words:
    Input:
        words: the list of words
    Output:
        A trie represented as a adjacency dictionary.
    """
    t = {0: []}
    for word in words:
        GrowTrie(word, t)
    return t

def GrowTrie(new_word: str, t: dict[int, list[list[int, str]]]):
    """
    Uses a new word to grow a trie.
    Input:
        new_word: The word being used to grow the trie.
        t: the trie represented as an adjacency dictionary.
    """
    current_node = 0 
    for i, letter in enumerate(new_word):
        edge = LookForEdge(letter, t[current_node])
        if edge is not None:
            current_node = edge[0]
        else:        
            first_unmatched_idx = i
            break
    
    new_node = len(t)
    for i, letter in enumerate(new_word[first_unmatched_idx:]):
        t[current_node].append([new_node, letter])
        t[new_node] = []
        current_node = new_node
        new_node += 1

def LookForEdge(edge_label: str, edge_list: list[list[int, str]]) -> list[int, str]:
    """
    Checks if an edge with a specific label is in a list of edges.
    Input:
        edge_label: The edge label being querried for.
        edge_list: The list of edges.
    """

    final_edge = None
    for edge in edge_list:
        if edge_label is edge[1]    :
            final_edge = edge
            break     
    return final_edge

if __name__ == "__main__":
    main()