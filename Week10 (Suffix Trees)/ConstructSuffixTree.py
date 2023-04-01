from ConstructTrie import *


def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "inputs", "SuffixTree", "dataset_876287_4.txt")
    with open(filepath) as f:
        word = f.readline().strip()
    word_list = ConstructSuffixTree(word)

    # answerpath = os.path.join(dirname, "answer.txt")
    # with open(answerpath, 'w') as f:
    #     for idx, word in enumerate(word_list):
    #         if idx != 0:
    #             f.write(' ')
    #         f.write(word)


def ConstructSuffixTree(word: str) -> dict[int, list[list[int]]]:
    suffixes = [word[i:] for i in range(len(word))]
    trie = ConstructTrie(suffixes)
    word_list = TraverseTrie(trie, 0, "")
    # return sorted(word_list)


def TraverseTrie(
    t: dict[int, list[list[int]]], current_node: int, current_string: str = None
) -> list[str]:
    """
    Traverse down a suffix trie and get a list of "words" that represent non-branching chains of nodes.
    Input:
        t: the suffix trie
        current_node: the current node
        current_string: the current growing word.
    """
    word_list = []
    if current_string is None:
        current_string = ""
    
    if len(t[current_node]) == 1:
        edge = t[current_node][0]
        current_string += edge[1]
        word_list.extend(TraverseTrie(t, edge[0], current_string))
    else:
        if len(current_string) > 0:
            word_list.append(current_string)
        for edge in t[current_node]:
            word_list.extend(TraverseTrie(t, edge[0], edge[1]))

    return word_list


if __name__ == "__main__":
    main()
