from TwoBreakDistance import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "inputs\\TwoBreakSorting")
    filepath = os.path.join(dirpath, "dataset_876232_5.txt")
    genome1, genome2 = ReadData_TwoGenomes(filepath)
    answer = TwoBreakSort(genome1, genome2)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, "w") as f:
        for idx0, sequence in enumerate(answer):
            if idx0 != 0:
                f.write('\n')
            for sub_sequence in sequence:
                f.write("(")
                for idx1, block_num in enumerate(sub_sequence):
                    if idx1 != 0:
                        f.write(" ")

                    if block_num > 0:
                        block_num = "+" + str(block_num)
                    else:
                        block_num = str(block_num)

                    f.write(block_num)
                f.write(")")

def TwoBreakSort(p: list[list[int]], q: list[list[int]]) -> list[list[list[int]]]:
    """
    TwoBreakSort: Converts the breakpoint graph of one sequence into the breakpoint graph of another sequence. Thereby
    "sorting" the first sequence into the second.
    Input:
        p: The sequence that will be converted
        q: The target sequence
    Output:
        A list containing all the intermediate states from p being sorted into q.
    """
    sequence_list = [p]
    reverse = q[0][0] < 0
    # Generate the breakpoint graphs as an adjacency dictionaries.
    adj_dict1 = GenerateAdjacencyDict([p])
    adj_dict2 = GenerateAdjacencyDict([q])
    keys = set(adj_dict2.keys())
    
    # Set object to keep track of breakpoint graph verticies that are now part of a trivial cycle.
    while len(keys) > 0:
        vertex1 = list(keys)[0]
        vertex2 = adj_dict2[vertex1][0]
        if adj_dict1[vertex1][0] != vertex2:
            edge1 = [vertex1, adj_dict1[vertex1][0]]
            edge2 = [vertex2, adj_dict1[vertex2][0]]
            TwoBreak(adj_dict1, edge1, edge2)
            sequence_list.append(BreakpointToSequence(adj_dict1, list(adj_dict2.keys()), reverse))
        keys.remove(vertex1)
        keys.remove(vertex2)
    return sequence_list

def TwoBreak(dict1: dict[int, list[int]], edge1: list[int], edge2: list[int]):
    """
    TwoBreak: Takes a breakpoint graph and performs a two-break edge rearrangement. If (u1, u2) and (v1, v2) are the 
    two edges being broken in the original breakpoint graph. Then the resulting graph after TwoBreak() will contain the
    edges (u1, v1) and (u2, v2).
    Input:
        dict1: The breakpoint graph represented as an adjacency dictionary.
        edge1: The first edge in the breakpoint graph being broken. (u1, u2)
        edge2: The second edge in the breakpoint graph being broken. (v1, v2)
    Output:
        The new breakpoint graph with modified edges.
    """
    dict1[edge1[0]] = [edge2[0]]
    dict1[edge2[0]] = [edge1[0]]
    
    dict1[edge1[1]] = [edge2[1]]
    dict1[edge2[1]] = [edge1[1]]

def BreakpointToSequence(adj_dict: dict[int, list[int]], keys: list[int] = None, reverse: bool = False) -> list[list[int]]:
    """
    BreakpointToSequence: Takes a breakpoint graph, represented as an adjacency dictionary and returns the sequence
    that the breakpoint graph represents.
    Input:
        adj_dict: The break point graph, represented as an adjacency dict.
        keys: The keys in the adjacency dict.
    Output:
        The sequence that represents the breakpoint graph.
    """
    result = list()
    subsequence = list()
    result.append(subsequence)
    if keys == None:
        keys = list(adj_dict.keys())
    
    vertex1 = keys[0]
    if reverse:
        vertex1 = keys[1]

    while len(keys) > 0:
        if vertex1 not in keys:
            subsequence = list()
            result.append(subsequence)
            vertex1 = list(keys)[0]
        
        if vertex1%2 == 0:
            subsequence.append(-vertex1//2)
            vertex2 = vertex1-1        
        else:
            subsequence.append((vertex1+1)//2)
            vertex2 = vertex1+1      

        keys.remove(vertex1)
        keys.remove(vertex2)
        vertex1 = adj_dict[vertex2][0]

    return result

if __name__ == "__main__":
    main()