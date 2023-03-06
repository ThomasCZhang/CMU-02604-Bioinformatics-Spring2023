import os, re
import numpy as np
from typing import Iterable

def main():
    dirpath = os.path.join(os.path.dirname(__file__), "inputs\\TwoBreakDistance")
    filepath = os.path.join(dirpath, "dataset_876232_4.txt")
    genome1, genome2 = ReadData_TwoBreakDistance(filepath)
    answer = TwoBreakDist(genome1, genome2)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, "w") as f:
        f.write(str(answer))


def ReadData_TwoBreakDistance(filepath):
    genomes = []
    with open(filepath) as f:
        for line in f:
            genome = []
            line = re.split("[\(\)\n]", line)
            line = list(filter(None, line))
            for chromasome in line:
                chromasome = chromasome.split()
                sequence = tuple([int(x) for x in chromasome])
                genome.append(sequence)
            genomes.append(genome)
    return tuple(genomes)

def TwoBreakDist(p: list[tuple[int]],q: list[tuple[int]]) -> int:
    """
    TwoBreakDist: Determines the two break distance between two genome synteny block sequences.
    Input:
        p, q: The two input genomes.
    Output:
        The two break distance.
    """
    adj_dict = GenerateAdjacencyDict([p ,q])
    num_blocks = int(len(adj_dict)/2)
    keys = list(adj_dict.keys())
    num_cycles = 0
    curr_vertex = keys[0]

    while len(adj_dict) > 0:
        if curr_vertex not in adj_dict:
            curr_vertex = keys[0]

        next_vertex = adj_dict[curr_vertex][0]
        
        del adj_dict[curr_vertex]
        keys.remove(curr_vertex)
        
        if next_vertex in adj_dict:
            adj_dict[next_vertex].remove(curr_vertex)
        else:
            num_cycles += 1
        curr_vertex = next_vertex
        
    return num_blocks - num_cycles

def GenerateAdjacencyDict(genomes: Iterable[list[tuple[int]]]) -> dict[int, list[int]]:
    """
    Generates the adjacency lists for a list of sequences of synteny blocks.
    Input:
        sequences: An iterable containing the sequences.
    Output:
        Adjacency list as a dictionary.
    """
    adj_dict = dict()
    for genome in genomes:
        new_adj_dict = GenomeAdjacencyList(genome)
        CombineDict(adj_dict, new_adj_dict)
    return adj_dict

def GenomeAdjacencyList(p: list[tuple[int]]) -> dict[int, list[int]]:
    """
    GenomeAdjacnecyList: Generates an adjacency list for a genomes. Each genome is represented as a sequence
    of circular chromasomes. Chromsomes are a sequence synteny blocks.
    Input:
        p: The genome represented as a sequence of chromasomes. Assumes that synteny blocks are numbered from 1 to 
        the number of synteny blocks. Negative numbers mean a synteny block is reversed.
    Output:
        An adjacency list represented as a dictionary. Key = one vertex of an edge, value = list of vertices that 
        share an edge with the key vertex.
    """ 
    adj_dict = dict()
    for chromasome in p:
        new_dict = ChromasomeAdjacencyList(chromasome)
        CombineDict(adj_dict, new_dict)
    return adj_dict

def ChromasomeAdjacencyList(p: tuple[int]) -> dict[int, list[int]]:
    """
    Chromasome: Generates an adjacency list for circular chromasomes. Each genome is represented as a sequence
    of synteny blocks.
    Input:
        p: The genome represented as a sequence of synteny blocks. Assumes that blocks are numbered from 1 to the number
        of synteny blocks. Negative numbers mean a synteny block is reversed.
    Output:
        An adjacency list represented as a dictionary. Key = one vertex of an edge, value = list of vertices that 
        share an edge with the key vertex.
    """ 
    adj_dict = {}
    for block in p:
        adj_dict[abs(block)*2-1] = list()
        adj_dict[abs(block)*2] = list()

    for idx in range(len(p)):
        curr_block = p[idx]
        prev_block = p[idx-1]

        if curr_block > 0:
            curr_vertex = abs(curr_block)*2-1
        else:
            curr_vertex = abs(curr_block)*2
        if prev_block > 0:
            prev_vertex = abs(prev_block)*2
        else:
            prev_vertex = abs(prev_block)*2-1

        adj_dict[curr_vertex].append(prev_vertex)
        adj_dict[prev_vertex].append(curr_vertex)
    
    return adj_dict

def CombineDict(dict1: dict[int, list[int]], dict2: dict[int, list[int]]):
    """
    Combines two dictionaries with integer keys and lists as values.
    If dict1 shares a key with dict2, then the value of the shared key in dict2 is extended to that key in dict1.
    Input:
        dict1: dictionary 1. This will be modified.
        dict2: dictionary 2. The contents of this dictionary will be added to dict1.
    """
    for key, value in dict2.items():
        if key in dict1:
            dict1[key].extend(value)
        else:
            dict1[key] = value


if __name__ == "__main__":
    main()