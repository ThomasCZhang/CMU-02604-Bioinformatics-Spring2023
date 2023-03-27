import os
import numpy as np
from graph import *

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "SimpleParsimony", "dataset_876251_10.txt")
    _, adjacency_dict = ReadData_SmallParsimony(filepath)
    t = SmallParsimony(adjacency_dict)

    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(t.GetScore()))
        for vertex in t.vertices.values():
            for child, weight in vertex.children:
                f.write("\n"+ vertex.sequence + "->" + child.sequence + ":" + str(weight))
            for parent, weight in vertex.parent:
                f.write("\n" + vertex.sequence + "->" + parent.sequence + ":" + str(weight))   

def ReadData_SmallParsimony(filepath):
    with open(filepath) as f:
        num_leaves = f.readline().strip()
        adjacency_dict = {}
        for line in f:
            line = line.strip().split("->")
            key = line[0].upper()
            line[1] = line[1].upper() 
            if key in adjacency_dict:
                adjacency_dict[key].append(line[1])
            else:
                adjacency_dict[key] = [line[1]]
            if line[1] in adjacency_dict:
                adjacency_dict[line[1]].append(key)
            else:
                adjacency_dict[line[1]] = [key]
                
    return num_leaves, adjacency_dict

def SmallParsimony(adjacency_dict: dict[str, list[str]]) -> Tree:
    """
    Solves the small parsimony problem on a tree. Leaves of the tree have DNA sequences.
    Input:
        adjacency_dict: An evolutionary tree represented as an adjacency dictionary.
    """
    t = Tree(adjacency_dict)
    leaves = t.GetLeaves()
    ProcessLeaves(leaves)
    num_leaves = 0
    seq_length = len(list(leaves)[0].sequence)

    vertices_remaining = len(t.vertices) - num_leaves
    while vertices_remaining > 0:
        for v in t.vertices.values(): # Get "ripe" vertex"
            if CheckRipe(v):
                break
        FillScoreMatrix(v, seq_length)
        v.processed = True
        vertices_remaining -= 1

    UpdateTreeRootSequences(t)
    UpdateInternalTreeSequences(t, leaves)
    UpdateTreeEdgeWeights(t)
    return t
    
def UpdateTreeEdgeWeights(t: Tree):
    """
    Updates the edge weights of the tree.
    Input:
        t: Tree object
    """
    for v in t.vertices.values():
        for index, (child, _) in enumerate(v.children):
            weight = CountNumDifferentCharacters(v.sequence, child.sequence)
            v.children[index][1] = weight
            child.parent[0][1] = weight

def CountNumDifferentCharacters(s1: str, s2: str):
    """
    Count the number of different characters between two strings.
    """
    count = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            count += 1
    return count

def ProcessLeaves(leaves: set[Vertex]):
    """
    Initializes a set of leaf vertices with score matricies based on the sequences of the leaves.
    Also sets the leaves to a processed state.
    Input:
        leaves: a set of leaf vertices.
    """
    for leaf in leaves:
        n = len(leaf.sequence)
        FillScoreMatrix(leaf, n)
        leaf.processed = True

def UpdateTreeRootSequences(t: Tree):
    """
    Updates the sequence of the root vertex of a tree.
    Input:
        t: The tree.
    """
    root = t.root
    root.sequence = ""
    for score_dict in root.scores:
        min_score = np.inf
        keys = sorted(list(score_dict.keys()))
        for key in keys:
            if score_dict[key] < min_score:
                letter = key
                min_score = score_dict[key]
        root.sequence += letter

def UpdateInternalTreeSequences(t: Tree, leaves: set[Vertex]):
    """
    Updates the sequences corresponding to the tree's internal nodes.
    Input:
        t: The tree
        leaves: The leaves of the tree
    """
    stack = [t.root]
    while len(stack) > 0:
        current_vertex = stack.pop(-1)
        for index, (child, _) in enumerate(current_vertex.children):
            if child in leaves:
                break
            stack.append(child)
            child.sequence = ""
            for index_1, letter in enumerate(current_vertex.sequence):
                best_letter,_= CalculateBestScoreOfChild(letter, child.scores[index_1])
                child.sequence += best_letter


def FillScoreMatrix(v: Vertex, n: int):
    """
    Fills in the score matrix of the current vertex:
    Input:
        v: The current vertex.
        n: sequence length
    """
    if len(v.children) == 0: # Leaf node, no need for recurrence relationship.
        v.InitializeScoreMatrix(n, leaf = True)
        for index, l in enumerate(v.sequence):
            v.scores[index][l] = 0
    else: # Use recurrance relationship to fill in score matrix
        v.InitializeScoreMatrix(n, leaf = False)
        for child,_ in v.children:
            UpdateScoreUsingChild(v, child, n)            

def UpdateScoreUsingChild(v: Vertex, c: Vertex, n: int):
    """
    Calculates the parsimony scores contribution of a child vertex and updates the score matrix of the parent
    vertex using the calculated scores.
    Input:
        v: the parent vertex.
        c: the child vertex.
        n: the length of the sequence.
    """
    for index_0 in range(n): # Range over all the characters in the sequence.
        for key in v.scores[index_0]: # Range over the alphabet.
            _, score = CalculateBestScoreOfChild(key, c.scores[index_0])
            v.scores[index_0][key] += score

def CalculateBestScoreOfChild(letter: str, child_scores: dict) -> int:
    """
    Calculates the best parsimony score that a child vertex can contribute for a specific letter.
    Input:
        letter: the letter being compared against.
        c: score dictionary from the child vertex.
    """
    min_score = np.inf
    best_key = None
    for key in child_scores:
        score = child_scores[key]
        if key != letter:
            score += 1
        
        if score < min_score:
            min_score = score
            best_key = key
        
    return best_key, min_score

def CheckRipe(vertex: int) -> bool:
    """
    Checks if a vertex in a tree is "ripe" by checking if all its children verticies are processed.
    """
    if vertex.processed is True:
        return False
    for child,_ in vertex.children:
        if child.processed is False:
            return False
    return True

if __name__ == "__main__":
    main()
