import os
from SmallParsimony import *

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "UnrootedSimpleParsimony", "dataset_876251_12.txt")
    _, adjacency_dict = ReadData_UnrootedSimpleParsimony(filepath)
    t = UnrootedSmallParsimony(adjacency_dict)

    answerpath = os.path.join(dirpath, "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(t.GetScore()))
        for vertex in t.vertices.values():
            for child, weight in vertex.children:
                f.write("\n"+ vertex.sequence + "->" + child.sequence + ":" + str(weight))
            for parent, weight in vertex.parent:
                f.write("\n" + vertex.sequence + "->" + parent.sequence + ":" + str(weight))   

def ReadData_UnrootedSimpleParsimony(filepath):
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
                
    return num_leaves, adjacency_dict

def UnrootedSmallParsimony(adjacency_dict: dict[str, list[str]]) -> Tree:
    """
    Solve the small parsimony problem for an unrooted tree.
    Input:
        adjacency_dict: A unrooted binary tree.
    """
    AddRootToAdjacencyDict(adjacency_dict)
    t = SmallParsimony(adjacency_dict)
    RemoveRoot(t)
    return t

def RemoveRoot(t: Tree):
    """
    Adds an edge connecting the children of the root of a binary tree. Then removes the root to create a unrooted tree.
    Input:
        t: The rooted Tree object.
    """
    children = []
    for child,_ in t.root.children:
        children.append(child)
    child1, child2 = children
    weight = CountNumDifferentCharacters(child1.sequence, child2.sequence)
    child1.parent = [[child2, weight]]
    child2.children.append([child1, weight])
    child2.parent = []
        
    for key, item in t.vertices.items():
        if item is t.root:
            break
    del t.vertices[key]
    t.root = None
    
def AddRootToAdjacencyDict(adjacency_dict: dict[str, list[str]]):
    """
    Adds a new node to a unrooted binary tree to create a rooted binary tree.
    Input:
        adjacnecy_dict: The unrooted binary tree.
    """
    v1, v2 = FindTwoInternalNodes(adjacency_dict)
    new_vertex = str(len(adjacency_dict))
    adjacency_dict[new_vertex] = [v1, v2]
    adjacency_dict[v1].remove(v2)
    adjacency_dict[v2].remove(v1)
    adjacency_dict[v1].append(new_vertex)
    adjacency_dict[v2].append(new_vertex)

def FindTwoInternalNodes(adjacency_dict: dict[str, list[str]]):
    """
    Finds two internal nodes in an unrooted tree.
    Input:
        adjacency_dict: An unrooted tree represented as an adjacency dictionary.
    """
    for key, neighbors in adjacency_dict.items():
        if len(neighbors) >= 3:
            for neighbor in neighbors:
                if len(adjacency_dict[neighbor]) >= 3:
                    return key, neighbor


if __name__ == "__main__":
    main()