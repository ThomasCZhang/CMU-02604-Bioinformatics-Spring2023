import numpy as np

class NoNeighbors(Exception):
    pass

class Tree:
    def __init__(self, adjacency_dict: dict[str, list[str]]) -> None:
        self.vertices = {}
        self.root = None
        self.CreateTreeFromAdjacencyDict(adjacency_dict)
    
    def CreateTreeFromAdjacencyDict(self, adjacency_dict: dict[str, list[str]]):
        """
        Creates a list of Vertex objects from an adjacency list.
        """
        vertex_degree = {}
        for key in adjacency_dict:
            vertex_degree[key] = len(adjacency_dict[key])
        
        leaf_degree = 1

        while len(adjacency_dict) > 0:
            # Get vertices with degree = 1
            leaf_verticies = []
            for key in vertex_degree:
                if vertex_degree[key] == 1:
                    leaf_verticies.append(key)

            # If tree is rooted, the root will have degree 0 once all the children are processed.
            if len(leaf_verticies) == 0: 
                self.root = self.vertices[list(adjacency_dict.keys())[0]]
                break

            for current_name in leaf_verticies:
                if len(adjacency_dict[current_name]) == 0:
                    raise NoNeighbors("Current node has no neighbors in the adjacency list.")
                parent_name = adjacency_dict[current_name][0] # There should be only 1 neighbor if key is a leaf.

                if current_name not in self.vertices:
                    self.vertices[current_name] = Vertex(current_name)
                if parent_name not in self.vertices:
                    self.vertices[parent_name] = Vertex(parent_name)

                self.vertices[current_name].AddParent(self.vertices[parent_name])
                self.vertices[parent_name].AddChild(self.vertices[current_name])

                del adjacency_dict[current_name]
                adjacency_dict[parent_name].remove(current_name)
                vertex_degree[current_name] -= 1
                vertex_degree[parent_name] -= 1

    def GetLeaves(self)-> set[type['Vertex']]:
        """
        Finds and returns the leaf vertices of the tree.
        """
        leaves = set()
        for v in self.vertices.values():
            if len(v.children) == 0:
                leaves.add(v)
        return leaves
    
    def GetScore(self) -> int:
        """
        Returns the score of the tree. Score = sum of all edges in the tree.
        """
        score = 0
        for v in self.vertices.values():
            for _, weight in v.children:
                score += weight
        return score
        
class Vertex:
    def __init__(self, seq) -> None:
        self.sequence = seq
        self.children = []
        self.parent = []
        self.scores = []
        self.processed = False
    
    def AddChild(self, v: type["Vertex"], weight: float=0):
        """
        Add a child vertex to the current vertex.
        Input:
            v: The child vertex.
            weight: The weight of the child vertex.
        """
        self.children.append([v, weight])
    
    def AddParent(self, v: type["Vertex"], weight: float=0):
        """
        Add a child vertex to the current vertex.
        Input:
            v: The child vertex.
            weight: The weight of the child vertex.
        """
        self.parent.append([v, weight])

    def SetNewSequence(self, s: str):
        """
        Overwrites the vertex seqeunce with a new string:
        s: the new string.
        """
        self.sequence = s

    def ExtendNewSequence(self, s: str):
        """
        Extends the vertex sequence with a string.
        s: the extension string.
        """
        self.sequence += s

    def InitializeScoreMatrix(self, n: int, leaf: bool = False):
        """
        Initializes the score matrix for the vertex.
        """
        val = 0
        if leaf is True:
            val = np.inf
        self.scores = [{"A": val, "C": val, "G": val, "T": val} for _ in range(n)]

