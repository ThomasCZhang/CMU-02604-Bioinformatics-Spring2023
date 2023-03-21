class Graph:
    def __init__(self) -> None:
        self.nodes = dict()

    def AddNode(self, n):
        "n: Node object"
        if n in self.nodes:
            raise Exception("Not Allowed to have two nodes  with the same name.")
        self.nodes[n.name] = n

class Node:
    def __init__(self, name) -> None:
        self.name = name
        self.edges = []
    
    def AddEdge(self, edge):
        """
        edge: edge object
        """
        self.edge.append = edge

class Edge:
    def __init__(self, start, end, weight = 1) -> None:
        self.weight = weight
        self.start = start
        self.end = end
