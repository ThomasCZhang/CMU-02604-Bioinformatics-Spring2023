import os
from SmallParsimony import *

def main():
    dirpath = os.path.dirname(__file__)
    filepath = os.path.join(dirpath, "inputs", "UnrootedSimpleParsimony", "input_0.txt")
    _, adjacency_dict = ReadData_UnrootedSimpleParsimony(filepath)
    t = Tree(adjacency_dict)
    for v in t.vertices.values():
        for child,_ in v.children:
            print(f"{v.sequence}->{child.sequence}")
        print(f"{v.sequence}->{v.parent[0][0].sequence}")
    # SmallParsimony(t)

    # answerpath = os.path.join(dirpath, "answer.txt")
    # with open(answerpath, 'w') as f:
    #     f.write(str(t.GetScore()))
    #     for vertex in t.vertices.values():
    #         for child, weight in vertex.children:
    #             f.write("\n"+ vertex.sequence + "->" + child.sequence + ":" + str(weight))
    #         for parent, weight in vertex.parent:
    #             f.write("\n" + vertex.sequence + "->" + parent.sequence + ":" + str(weight))   

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
            # if line[1] in adjacency_dict:
            #     adjacency_dict[line[1]].append(key)
            # else:
            #     adjacency_dict[line[1]] = [key]
                
    return num_leaves, adjacency_dict

if __name__ == "__main__":
    main()