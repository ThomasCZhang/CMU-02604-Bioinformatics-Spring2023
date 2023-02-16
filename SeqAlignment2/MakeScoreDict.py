import os

dirpath = os.path.join(os.path.dirname(__file__), "blossum.txt")
pam250 = {}
with open(dirpath) as f:
    for idx, line in enumerate(f):
        if idx == 0:
            all_keys = line.strip().split()
            for val in all_keys:
                pam250[val] = {}
        else:
            line = line.strip().split()
            key = line[0]
            for idx1, val in enumerate(all_keys):
                pam250[key][val] = int(line[idx1+1])

newpath = os.path.join(os.path.dirname(__file__), "ScoreDict.txt")
with open(newpath, "w") as f:
    f.write("{\n")
    for idx, key in enumerate(pam250):
        if idx != 0:
            f.write("\n")
        f.write('"' + key + '"' + ":{")
        for idx1, key1 in enumerate(pam250[key]):
            f.write('"' + key1 + '": '+ f"{pam250[key][key1]: 2d}")
            if idx1 != len(pam250[key])-1:
                f.write(", ")
            if (idx1+1) % 5 == 0:
                f.write("\n\t ")
        f.write("},")
    f.write("\n}")