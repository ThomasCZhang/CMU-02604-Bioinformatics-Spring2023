def main():
    import os
    from glob import glob
    
    dirpath = os.path.join(os.path.dirname(__file__), "Files")
    filepaths = glob(dirpath + "\\*.txt")
    for path in filepaths:
        with open(path) as f:
            for line in f:
                int1, int2 = line.split(sep =" ")
                print(f"From file {os.path.basename(f.name)}, the sum is: {SumTwoInts(int(int1), int(int2))}")          
                
# SumTwoInts takes two integers as an input and returns the sum.
def SumTwoInts(int1, int2):
    return int1 + int2

if __name__ == "__main__":
    main()