import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\inputs\\DPSeq")
    filepaths = glob(dirpath + "\\data*.txt")
    for filepath in filepaths:
        money, coins = ReadTestInputs_DPChange(filepath)
        answer = DPChange(money, coins)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(answer))

def ReadTestInputs_DPChange(FilePath: str) -> tuple[int, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        money: A quantity of money. Integer.

        coins: The available currency denominations.
    """
    with open(FilePath) as f:
        money = int(f.readline().strip())
        coins = [int(i) for i in f.readline().strip().split(" ")]

    return money, coins

def DPChange(money: int, coins: list[int])-> int:
    """
    DPChange uses dynamic programing (DP) to determine the minimum number of coins required to sum to a given
    money value. Assumes that all coins have values greater than or equal to 1.

    Input:
        money: The money (sum) we are trying to reach.

        coins: A list of currency values that can be used.
    Output:
        min_num: The minimum number of coins that have value equal to money.
    """
    MinNumCoins = [0 for i in range(money+1)]
    for m in range(1, money+1):
        MinNumCoins[m] = m
        for k in range(len(coins)):
            if m >= coins[k]:
                if MinNumCoins[m - coins[k]]+1 < MinNumCoins[m]:
                    MinNumCoins[m] = MinNumCoins[m- coins[k]]+1
    return MinNumCoins[money]

if __name__ == "__main__":
    main()