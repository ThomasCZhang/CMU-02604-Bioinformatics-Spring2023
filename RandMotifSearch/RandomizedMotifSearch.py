import os
import random
from ProfileMostProbableKmer import *
from glob import glob


def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\Inputs\RandomizedMotifSearch")
    filepaths = glob(dirpath + "\\*.txt")

    for path in filepaths:
        k, t, Dna = ReadTestFiles_RandomizedMotifSearch(path)
        BestMotifs = RandomizedMotifSearch(Dna, k, t)
        print(f"File {os.path.basename(path)}.")
        for motif in BestMotifs:
            print(motif)

def ReadTestFiles_RandomizedMotifSearch(filepath: str) -> tuple[int, int, list[dict[str, float]]]:
    """
    ReadData(): Reads data from a .txt file and formats that data for RandomizedMotifSearch.
    .txt Format should be as follows:\n
    Line 1: Two integers separated by a space. The first number is the  length of k-mers. The second is the number of
    DNA sequences.\n
    Line 2: The DNA sequences separated by spaces.
    Input:\n
    filepath (string): The path to the .txt file.\n
    Output:\n
    k (int): the number of letters in the motif.\n
    t (int): the number of sequences.\n
    Dna: (list of strings): the scoring profile for the motifs.
    """
    with open(filepath) as f:
        Dna = []
        for ind, line in enumerate(f):
            if ind == 0:
                string = line.strip()
                k, t = [int(x) for x in string.split(" ")]
            else:
                string = line.strip()
                Dna.extend(string.split(" "))
    return (k, t, Dna)

# Insert your RandomizedMotifSearch function here, along with any subroutines you need
def RandomizedMotifSearch(Dna: list[str], k: int, t: int) -> list[str]:
    """
    RandomizedMotifSearch(Dna, k, t): Finds the optimal k-length motifs from t Dna sequences by running
    RandomizedMotifSearch many times.\n
    Input:\n
    Iterations (int): The number of times to run RandomizedMotifSearch
    Dna (list of str): DNA sequences to find motifs from.\n
    k (int): the length of the motifs.\n
    t (int): the number of DNA sequences.\n
    Output:\n
    BestMotifs (list of str): The best motifs. 
    """
    Iterations = 1000
    BestMotifs = ["" for i in range(t)]
    BestScore = k*t
    for i in range(Iterations):
        Motifs, NewScore = SingleRandomizedMotifSearch(Dna, k, t)
        if NewScore < BestScore:
            BestScore = NewScore
            BestMotifs = Motifs
    return BestMotifs


def SingleRandomizedMotifSearch(Dna: list[str], k: int, t: int) -> tuple[list[str], int]:
    """
    SingleRandomizedMotifSearch(Dna, k, t): Finds the optimal k-length motifs from t Dna sequences.\n
    Input:\n
    Dna (list of str): DNA sequences to find motifs from.\n
    k (int): the length of the motifs.\n
    t (int): the number of DNA sequences.\n
    Output:\n
    BestMotifs (list of str): The best motifs. \n
    BestScore (float): The best score.
    """
    BestMotifs = ChooseRandomMotifs(Dna, k)
    while True:
        Profile = CreateProfile(BestMotifs)
        Motifs = FindMostProbableMotifs(Dna, k, Profile)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
        else:
            BestScore = Score(BestMotifs)
            return BestMotifs, BestScore


def ChooseRandomMotifs(Dna: list[str], k: int) -> list[str]:
    """
    ChooseRandommotifs(Dna): Chooses random motifs from a list of DNA sequences.\n
    Input:\n
    Dna (list of strings): List of DNA sequences.\n
    k (int): The length of the motifs.\n
    Output:\n
    Motifs (list of strings): List of motifs
    """
    Motifs = ["" for i in range(len(Dna))]
    for i, Sequence in enumerate(Dna):
        Start = random.randrange(len(Sequence)-k+1)
        Stop = Start + k
        Motifs[i] = Sequence[Start:Stop]
    return Motifs


def CreateProfile(Motifs: list[str]) -> list[dict[str, int]]:
    """
    CreateProfile: Creates a profile matrix based on a list of motifs.\n
    Input:\n
    Motifs (list of strings): The list of motifs used to generate the profile matrix.\n
    k (int): length of the motifs.\n
    Output:\n
    Profile (list of dictionaries(string to int)): The profile matrix represented as a list of dictionaries.
    """
    Pseudocount = 1
    Profile = CountLettersAtEachPosition(Motifs)
    for Dictionary in Profile:
        AddConstantToDictionaryValues(Pseudocount, Dictionary)

    NormalizationFactor = len(Motifs) + 4*Pseudocount
    for Dictionary in Profile:
        for Key in Dictionary:
            Dictionary[Key] /= NormalizationFactor

    return Profile


def InitializeProfileWithPseudoCounts(PseudocountValue: float, k: int) -> list[dict[str, int]]:
    """
    InitializeProfileWithPseduoCounts(): Creates a profile matrix (represented as a list of dictionaries of string
    to int) with pseudo counts.\n
    Input:\n
    PsuedocountValue (float): value to use as pseudocount for profile matrix.\n
    k (int): the length of the k-mers.
    Output:\n
    Profile (list of dictionaries(string to int)): The initialized profile matrix.
    """
    Profile = [{} for i in range(k)]
    for i in range(k):
        Profile[i]["A"] = PseudocountValue
        Profile[i]["C"] = PseudocountValue
        Profile[i]["G"] = PseudocountValue
        Profile[i]["T"] = PseudocountValue
    return Profile


def FindMostProbableMotifs(Dna: list[str], k: int, Profile: list[dict[str, int]]) -> list[str]:
    """
    FindMostProbableMotifs(Dna, k, Profile): Find the most probable k-length motifs from a list of Dna sequences using
    Profile.\n
    Input:\n
    Dna (list of strings): The list of Dna sequences.\n
    k (int): The length of the motifs. \n
    Profile (list of dictionaries(string to int)): The Profile matrix used to find the most likely motifs.\n
    Output:\n
    BestMotifs (list of strings): The best motifs based on the profile matrix.
    """
    BestMotifs = ["" for i in range(len(Dna))]
    for i, Text in enumerate(Dna):
        BestMotifs[i] = ProfileMostProbableKmer(Text, k, Profile)
    return BestMotifs


def Score(Motifs: list[str]) -> float:
    """
    Score(Motifs, Profile). Scores a list of motifs according to a profile.\n
    Input:\n
    Motifs (list of strings): Motifs to be scored.\n
    Profile (list of dictionaries(string to int). Profile used to score the motifs.\n
    Output:\n
    FinalScore (float): The score of the motifs.
    """
    FinalScore = 0
    CountProfile = CountLettersAtEachPosition(Motifs)
    for Dictionary in CountProfile:
        MaxKey = FindMaxKey(Dictionary)
        for Key in Dictionary:
            if Key != MaxKey:
                FinalScore += Dictionary[Key]
    return FinalScore



def CountLettersAtEachPosition(Motifs: list[str]) -> list[dict[str, int]]:
    """
    Counts the number of times a letter appears in each position of the motifs.\n
    Input: \n
    Motifs (list of strings): The motifs being analyzed.\n
    Output: \n
    Counts (list[dict[str, int]]): The number of times a letter appears in a position for a collection of motifs.
    """
    k = len(Motifs[0])
    Counts = [{"A": 0, "C": 0, "G": 0, "T": 0,} for i in range(k)]
    for Text in Motifs:
        for j, Letter in enumerate(Text):
            Counts[j][Letter.upper()] += 1
    return Counts

def FindMaxKey(Dictionary: dict[str, int]) -> str:
    """
    Returns the key containing the max value in a dictionary of strings to ints.\n
    Input:\n
    Dictionary (dictionary of string to int): The dictionary being analyzed.\n
    Output:\n
    MaxKey (str): The key of the max value in Dictionary.
    """
    MaxKey = list(Dictionary.keys())[0]
    for Key in Dictionary:
        if Dictionary[Key] > Dictionary[MaxKey]:
            MaxKey = Key
    return MaxKey

def AddConstantToDictionaryValues(Constant: float, Dictionary: dict[str, int]):
    """
    Adds some constant value to all values in a Dictionary of strings to ints.\n
    Input:\n
    Constant (float): The value being added to all elements in the dictionary.\n
    Dictionary (dictionary of strings to ints): The dictionary being modified.
    """
    for Key in Dictionary:
        Dictionary[Key] += Constant

if __name__ == "__main__":
    main()


# def Score(Motifs: list[str], Profile: list[dict[str,int]]) -> float:
#     """
#     Score(Motifs, Profile). Scores a list of motifs according to a profile.\n
#     Input:\n
#     Motifs (list of strings): Motifs to be scored.\n
#     Profile (list of dictionaries(string to int). Profile used to score the motifs.\n
#     Output:\n
#     FinalScore (float): The score of the motifs.
#     """
#     for Index, Motif in enumerate(Motifs):
#         if Index == 0:
#             FinalScore = ScoreKMer(Motif, Profile)
#         else:
#             FinalScore += ScoreKMer(Motif, Profile)
#     return FinalScore
