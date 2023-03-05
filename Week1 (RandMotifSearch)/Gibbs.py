import os
from glob import glob
from RandomizedMotifSearch import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\Inputs\Gibbs")
    filepaths = glob(dirpath + "\\data*.txt")
    for path in filepaths:
        k, t, N, Dna = ReadTestFilesGibbsSampler(path)
        BestMotifs = GibbsSampler(Dna, k, t, N)

    AnswerPath = os.path.join(os.path.dirname(__file__), "GibbsAnswer.txt")
    with open(AnswerPath, "w") as f:
        for Line in BestMotifs:
            f.write(Line)
            f.write("\n")

def ReadTestFilesGibbsSampler(FilePath: str) -> tuple[int, int, int, list[str]]:
    """
    ReadTestFilesGibbsSampler: Reads test files for Gibbs Sampler.
    Line 1: k, t, N
    Line 2 to the end: The Dna strands.

    Input:
        FilePath (str): Path to the test file.

    Output:
        k (int): the number of characters in the KMer
        
        t (int): the number of dna strings
        
        N (int): the number of time 
    """
    with open(FilePath) as f:
        Dna = []
        for index, line in enumerate(f):
            line = line.strip()
            if index == 0:
                k, t, N = [int(i) for i in line.split(" ")]
            else:
                Dna.extend(line.split(" "))
    return k, t, N, Dna

def GibbsSampler(Dna: list[str], k: int, t: int, N: int) -> list[str]:
    """
    Performs the Gibbs sampling to find the best motifs from a collection of DNA sequences.
    Input:
        Dna (list[str]): Collection of DNA sequences.
        
        k (int): number of characters in a K_mer.
        
        t (int): number of DNA.
        
        N (int): number of times a motif is replaced in one run of GibbsSampler from a random set of start motifs.

    Output: 
        BestMotif (list[str]): the best motifs.
    """
    NumIterations = 20
    BestMotifs = ["" for i in range(t)]
    BestScore = k*t # Arbitrary large starting score 
    for i in range(NumIterations):
        Motifs = GibbsSamplerOneIteration(Dna, k, t, N)
        NewScore = Score(Motifs)
        if NewScore < BestScore:
            BestMotifs = Motifs
            BestScore = NewScore
    return BestMotifs


def GibbsSamplerOneIteration(Dna: list[str], k: int, t: int, N: int) -> list[str]:
    """
    Finds the best motif reachable from a single random starting collection of motifs from Dna using Gibbs Sampling.

    Input:

        Dna (list[str]): List of strings. Each string represents a separate Dna strand.
    
        k (int): the number of bases in each motif.
    
        t (int): the number of Dna strands.
    
        N (int): the number of times to run the Gibbs iteration.

    Output:
        BestMotifs (list[str]): A list of strings, each string contains the best motif from each Dna strand as determined
        by Gibbs sampling.
    """
    BestMotifs = ChooseRandomMotifs(Dna, k)
    for i in range(N):
        ind = random.randrange(t) # The index of the motif that will be ignored when making the profile.
        Profile = CreateProfile([*BestMotifs[:ind], *BestMotifs[ind+1:]])
        NewKMer = ChooseKMer(Dna[ind], k, Profile)
        NewMotifs = [*BestMotifs[:ind], NewKMer, *BestMotifs[ind+1:]]
        if Score(NewMotifs) < Score(BestMotifs):
            BestMotifs = NewMotifs
    return BestMotifs

        
def ChooseKMer(Sequence: str, k: int, Profile: list[dict[str,int]]) -> str:
    """
    Chooses which Kmer to put back to the best motif list based on Gibbs sampling.
    """
    Probabilities = KMerProbabilitiesOfStrand(Sequence, k, Profile)
    KMerStartIndex = ChooseKMerIndex(Probabilities)
    Motif = Sequence[KMerStartIndex: KMerStartIndex+k]
    return Motif
        
def KMerProbabilitiesOfStrand(Dna: str, k: int, Profile: list[dict[str, int]]) -> list[float]:
    """
    Calculates the chance of each KMer in a string occuring based on a profile.
    Returns the probabilities as a list of ints.

    Input:
        Dna (str): The string being analyzed.
        
        k (int): the length of the KMers.
    
        profile (list[dict[str,int]]): The profile used to calculate probabilites.
    
    Output:
        Probabilities list[int]: A list containing the probabilites of each kmer in Dna.
    """
    NumKMers = len(Dna)-k+1
    Probabilities = [0 for i in range(NumKMers)]
    for i in range(NumKMers):
        KMer = Dna[i:i+k]
        Probabilities[i] = ScoreKMer(KMer, Profile)
    
    # Normalize the Probabilities so they sum to 1
    NormalizationFactor = sum(Probabilities)
    Probabilities = [i/NormalizationFactor for i in Probabilities]
    return Probabilities


def ChooseKMerIndex(Probabilities: list[float]) -> int:
    """
    Chooses which Kmer to put back to the best motif list
    """
    DiceRoll = random.random()
    sum = 0
    for index, val in enumerate(Probabilities):
        sum += val
        if sum > DiceRoll:
            return index

if __name__ == "__main__":
    main()