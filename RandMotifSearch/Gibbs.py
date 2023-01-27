import os
from glob import glob
from RandomizedMotifSearch import *


def GibbsSampler(Dna: list[str], k: int, t: int, N: int) -> list[str]:
    """
    Input: Integers k, t, and N, followed by a collection of strings Dna.
    Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
    """
    NumIterations = 20
    BestMotifs = ["" for i in range(t)]
    BestScore = k*t # Arbitrary large starting score 
    for i in range(NumIterations):
        Motifs = GibbsSamplerOneIteration(Dna, k, t, N)
        NewScore = Score(Motifs)
        if NewScore < BestScore:
            BestMotifs = Motifs
    return BestMotifs

def GibbsSamplerOneIteration(Dna: list[str], k: int, t: int, N: int) -> list[str]:
    """
    Finds the best motif reachable from a single random starting collection of motifs from Dna using Gibbs Sampling.\n
    Input:\n
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
        ind = random.randrange(k) # The index of the motif that will be ignored when making the profile.
        Profile = CreateProfile([*BestMotifs[:ind], *BestMotifs[ind+1:]])
        ChooseKMer(Dna[ind], k, Profile)

        
def ChooseKMer(Sequence: str, k: int, Profile: list[dict[str,int]]) -> str:
    """
    Chooses which Kmer to put back to the best motif list based on Gibbs sampling.
    """
    Probabilities = KMerProbabilitiesOfStrand(Sequence, k, Profile)
        
def KMerProbabilitiesOfStrand(Dna: str, k: int, Profile: list[dict[str, int]]) -> list[float])
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
        Probabilities[i] = ScoreKMer(Kmer, Profile)
    
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
        if val > DiceRoll:
            return index