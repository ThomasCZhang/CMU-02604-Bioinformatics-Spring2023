from RandomizedMotifSearch import *

def CreateHiddenMatrix(Dna, Profile) -> list[list[float]]:
    """
    CreateHiddenMatrix Creates a hidden matrix from a list of Dna strings and a
    Profile matrix.
    """
    k = len(Profile)
    n = len(Dna)
    HiddenMatrix = [[] for i in Dna]
    for i in range(n):
        HiddenMatrix[i] = CreateHiddenMatrixRow(Dna[i], Profile)
    NormalizeHiddenMatrix(HiddenMatrix)
    return HiddenMatrix

def NormalizeHiddenMatrix(Matrix: list[list[float]]):
    for i in range(len(Matrix)):
        NormalizationFactor = sum(Matrix[i])
        Matrix[i] = [x/NormalizationFactor for x in Matrix[i]]

def CreateHiddenMatrixRow(Sequence, Profile) -> list[float]:
    n = len(Sequence)
    k = len(Profile)
    Hidden_Matrix_Row = [0.0 for i in range(n-k+1)]
    for i in range(0, n-k+1):
        Motif = Sequence[i:i+k]
        Hidden_Matrix_Row[i] = ScoreKMer(Motif, Profile)
    return Hidden_Matrix_Row

def CreateProfileWithHiddenMatrix(Dna, HiddenMatrix, k) -> list[dict[str, int]]:
    Profile = [{"A": 0, "C": 0, "G": 0, "T": 0} for i in range(k)]
    for i, Strand in enumerate(Dna):
        for j in range(len(Strand)-k+1):
            K_mer = Strand[j:j+k]
            WeightedAdditionToProfile(K_mer, HiddenMatrix[i][j], Profile)

    Pseudocount = 1
    for Dictionary in Profile:
        AddConstantToDictionaryValues(Pseudocount, Dictionary)
        NormalizeDictionaryOfNumbers(Dictionary)
    
    return Profile

def WeightedAdditionToProfile(K_mer, Weight, Profile):
    for ind, letter in enumerate(K_mer):
        Profile[ind][letter.upper()] += Weight

def NormalizeDictionaryOfNumbers(Dictionary):
    total = 0
    for Key in Dictionary:
        total += Dictionary[Key]
    
    for Key in Dictionary:
        Dictionary[Key] /= total