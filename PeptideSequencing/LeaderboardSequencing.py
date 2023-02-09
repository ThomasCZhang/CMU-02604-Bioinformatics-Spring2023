import os, copy
from CyclopeptideScoring import *
from CyclopeptideSequencing import *
from Trim import Trim
from custom_classes import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\LeaderboardSeq")
    filepaths = glob(dirpath + "\\input_0.txt")
    for filepath in filepaths:
        N, spectrum = ReadTestInputs_LeaderSeq(filepath)
        answer = LeaderboardSequencing(spectrum, N)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        for idx0, sp in enumerate(answer):
            if idx0 != 0:
                f.write(" ")
            for idx1, val in enumerate(sp.protein.peptide):
                if idx1 != 0:
                    f.write("-")
                f.write(str(val))

def ReadTestInputs_LeaderSeq(FilePath: str) -> tuple[int, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        N: the number of leader peptides to keep.

        spectrum: the mass spectrum.
    """
    with open(FilePath) as f:
        N = int(f.readline().strip())
        spectrum = [int(i) for i in f.readline().strip().split(" ")]

    return N, spectrum


def LeaderboardSequencing(spectrum: list[int], N: int) -> list[list[int]]:
    """
    CyclopeptideSequencing: Returns the sequence that generates.

    Input:
        spectrum: a list of masses.

    Ouput:
        final_peptides: a list of peptides that can generate specturm.
    """
    amino_acid_mass, _ = GenerateAAInfo()
    masses = set(list(amino_acid_mass.values()))

    leaderboard = [Scored_Protein([],"")]
    leader_peptides = []
    best_score = 0
    parent_mass = ParentMass(spectrum)
    gen = 0
    while leaderboard:
        print(f"\nCurrentGen: {gen}.")
        leaderboard = ExpandProtein(leaderboard, masses, spectrum)
        print(f"Leader Board Size: {len(leaderboard)}.")
        remove_list = []
        for ind, p in enumerate(leaderboard):
            print(f"\rLeaderboard while loop index: {ind}.", end = "")
            if p.protein.mass == parent_mass:
                current_score = CyclopeptideScoring(p.protein.peptide, spectrum)
                if current_score > best_score:
                    best_score = current_score
                    leader_peptides = [p]
                elif current_score == best_score:
                    leader_peptides.append(p)
            elif p.protein.mass > parent_mass:
                remove_list.append(ind)
        
        if len(remove_list) > 0:
            print()
        for i in range(len(remove_list)-1, -1, -1):
            print(f"\rRemoving Idx: {len(remove_list)-i-1}", end = "")
            idx = remove_list[i]
            del leaderboard[idx]

        leaderboard = Trim(leaderboard, spectrum, N)
        if len(leaderboard) != 0:
            print(f"\nMass of last peptide: {leaderboard[-1].protein.mass}")
        gen += 1
        
    
    print(f"\n\nNumber of peptides: {len(leader_peptides)}, \nBest Score: {best_score}")
    return leader_peptides


def ExpandProtein(proteins: list[Scored_Protein], masses: set[int], spectrum: list[int]) -> set[str]:
    """
    ExpandPeptide: Creates a list of all possible single letter peptide extensions of the all the peptide chains in 'peptides'.

    Input:
        protein: A protein object. These will be extended by each of the unique amino acid masses.

        masses: A set of unique amino acid masses.
        
    Output:
        new_peptides: The new peptides.
    """
    new_proteins=[]
    for prot in proteins:
        for mass in masses:
            temp_protein = copy.deepcopy(prot)
            temp_protein.AddAminoAcid(mass, spectrum)
            new_proteins.append(temp_protein)

    return new_proteins

if __name__ == "__main__":
    main()

