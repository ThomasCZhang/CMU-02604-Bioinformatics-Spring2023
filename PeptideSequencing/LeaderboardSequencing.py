import os
from TheoreticalSpectrum import *
from CycloPeptideScoring import *
from CyclopeptideSequencing import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\LeaderboardSeq")
    filepaths = glob(dirpath + "\\input_0.txt")
    for filepath in filepaths:
        N, spectrum = ReadTestInputs_LeaderSeq(filepath)
        answer = LeaderboardSequencing(spectrum, N)

    answerpath = os.path.join(os.path.dirname(__file__), "answer2.txt")
    with open(answerpath, 'w') as f:
        for idx0, peptide in enumerate(answer):
            if idx0 != 0:
                f.write(" ")
            for idx1, val in enumerate(peptide):
                if idx1 != 0:
                    f.write("-")
                if idx1 < 38:
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
    leaderboard = [[0]]
    leader_peptides = [[]]
    best_score = 0
    parent_mass = ParentMass(spectrum)
    gen = 0
    while leaderboard:
        print(f"\nCurrentGen: {gen}.")
        leaderboard = ExpandPeptide(leaderboard)
        print(f"Leader Board Size: {len(leaderboard)}.")
        remove_list = []
        for ind, peptide in enumerate(leaderboard):
            print(f"\rLeaderboard while loop index: {ind}.", end = "")
            peptide_mass = CalculatePeptideMass(peptide)
            if peptide_mass == parent_mass:
                current_score = CycloScore(peptide, spectrum)
                if current_score > best_score:
                    best_score = current_score
                    leader_peptides = [peptide]
                elif current_score == best_score:
                    leader_peptides.append(peptide)
                    
            elif peptide_mass > parent_mass:
                remove_list.append(ind)
        if len(remove_list) > 0:
            print()
        for i in range(len(remove_list)-1, -1, -1):
            print(f"\rRemoving Idx: {len(remove_list)-i-1}", end = "")
            idx = remove_list[i]
            del leaderboard[idx]

        leaderboard = Trim(leaderboard, spectrum, N)
        if len(leaderboard) != 0:
            print(f"\nMass of last peptide: {CalculatePeptideMass(leaderboard[-1])}")
        gen += 1
        
    
    print(f"\nNumber of peptides: {len(leader_peptides)}, \nBest Score: {best_score}")
    return leader_peptides

def Trim(leaderboard: list[list[int]], spectrum: list[int],N: int) -> list[list[int]]:
    """
    Trim: Takes a list of peptides, a spectrum and an integer N. Returns the top N scoring peptides based on the
    given spectrum.

    Input:
        leaderboard: list of peptides. each peptide represented as a list of its amino acid masses in order.

        specturm: mass spectrum as a list of ints.

        N: the number of candidates to keep.
    
    Output:
        top_peptides: the top N scoring peptides.
    """
    
    # top_scores = [-1 for i in range(N)]
    top_scores = {} # key = score. value = index of peptide in leaderboard.
    print("")
    for idx, peptide in enumerate(leaderboard):
        print(f"\rTrimming Index: {idx}", end = "")
        current_score = LinearScore(peptide, spectrum)
        if current_score in top_scores:
            top_scores[current_score].append(idx)
        elif len(top_scores) < N:
            top_scores[current_score] = [idx]
        else:
            for score in top_scores:
                if current_score > score:
                    del top_scores[score]
                    top_scores[current_score] = [idx]
                    break
    
    best_scores = list(top_scores)
    best_scores.sort(reverse=True)

    top_peptides = [leaderboard[i] for score in best_scores for i in top_scores[score]]
    if len(top_peptides) > 6000:
        top_peptides = top_peptides[:3000]

    return top_peptides

def LinearScore(peptide: list[int], spectrum: list[int]) -> int:
    """
    Score: Calculates the score between the linear spectrum of a peptide and a given spectrum.

    Input:
        peptide: a peptide represented as a list of peptide masses.

    Output:
        spectrum: the given mass spectrum.
    """
    peptide_spectrum = LinearSpectrumMass(peptide)
    score = 0
    used_indicies = set() # indexes of masses in spectrum that have already been matched with a mass in peptide_spectrum.
    for pep_mass in peptide_spectrum:
        for ind, spec_mass in enumerate(spectrum):
            if (pep_mass == spec_mass) and (ind not in used_indicies):
                used_indicies.add(ind)
                score += 1
                break
    return score

def CycloScore(peptide: list[int], spectrum: list[int]) -> int:
    """
    Score: Calculates the score between the linear spectrum of a peptide and a given spectrum.

    Input:
        peptide: a peptide represented as a list of peptide masses.

    Output:
        spectrum: the given mass spectrum.
    """
    peptide_spectrum = CycloSpectrumMass(peptide)
    score = 0
    used_indicies = set() # indexes of masses in spectrum that have already been matched with a mass in peptide_spectrum.
    for pep_mass in peptide_spectrum:
        for ind, spec_mass in enumerate(spectrum):
            if (pep_mass == spec_mass) and (ind not in used_indicies):
                used_indicies.add(ind)
                score += 1
                break
    return score

if __name__ == "__main__":
    main()