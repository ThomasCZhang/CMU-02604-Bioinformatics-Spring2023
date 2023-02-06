import os
from TheoreticalSpectrum import *
from CycloPeptideScoring import *
from CyclopeptideSequencing import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\LeaderboardSeq")
    filepaths = glob(dirpath + "\\input*.txt")
    for filepath in filepaths:
        N, spectrum = ReadTestInputs_LeaderSeq(filepath)
        answer = LeaderboardSequencing(spectrum, N)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        # for idx0, peptide in enumerate(answer):
        #     if idx0 != 0:
        #         f.write(" ")
            for idx, val in enumerate(answer):
                if idx != 0:
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


def LeaderboardSequencing(spectrum, N) -> list[int]:
    """
    CyclopeptideSequencing: Returns the sequence that generates.

    Input:
        spectrum: a list of masses.

    Ouput:
        final_peptides: a list of peptides that can generate specturm.
    """
    leaderboard = [[0]]
    leader_peptides = []
    while leaderboard:
        leaderboard = ExpandPeptide(leaderboard)
        remove_list = []
        for ind, peptide in enumerate(leaderboard):
            if CalculatePeptideMass(peptide) == ParentMass(spectrum):
                if CycloScore(peptide, spectrum) > CycloScore(leader_peptides, spectrum):
                    leader_peptides = peptide
                # elif CycloScore(peptide, spectrum) == CycloScore(leader_peptides, spectrum):
                #     leader_peptides.append(peptide)
            elif CalculatePeptideMass(peptide) > ParentMass(spectrum):
                remove_list.append(ind)

        for i in range(len(remove_list)-1, -1, -1):
            idx = remove_list[i]
            leaderboard = leaderboard[:idx] + leaderboard[idx+1:]

        leaderboard = Trim(leaderboard, spectrum, N)
                
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
    
    top_peptides_dict = {} # key = score. Values = peptides.
    top_scores = [-1 for i in range(N)]    

    for peptide in leaderboard:
        current_score = LinearScore(peptide, spectrum)
        for idx, ele in enumerate(top_scores):
            if current_score > ele:
                # Delete the tracker for the old best score
                if top_scores[idx] in top_peptides_dict:
                    del top_peptides_dict[top_scores[idx]]
                top_scores[idx] = current_score
                top_peptides_dict[current_score] = [peptide]

            elif (current_score == ele) and (peptide not in top_peptides_dict[ele]):
                top_peptides_dict[ele].append(peptide)

    top_peptides = []
    for peptide_lists in top_peptides_dict.values():
        top_peptides = top_peptides + peptide_lists

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