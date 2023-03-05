import os
from glob import glob
from custom_classes import *
from TheoreticalSpectrum import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\Trim")
    filepaths = glob(dirpath + "\\data*.txt")
    aa_masses,_ = GenerateAAInfo()
    for filepath in filepaths:
        peptides, spectrum, N = ReadTestInputs_Trim(filepath)
        sp = []
        for p in peptides:
            p_masses = AAStringToMass(p, aa_masses)
            protein = Scored_Protein(p_masses, p)
            protein.LinearScore(spectrum)
            sp.append(protein)
        answer = Trim(sp, spectrum, N)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        for idx0, sp in enumerate(answer):
            if idx0 != 0:
                f.write(" ")
            f.write(sp.protein.aa)

def ReadTestInputs_Trim(FilePath: str) -> tuple[int, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        peptides: A list of peptide strings.

        spectrum: the mass spectrum.

        N: the number of leader peptides to keep.
    """
    with open(FilePath) as f:
        peptides = f.readline().strip().split(" ")
        spectrum = [int(i) for i in f.readline().strip().split(" ")]
        N = int(f.readline().strip())

    return peptides, spectrum, N

def Trim(leaderboard: list[Scored_Protein], spectrum: list[int],N: int) -> list[list[int]]:
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
    
    all_scores = {} # key = score. value = index of peptide in leaderboard.
    for idx, p in enumerate(leaderboard):
        if p.score in all_scores:
            all_scores[p.score].append(idx)
        else:
            all_scores[p.score] = [idx]
    
    sorted_scores = list(all_scores)
    sorted_scores.sort(reverse=True)

    top_peptides = []
    for score in sorted_scores:
        for idx in all_scores[score]:
            top_peptides.append(leaderboard[idx])
        if len(top_peptides) >= N:
            break

    return top_peptides

if __name__ == "__main__":
    main()