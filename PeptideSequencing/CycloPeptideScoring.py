from TheoreticalSpectrum import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\scoring")
    filepaths = glob(dirpath + "\\dataset_876195_1.txt")
    for filepath in filepaths:
        peptide, spectrum = ReadTestInputs_PeptideScoring(filepath)
        score = LinearScoring(peptide, spectrum)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        f.write(str(score))

def ReadTestInputs_PeptideScoring(FilePath: str) -> tuple[str, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath: Path to the test file.
    
    Output:
        text: The peptide sequence (single letter amino acid).
        
        spectrum: The mass spectrum.
    """
    with open(FilePath) as f:
        text = f.readline().strip()
        spectrum = [int(i) for i in f.readline().strip().split(" ")]

    return text, spectrum

def CyclopeptideScoring(peptide: str, spectrum: list[int]) -> int:
    """
    CyclopeptideScoring: Returns the score of a peptides theoretical spectrum compared to a given spectrum. Assumes
    that the peptide is a cyclical peptide.

    Input:
        peptide: The given peptide. String of single letter amino acids.

        spectrum: The given spectrum.
    
    Output:
        score: the score.
    """
    aa_mass, alphabet = GenerateAAInfo()
    peptide_spectrum = CycloSpectrum(peptide, alphabet, aa_mass)

    score = 0
    used_indicies = [] # indexes of masses in spectrum that have already been matched with a mass in peptide_spectrum.
    for pep_mass in peptide_spectrum:
        for ind, spec_mass in enumerate(spectrum):
            if (pep_mass == spec_mass) and (ind not in used_indicies):
                used_indicies.append(ind)
                score += 1
                break
    return score

def LinearScoring(peptide: str, spectrum: list[int]) -> int:
    """
    CyclopeptideScoring: Returns the score of a peptides theoretical spectrum compared to a given spectrum.
    Assumes linear peptide

    Input:
        peptide: The given peptide. String of single letter amino acids.

        spectrum: The given spectrum.
    
    Output:
        score: the score.
    """
    aa_mass, alphabet = GenerateAAInfo()
    peptide_spectrum = LinearSpectrum(peptide, alphabet, aa_mass)

    score = 0
    used_indicies = [] # indexes of masses in spectrum that have already been matched with a mass in peptide_spectrum.
    for pep_mass in peptide_spectrum:
        for ind, spec_mass in enumerate(spectrum):
            if (pep_mass == spec_mass) and (ind not in used_indicies):
                used_indicies.append(ind)
                score += 1
                break
    return score

if __name__ == "__main__":
    main()