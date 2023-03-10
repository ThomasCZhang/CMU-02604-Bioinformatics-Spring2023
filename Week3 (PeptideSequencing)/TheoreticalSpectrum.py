import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\Spectrum")
    filepaths = glob(dirpath + "\\data*2.txt")
    
    AminoAcidMassDict, AminoAcids = GenerateAAInfo()

    for path in filepaths:
        peptide_sequence = ReadTestInputs_Spectrum(path)
        peptide_mass = AAStringToMass(peptide_sequence, AminoAcidMassDict)
        spectrum = LinearSpectrum(peptide_mass)
        cycle_spectrum = CycloSpectrum(peptide_mass)

    AnswerPath = os.path.join(os.path.dirname(__file__), "LinearSpectrum.txt")
    with open(AnswerPath, "w") as f:
        for i, ele in enumerate(spectrum):
            if i == 0:
                f.write(str(ele))
            else:   
                f.write(" " + str(ele))


    AnswerPath = os.path.join(os.path.dirname(__file__), "CycloSpectrum.txt")

    with open(AnswerPath, "w") as f:
        for i, ele in enumerate(cycle_spectrum):
            if i == 0:
                f.write(str(ele))
            else:
                f.write(" " + str(ele))

def AAStringToMass(peptide_string: str, AminoAcidMass: dict[str, float]) -> list[int]:
    peptide = []
    for character in peptide_string:
        peptide.append(AminoAcidMass[character])
    return peptide


def GenerateAAInfo()-> tuple[dict[str, int], list[str]]:
    """
    CreateAADictionary: Creates and returns a dictionary for amino acid mass. 
    
    Output:
        AminoAcidMassDict: Dictionary containing masses of amino acid weights. Keys are single letter representation
        of amino acids.

        AminoAcids: A list of the single letter representations of amino acids.
    """
    AminoAcidMassDict = {
        "G": 57,
        "A": 71,
        "S": 87,
        "P": 97,
        "V": 99,
        "T": 101,
        "C": 103,
        "I": 113,
        "L": 113,
        "N": 114,
        "D": 115,
        "K": 128,
        "Q": 128,
        "E": 129,
        "M": 131,
        "H": 137,
        "F": 147,
        "R": 156,
        "Y": 163,
        "W": 186,
    }

    AminoAcids = [
        "G", "A", "S", "P", "V",
        "T", "C", "I", "L", "N",
        "D", "K", "Q", "E", "M",
        "H", "F", "R", "Y", "W"
    ]
    return AminoAcidMassDict, AminoAcids

def ReadTestInputs_Spectrum(FilePath: str) -> str:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        text(str): The peptide sequence (single letter amino acid).
    """
    with open(FilePath) as f:
        text = f.readline().strip()
    return text


def CycloSpectrum(Peptide: list[int]) -> list[float]:
    """
    CycloSpectrum takes a peptide sequence (str) as input and returns a list of all
    possible fragment lengths assuming the peptide is cyclical.

    Input: An amino acid string Peptide.

    Output: The cyclic spectrum of Peptide.
    """

    PrefixMass = [0 for i in range(len(Peptide)+1)]
    for i in range(1, len(Peptide)+1):
        PrefixMass[i] = PrefixMass[i-1] + Peptide[i-1]
    
    peptideMass = PrefixMass[len(Peptide)]
    cycle_spectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            cycle_spectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(Peptide):
                cycle_spectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
    return sorted(cycle_spectrum)
                


def LinearSpectrum(Peptide: list[int]) -> list[int]:
    """
    LinearSpectrum: Generates the theoretical mass spectrum generated by a peptide.

    Input:
        Peptide: An amino acid string Peptide.
    
        Alphabet: The single letter representations of the amino acids.

        AminoAcidMass: The masses of the amino acids stored in a dictionary. Keys are the single letter representations
        of amino acids.

    Output: 
        lin_spectrum: The linear spectrum of Peptide.
    """
    PrefixMass = [0 for i in range(len(Peptide)+1)]
    for i in range(1, len(Peptide)+1):
        PrefixMass[i] = PrefixMass[i-1] + Peptide[i-1]
    lin_spectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            lin_spectrum.append(PrefixMass[j]-PrefixMass[i])
    return sorted(lin_spectrum)

if __name__ == "__main__":
    main()