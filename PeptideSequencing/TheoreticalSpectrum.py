import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\Inputs\Spectrum")
    filepaths = glob(dirpath + "\\input*.txt")
    
    for path in filepaths:
        peptide_sequence = ReadTestInputs_Spectrum(path)
        spectrum = LinearSpectrum(peptide_sequence, AminoAcids, AminoAcidMassDict)
        cycle_spectrum = CycloSpectrum(peptide_sequence, AminoAcids, AminoAcidMassDict)
        
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


def CycloSpectrum(Peptide: str, Alphabet: list[str], AminoAcidMass: dict[str, float]) -> list[float]:
    """
    CycloSpectrum takes a peptide sequence (str) as input and returns a list of all
    possible fragment lengths assuming the peptide is cyclical.

    Input: An amino acid string Peptide.

    Output: The cyclic spectrum of Peptide.
    """
    PrefixMass = [0 for i in range(len(Peptide)+1)]
    for i in range(1, len(Peptide)+1):
        for _, ele in enumerate(Alphabet):
            if ele == Peptide[i-1]:
                PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[ele]
    peptideMass = PrefixMass[len(Peptide)]
    cycle_spectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            cycle_spectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(Peptide):
                cycle_spectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
    return sorted(cycle_spectrum)
                


def LinearSpectrum(Peptide: str, Alphabet: list[str], AminoAcidMass: dict[str, float]) -> list[int]:
    """
    Input: An amino acid string Peptide.
    Output: The linear spectrum of Peptide.
        """
    PrefixMass = [0 for i in range(len(Peptide)+1)]
    for i in range(1, len(Peptide)+1):
        for _, ele in enumerate(Alphabet):
            if ele == Peptide[i-1]:
                PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[ele]
    lin_spectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            lin_spectrum.append(PrefixMass[j]-PrefixMass[i])
    return sorted(lin_spectrum)

if __name__ == "__main__":
    main()