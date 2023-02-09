import os
from glob import glob
from TheoreticalSpectrum import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\CyclopeptideSequencing")
    filepaths = glob(dirpath + "\\dataset*.txt")
    for filepath in filepaths:
        spectrum = ReadTestInputs_PeptideEncoding(filepath)
        peptideList = CyclopeptideSequencing(spectrum)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        for idx1, peptide in enumerate(peptideList):
            if idx1 != 0:
                f.write(" ")
            for idx2, mass in enumerate(peptide):
                if idx2 != 0:
                    f.write("-")
                f.write(str(mass))

def ReadTestInputs_PeptideEncoding(FilePath: str) -> list[int]:
    """
    Reads the test files for PeptideEncoding.
    Returns a Dna sequence (text) and a peptide sequence (peptide) as strings.
    """
    with open(FilePath) as f:
        text = [int(i) for i in f.readline().strip().split()]
    return text

def CyclopeptideSequencing(spectrum: list[int]) -> list[str]:
    """
    CyclopeptideSequencing: Returns the sequence that generates.

    Input:
        spectrum: a list of masses.

    Ouput:
        final_peptides: a list of peptides that can generate specturm.
    """
    parent_mass = ParentMass(spectrum)
    candidate_peptides = [[0]]
    final_peptides = []
    while candidate_peptides:
        candidate_peptides = ExpandPeptide(candidate_peptides)
        remove_list = []
        for ind, peptide in enumerate(candidate_peptides):
            peptide_spectrum = LinearSpectrum(peptide)
            if CalculatePeptideMass(peptide) == parent_mass:
                peptide_spectrum = CycloSpectrum(peptide)
                if CheckSpectrumCompatibility(peptide_spectrum, spectrum) and (peptide not in final_peptides):
                    final_peptides.append(peptide)
                remove_list.append(ind)
            elif not CheckSpectrumCompatibility(peptide_spectrum, spectrum):
                remove_list.append(ind)
        
        for i in range(len(remove_list)-1, -1, -1):
            j = remove_list[i]
            candidate_peptides = candidate_peptides[:j] + candidate_peptides[j+1:]
            
    return final_peptides


def ExpandPeptide(mass_chains: set[str]) -> set[str]:
    """
    ExpandPeptide: Creates a list of all possible single letter peptide extensions of the all the peptide chains in 'peptides'.

    Input:
        peptides: A list of peptide chains. These will be extended by each of the 20 amino acids.

    Output:
        new_peptides: The new peptides.
    """
    amino_acid_mass, _ = GenerateAAInfo()
    masses = set(list(amino_acid_mass.values()))
    new_masses=[]
    for mass_chain in mass_chains:
        for mass in masses:
            if mass_chain == [0]:
                new_masses.append([mass])
            else:
                new_masses.append([*mass_chain, mass])

    return new_masses

def CalculatePeptideMass(peptide: list[int]) -> int:
    """
    CalculatePeptideMass. calculates the mass of a peptide string. Assumes peptides are represented by single letters.

    Input:
        peptide: A peptide represented as a list of aa masses

    Output:
        mass: the mass of peptide.
    """
    return sum(peptide)

def ParentMass(spectrum: list[int]) -> int:
    """
    ParentMass: Returns the parent mass of a mass spectrum.

    Input:
        spectrum: mass spectrum data as a list of ints.

    Output:
        parent_mass: The parent mass of spectrum.
    """
    parent_mass = 0
    for mass in spectrum:
        if mass > parent_mass:
            parent_mass = mass
    return parent_mass

def CheckSpectrumCompatibility(peptide_spectrum: list[int], spectrum: list[int]) -> bool:
    """
    CheckSpectrumCompatibility: Checks if the spectrum generated by peptide is compatible with a given spectrum.

    Input:
        peptide: The peptide as a string of single letter amino acids.

        spectrum: The mass spectrum being used to determine whether peptide is compatible or not.
    
    Output:
        bool: True if peptide is compatible. False otherwise.
    """

    used_indicies = [] # indexes of masses in spectrum that have already been matched with a mass in peptide_spectrum.
    for pep_mass in peptide_spectrum:
        found_match = False
        for ind, spec_mass in enumerate(spectrum):
            if (pep_mass == spec_mass) and (ind not in used_indicies):
                used_indicies.append(ind)
                found_match = True
                break
        if not found_match: # If there is a unmatched mass between the theoretical spectrum and given spectrum.
            return False

    return True

if __name__ == "__main__":
    main()