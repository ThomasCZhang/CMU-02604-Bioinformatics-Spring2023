from TheoreticalSpectrum import *

def CyclopeptideSequencing(Spectrum list[int]) -> list[str]:
    candidate_peptides = [""]
    final_peptides = []
    while not candidate_peptides:

    pass

def Expand(peptides: list[str]) -> list[str]:
    """
    Expand: Creates a list of all possible single letter peptide extensions of the all the peptide chains in 'peptides'.

    Input:
        peptides: A list of peptide chains. These will be extended by each of the 20 amino acids.

    Output:
        new_peptides: The new peptides.
    """
    _, AminoAcids = GenerateAAInfo()
    new_peptides=[]
    for peptide in peptides:
        for aa in AminoAcids:
            new_peptides.append(peptide+aa)
    return new_peptides



CyclopeptideSequencing(Spectrum)
    CandidatePeptides ← a set containing only the empty peptide
    FinalPeptides ← empty list of strings
    while CandidatePeptides is nonempty
        CandidatePeptides ← Expand(CandidatePeptides)
        for each peptide Peptide in CandidatePeptides
            if Mass(Peptide) = ParentMass(Spectrum)
                if Cyclospectrum(Peptide) = Spectrum and Peptide is not in FinalPeptides
                    append Peptide to FinalPeptides
                remove Peptide from CandidatePeptides
            else if Peptide is not consistent with Spectrum
                remove Peptide from CandidatePeptides
    return FinalPeptides

Expand