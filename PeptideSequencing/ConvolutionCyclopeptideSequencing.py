from SpectralConvolution import *
from LeaderboardSequencing import *


def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\ConvolutionSeq")
    filepaths = glob(dirpath + "\\input_1.txt")
    for filepath in filepaths:
        M, N, spectrum = ReadTestInputs_ConvolutionSeq(filepath)
        answer = ConvolutionCyclopeptideSequencing(spectrum, M, N)

    answer = answer[1:2]
    answerpath = os.path.join(os.path.dirname(__file__), "new_answer.txt")
    with open(answerpath, 'w') as f:
        for idx0, peptide in enumerate(answer):
            if idx0 != 0:
                f.write(" ")
            for idx1, val in enumerate(peptide):
                if idx1 != 0:
                    f.write("-")
                if idx1 < 38:
                    f.write(str(val))

def ReadTestInputs_ConvolutionSeq(FilePath: str) -> tuple[int, int, list[int]]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath(str): Path to the test file.
    
    Output:
        M: The number of most frequent masses to keep. 

        N: the number of leader peptides to keep.

        spectrum: the mass spectrum.
    """
    with open(FilePath) as f:
        M = int(f.readline().strip())
        N = int(f.readline().strip())
        spectrum = [int(i) for i in f.readline().strip().split(" ")]

    return M, N, spectrum


def ConvolutionCyclopeptideSequencing(spectrum: list[int], M: int, N: int) -> list[int]:
    """
    ConvolutionCyclopeptideSequencing: Performs LeaderboardCyclopeptideSequecing but with 
    the spectral convolution instead of the raw spectrum.

    Input:
        spectrum: The mass spectrum generated from data.

        M: the number of most frequent masses to be kept from the spectral convolution of spectrum.

        N: how many scores to keep for trimming in leaderboard.

    Output:
    """
    convolution_spectrum = SpectralConvolution(spectrum)
    mass_frequencies = CountMassFrequencies(convolution_spectrum)

    # aa_mass_dict, _ = GenerateAAInfo() # Not sure if we only want masses that correspond to the 20 common amino acids.
    del_list = []
    for mass in mass_frequencies:
        # if (mass < 57) or (mass > 200) or (mass not in aa_mass_dict.values()):
        if (mass < 57) or (mass > 200):
            del_list.append(mass)
    
    for mass in del_list:
        del mass_frequencies[mass]

    potential_masses = GetHighestFrequencyMasses(mass_frequencies, M)
    potential_masses.sort()
    peptides = ConvolutionLeaderboardSequencing(spectrum, N, potential_masses)
    return peptides

def GetHighestFrequencyMasses(frequency_dict: dict[int, int], M):
    """
    Returns a list of the most frequent M masses given a frequency dictionary where key = mass and value = frequency.

    Input:
        frequency_dict: The frequency dictionary.

    Ouput:
        best_masses: The M most frequent masses. Ties count for a single position.
    """
    highest_frequencies = {}
    for mass in frequency_dict:
        current_freq = frequency_dict[mass]
        if current_freq in highest_frequencies:
            highest_frequencies[current_freq].append(mass)
        elif len(highest_frequencies) < M: 
            highest_frequencies[current_freq] = [mass]
        else:
            for freq in highest_frequencies:
                if current_freq > freq:
                    del highest_frequencies[freq]
                    highest_frequencies[current_freq] = [mass]
                    break
    
    best_masses = [i for key in highest_frequencies for i in highest_frequencies[key]]
    return best_masses

def CountMassFrequencies(spectrum: list[int]) -> dict[int, int]:
    """
    CountMassFrequencies counts the number of times a mass appears in a mass spectrum.

    Input:
        spectrum: the given spectrum.

    Output:
        frequencies: the number of times each mass appears in spectrum. Stored as a dictionary
        where key = mass and value = # of appearences.
    """
    frequencies = {}
    for mass in spectrum:
        if mass in frequencies:
            mass += 1
        else:
            frequencies[mass] = 0
    return frequencies

def ConvolutionLeaderboardSequencing(spectrum: list[int], N: int, aa_list: list[int]) -> list[list[int]]:
    """
    CyclopeptideSequencing: Returns the sequence that generates.

    Input:
        spectrum: a list of masses.

        N: The number of top scores to keep when trimming

        aa_list: the allowed aa_masses when expanding the peptides.

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
        leaderboard = ConvolutionExpandPeptide(leaderboard, aa_list)
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

def ConvolutionExpandPeptide(mass_chains: set[str], aa_list: list[int]) -> set[str]:
    """
    ExpandPeptide: Creates a list of all possible single letter peptide extensions of the all the peptide chains in 'peptides'.

    Input:
        mass_chain: A list of peptide chains. These will be extended by each of amino acids in aa_list.

        aa_list: The masses that can be used to extend the peptide chains.

    Output:
        new_peptides: The new peptides.
    """
    new_masses=[]
    for mass_chain in mass_chains:
        for mass in aa_list:
            if mass_chain == [0]:
                new_masses.append([mass])
            else:
                new_masses.append([*mass_chain, mass])

    return new_masses

if __name__ == "__main__":
    main()