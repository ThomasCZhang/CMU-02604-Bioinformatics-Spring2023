from SpectralConvolution import *
from LeaderboardSequencing import *

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\ConvolutionSeq")
    filepaths = glob(dirpath + "\\real*.txt")
    for filepath in filepaths:
        M, N, spectrum = ReadTestInputs_ConvolutionSeq(filepath)
        answer = ConvolutionCyclopeptideSequencing(spectrum, M, N)

    answerpath = os.path.join(os.path.dirname(__file__), "answer.txt")
    with open(answerpath, 'w') as f:
        for idx0, sp in enumerate(answer):
            if idx0 != 0:
                f.write("\n")
            for idx1, val in enumerate(sp.protein.peptide):
                if idx1 != 0:
                    f.write(" ")
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
        spectrum = [int(round(float(i))) for i in f.readline().strip().split(" ")]

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

    del_list = []
    for mass in mass_frequencies:
        if (mass < 57) or (mass > 200):
            del_list.append(mass)
    
    for mass in del_list:
        del mass_frequencies[mass]

    potential_masses = GetHighestFrequencyMasses(mass_frequencies, M)
    potential_masses.sort()
    peptides = ConvolutionLeaderboardSequencing(spectrum, N, potential_masses)
    return peptides

def GetHighestFrequencyMasses(mass_frequency_dict: dict[int, int], M):
    """
    Returns a list of the most frequent M masses given a frequency dictionary where key = mass and value = frequency.

    Input:
        frequency_dict: The frequency dictionary. Key = mass, value = frequency.

    Ouput:
        best_masses: The M most frequent masses. 
    """
    frequencies = {}
    for mass in mass_frequency_dict:
        current_freq = mass_frequency_dict[mass]
        if current_freq in frequencies:
            frequencies[current_freq].append(mass)
        else: 
            frequencies[current_freq] = [mass]

    all_freqs = list(frequencies)
    all_freqs.sort(reverse = True)

    best_masses = []
    for freq in all_freqs:
        best_masses.extend(frequencies[freq])
        if len(best_masses) >= M:
            break    
    
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
            frequencies[mass] += 1
        else:
            frequencies[mass] = 1
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
    leaderboard = [Scored_Protein([], "")]
    leader_peptides = []
    best_score = 0
    parent_mass = ParentMass(spectrum)
    gen = 0
    while leaderboard:
        print(f"\rCurrentGen: {gen}.", end = " ")
        leaderboard = ExpandProtein(leaderboard, aa_list, spectrum)
        remove_list = []
        for ind, p in enumerate(leaderboard):
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
            idx = remove_list[i]
            del leaderboard[idx]

        leaderboard = Trim(leaderboard, spectrum, N)
        gen += 1
        
    
    print(f"\n\nNumber of peptides: {len(leader_peptides)}, \nBest Score: {best_score}")
    return leader_peptides

if __name__ == "__main__":
    main()