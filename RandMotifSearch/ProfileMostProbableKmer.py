import os
from glob import glob


def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\Inputs\ProfileMostProbableKmer")
    filepaths = glob(dirpath + "\\*.txt")
    for path in filepaths:
        Text, k, Profile = ReadTestFiles_ProfileMostProbableKmer(path)
        bestMotif = ProfileMostProbableKmer(Text, k, Profile)
        print(
            f"The best motif for {os.path.basename(path)} is {bestMotif} with a score of {ScoreKMer(bestMotif, Profile): 0.2e}.")


def ReadTestFiles_ProfileMostProbableKmer(filepath: str) -> tuple[str, int, list[dict[str, float]]]:
    """
    ReadData(): Reads data from a .txt file and formats that data for ProfileMostProbableKmer.
    .txt Format should be as follows:\n
    Line 1: string of "A", "T", "C", "G".\n
    Line 2: an integer (value of k).\n
    Line 3 to Line 6: A, C, G, T profiles (in that order).\n
    Input:\n
    filepath (string): The path to the .txt file.\n
    Output:\n
    Text (string): the text from which to find the best motif.\n
    k (int): the number of letters in the motif.\n
    Profile: (list of dictionaries): the scoring profile for the motifs.
    """
    switchDict = {
        0: "A",
        1: "C",
        2: "G",
        3: "T",
    }
    with open(filepath) as f:
        for ind, line in enumerate(f):
            if ind == 0:
                Text = line.strip()
            elif ind == 1:
                k = int(line.strip())
                # Profile = np.zeros((4, k))
                Profile = [{} for i in range(0, k)]
            else:
                values = [float(val) for val in line.strip().split(" ")]
                # Profile[ind-2, :] = np.array(values)
                for a, b in enumerate(values):
                    Profile[a][switchDict[ind-2]] = b
    return (Text, k, Profile)


def ProfileMostProbableKmer(Text: str, k: int, Profile: list[dict[str, float]]) -> str:
    """
    Takes a Text (string), a number k (int), and a 4 by k matrix profile (numpy matrix).\n
    Input:\n
    Text (string): The string from which to find the best motif.\n
    k (int): The number of letters in the motif.\n
    Profile (list of dictionaries): The scoring matrix.\n
    Output:\n
    bestMotif (string): The best motif.
    """
    # For i from 0 to last starting position
    #   score the k-mer starting from position i in Text using the scoring matrix "Profile"
    bestMotif = Text[0:k]
    currentBestScore = ScoreKMer(bestMotif, Profile)
    for i in range(1, len(Text)-k+1):
        currentMotif = Text[i:i+k]
        score = ScoreKMer(currentMotif, Profile)
        if score > currentBestScore:
            currentBestScore = score
            bestMotif = currentMotif
    return bestMotif


def ScoreKMer(motif: str, Profile: list[dict[str, float]]) -> float:
    """
    ScoreKmer() Takes a string of length k and then scores it using a 4 by k matrix Profile.\n
    Input:\n
    motif (string): The string to be scored.\n
    Profile (list of dictionaries): The matrix profile used to score the string.\n
    Output:\n
    score (float): The score of the motif (value should be between 0 and 1).
    """
    for index, element in enumerate(motif):
        if index == 0:
            currentScore = Profile[index][element.upper()]
        else:
            currentScore *= Profile[index][element.upper()]
    return currentScore


if __name__ == "__main__":
    main()
