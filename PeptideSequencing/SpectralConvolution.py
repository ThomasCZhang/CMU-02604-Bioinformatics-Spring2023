import os
from glob import glob

def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\\Inputs\\convolution")
    filepaths = glob(dirpath + "\\data*.txt")
    for filepath in filepaths:
        spectrum = ReadTestInputs_Convolution(filepath)
        convolution_spectrum = SpectralConvolution(spectrum)

    convolution_spectrum.sort()
    answerpath = os.path.join(os.path.dirname(__file__), "new_answer.txt")
    with open(answerpath, 'w') as f:
        for idx, mass in enumerate(convolution_spectrum):
            if idx != 0:
                f.write(" ")
            f.write(str(mass))

def ReadTestInputs_Convolution(FilePath: str) -> list[int]:
    """
    Reads peptide sequence for theoretical specturm.
    
    Input:
        FilePath: Path to the test file.
    
    Output:       
        spectrum: The mass spectrum.
    """
    with open(FilePath) as f:
        spectrum = [int(i) for i in f.readline().strip().split(" ")]

    return spectrum


def SpectralConvolution(spectrum: list[int]) -> list[int]:
    """
    SpectralConvolution takes a spectrum as an input and returns its spectral convolution

    Input:
        spectrum: the given mass spectrum

    Output:
        convolution: The spectral convolution as a list of integers.
    """
    convolution_spectrum = []
    for i in range(len(spectrum)-1):
        for j in range(i+1, len(spectrum)):
            diff = abs(spectrum[j]-spectrum[i])
            if diff != 0:
                convolution_spectrum.append(diff)

    return convolution_spectrum

if __name__ == "__main__":
    main()