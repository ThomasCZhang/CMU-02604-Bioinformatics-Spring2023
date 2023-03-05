import os
from glob import glob


def main():
    dirpath = os.path.join(os.path.dirname(__file__),
                           "Files\Inputs\PeptideEncoding")
    filepaths = glob(dirpath + "\\input*.txt")

    for path in filepaths:
        text, peptide = ReadTestInputs_PeptideEncoding(path)
        substrings = PeptideEncoding(text, peptide)
        print(f"File {os.path.basename(path)}.")
        total_starts = 0
        for substring in substrings:
            total_starts += 1
        print(total_starts)


codon_dict = {
    # Starting with U
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    # Starting with C
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    # Starting with A
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    # Starting with G
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

amino_acid_dict = {
    # Positive
    "Arg": "R", "His": "H", "Lys": "K",
    # Negative
    "Asp": "D", "Glu": "E",
    # Polar, Uncharged
    "Ser": "S", "Thr": "T", "Asn": "N", "Gln": "Q",
    # Hydrophobic
    "Ala": "A", "Val": "V", "Ile": "I", "Leu": "L", "Met": "M",
    "Phe": "F", "Tyr": "Y", "Trp": "W",
    # Special Cases
    "Cys": "C", "Gly": "G", "Pro": "P",
    # Stop
    "Stop": "*",
}


def ReadTestInputs_PeptideEncoding(FilePath: str) -> tuple[str, str]:
    """
    Reads the test files for PeptideEncoding.
    Returns a Dna sequence (text) and a peptide sequence (peptide) as strings.
    """
    with open(FilePath) as f:
        text = f.readline().strip()
        peptide = f.readline().strip()

    return text, peptide


def PeptideEncoding(text: str, peptide: str) -> list[str]:
    """
    Input: A DNA string Text, an amino acid string Peptide
    Output: All substrings of Text encoding Peptide (if any such substrings exist).
    """
    complement = CreateDnaComplement(text)
    peptide_length = len(peptide)
    # Number of starting positions possible on one strand.
    num_positions = len(text)-3*peptide_length+1
    results = []
    for i in range(num_positions):
        # Forward strand starting and ending index
        f_start = i
        f_end = f_start+3*peptide_length
        forward_match = CheckPeptideSequence(text[f_start:f_end], peptide)
        # Complementary strand starting and ending index
        c_start = num_positions-1-i
        c_end = c_start+3*peptide_length
        complement_match = CheckPeptideSequence(complement[c_start:c_end], peptide)
        if forward_match or complement_match:
            results.append(text[f_start:f_end])
    return(results)

def CheckPeptideSequence(text: str, peptide: str) -> bool:
    """
    Checks if a Dna sequence (text) codes for a amino acid sequence (peptide).
    Output: True if text encodes for peptide, false otherwise.
    """
    rna_text = ReplaceCharacter(text, "T", "U")
    for i in range(len(peptide)):
        codon_index = i*3
        # print(f"RNA length: {len(rna_text)} ; start: {codon_index} ; end: {codon_index+3}.")
        # print(f"RNA sequence: {rna_text}. Peptide Sequences: {peptide}.")
        codon = rna_text[codon_index: codon_index+3]
        amino_acid = TranslateRnaCodon(codon)
        if amino_acid != peptide[i]:
            return False
    return True

def CreateDnaComplement(text: str) -> str:
    """
    Creates the complementary DNA strand to a given DNA strand.
    Input:
    text (str): the DNA strand to find the complement of. Assumed to be written
    5' to 3'.
    Output:
    complement (str): The reverse DNA strand. Will be reversed so that the strand
    is written 3'to 5'.
    """
    complement_dict = {
        "A":"T", "T":"A", "G":"C", "C":"G"
    }

    n = len(text)-1
    complement = ""
    for i in range(n):
        complement += complement_dict[text[n-i].upper()]

    return complement



def ReplaceCharacter(text: str, old_char: str, new_char: str) -> str:
    """
    Replaces all instances of a character in a DNA string with another character. 
    Returns the new string
    """
    new_text_list = ["" for i in range(len(text))]
    for i in range(len(text)):
        if text[i] == old_char:
            new_text_list[i] = new_char
        else:
            new_text_list[i] = text[i]
    
    new_text = ""
    for ele in new_text_list:
        new_text += ele
    
    return new_text

def TranslateRnaCodon(codon: str) -> str:
    """
    Translates a three letter RNA codon to its amino acid.
    Returns the single letter code of the amino acid.
    """
    three_letter_aa = codon_dict[codon.upper()]
    single_letter_aa = amino_acid_dict[three_letter_aa]
    return single_letter_aa

def ReverseString(text: str) -> str:
    reversed_string = ""
    for i in range(len(text)):
        reversed_string += text[len(text)-1-i]
    return reversed_string

if __name__ == "__main__":
    main()