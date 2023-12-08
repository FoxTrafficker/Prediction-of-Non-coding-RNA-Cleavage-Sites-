import json


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of tuples containing the headers and sequences.

    Parameters
    ----------
    file_path : str
        The path to the FASTA file.

    Returns
    -------
    list of tuples
        A list where each tuple contains two elements:
        1. The header (description) of the sequence (str).
        2. The sequence itself (str), which could be RNA, DNA, or protein.

    Notes
    -----
    The FASTA file format is a text-based format for representing nucleotide sequences or peptide sequences,
    in which nucleotides or amino acids are represented using single-letter codes. Each sequence in the file
    is introduced by a line starting with '>', followed by lines of sequence data.

    Example
    -------
    fasta_file.fa :
        >101d_B mol:na length:12  DNA (5'-D(*CP*GP*CP*GP*AP*AP*TP*TP*(CBR)P*GP*CP*G)-3')
        CGCGAATTCGCG
        >101m_A mol:protein length:154  MYOGLOBIN
        MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKK
        GHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG

    >>> read_fasta("fasta_file.fa ")
    [
        (
            "101d_B mol:na length:12  DNA (5'-D(*CP*GP*CP*GP*AP*AP*TP*TP*(CBR)P*GP*CP*G)-3')",
            "AGCT"
        ),
        (
            "101m_A mol:protein length:154  MYOGLOBIN",
            "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG"
        )
    ]
    """
    with open(file_path, "r") as file:
        header_sequences_paires = []
        current_seq = ""
        for line in file:
            if line.startswith(">"):
                if current_seq:
                    header_sequences_paires.append((header, current_seq))
                    current_seq = ""
                header = line[1:].strip()
            else:
                current_seq += line.strip()
        if current_seq:
            header_sequences_paires.append((header, current_seq))
    return header_sequences_paires


def extract_structure_name(sequence):
    """
    Extract and return the first part of the structure name from a sequence description.

    This function takes the first word of the sequence description and then returns
    the part of the word before the underscore. If there is no underscore in the description,
    it returns the entire first word.

    Parameters
    ----------
    sequence : str
        A string describing the sequence, typically a line from a FASTA file.

    Returns
    -------
    str
        The first part of the extracted structure name. If there is no underscore in the description,
        the entire first word is returned.

    Examples
    --------
    >>> extract_structure_name(">101m_A mol:protein length:154  MYOGLOBIN")
    '101m'

    >>> extract_structure_name(">157d_A mol:na length:12  RNA (5'-R(*CP*GP*CP*GP*AP*AP*UP*UP*AP*GP*CP*G)-3')")
    '157d'

    >>> extract_structure_name(">SEQUENCE_ABC_DEF")
    'SEQUENCE'

    """
    first_word = sequence.split()[0]
    part_before_underscore = first_word.split("_")[0]
    return part_before_underscore


def identify_molecule(sequence):
    """
    Identify the type of molecule from a given fasta discription.

    Parameters
    ----------
    sequence : str
        A discription representing the molecular sequence, which should contain
        specific substrings ('mol:na' or 'mol:protein') to indicate the molecule type.

    Returns
    -------
    str
        The type of molecule identified from the discription. This can be either 'na'
        (indicating RNA or DNA), 'protein', or 'Unknown' if the molecule type
        cannot be determined from the sequence.

    Examples
    --------
    >>> identify_molecule(">157d_A mol:na length:12  RNA (5'-R(*CP*GP*CP*GP*AP*AP*UP*UP*AP*GP*CP*G)-3')")
    'na'

    >>> identify_molecule(">101m_A mol:protein length:154  MYOGLOBIN")
    'protein'

    >>> identify_molecule("RandomSequenceWithoutIdentifier")
    'Unknown'
    """
    if "mol:na" in sequence:
        return "na"
    elif "mol:protein" in sequence:
        return "protein"
    else:
        return "Unknown"


def build_dictionary_from_fasta(fasta_file):
    """
    Construct a dictionary from a FASTA file, categorizing sequences by their structure name and molecule type.

    Parameters
    ----------
    fasta_file : str
        The path to the FASTA file to be read.

    Returns
    -------
    dict
        A dictionary where each key is a structure name extracted from the FASTA file. The value for each key is another
        dictionary with keys 'na', 'protein', and 'Unknown', each of which is associated with a list of sequences
        belonging to that category.

    Notes
    -----
    The function relies on the `read_fasta` function to read sequences from the FASTA file and on two additional
    functions: `extract_structure_name` and `identify_molecule`. `extract_structure_name` is used to extract the
    structure name from the sequence header, and `identify_molecule` is used to determine the molecule type
    ('na' or 'protein'). If the molecule type cannot be determined, it is categorized as 'Unknown'.

    Example
    -------
    Assuming `fasta_file` is a path to a valid FASTA file:

    >>> result = build_dictionary(fasta_file)
    >>> print(result)
    {'structure1': {'na': ['sequence1', 'sequence3'], 'protein': ['sequence2'], 'Unknown': []},
     'structure2': {'na': [], 'protein': ['sequence4'], 'Unknown': ['sequence5']}
    }
    """
    sequences = read_fasta(fasta_file)
    result = {}
    for header, seq in sequences:
        structure_name = extract_structure_name(header)
        molecule_type = identify_molecule(header)

        if structure_name not in result:
            result[structure_name] = {"na": [], "protein": [], "Unknown": []}

        if molecule_type in result[structure_name]:
            result[structure_name][molecule_type].append(seq)
        else:
            result[structure_name]["Unknown"].append(seq)

    return result


def find_na_protein_pairs_id(data_dict):
    res = []
    for id, types in data_dict.items():
        print(f"Structure: {id}, Types: {types}")
        if len(types["na"]) > 0 and len(types["protein"]) > 0:
            res.append(id)
    return res


def find_protein_only_id(data_dict):
    res = []
    for id, types in data_dict.items():
        print(f"Structure: {id}, Types: {types}")
        if len(types["na"]) == 0 and len(types["protein"]) > 0:
            res.append(id)
    return res


def json_write_back(res, filename="id_list.json"):
    with open(filename, "w") as file:
        json.dump(res, file, indent=4)
    print(f"List has been saved to {filename}")


if __name__ == "__main__":
    fasta_file = "/Users/dingsongting/VscodeProjects/Prediction-of-Non-coding-RNA-Cleavage-Sites-/script/pdb_seqres.txt"
    data_dict = build_dictionary_from_fasta(fasta_file)
    res = finding_na_protein_pairs(data_dict)
    json_write_back(res)
