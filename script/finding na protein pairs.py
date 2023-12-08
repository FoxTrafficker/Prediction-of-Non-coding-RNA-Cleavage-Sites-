import json


def read_fasta(file_path):
    with open(file_path, "r") as file:
        sequences = []
        current_seq = ""
        for line in file:
            if line.startswith(">"):
                if current_seq:
                    sequences.append((header, current_seq))
                    current_seq = ""
                header = line[1:].strip()
            else:
                current_seq += line.strip()
        if current_seq:
            sequences.append((header, current_seq))
    return sequences


def extract_structure_name(sequence):
    first_word = sequence.split()[0]
    part_before_underscore = first_word.split("_")[0]
    return part_before_underscore


def identify_molecule(sequence):
    if "mol:na" in sequence:
        return "na"
    elif "mol:protein" in sequence:
        return "protein"
    else:
        return "Unknown"


def build_dictionary(fasta_file):
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


if __name__ == "__main__":
    fasta_file = "pdb_seqres.txt"  # 替换为你的 FASTA 文件路径
    dataset = build_dictionary(fasta_file)

    res = []

    for structure, types in dataset.items():
        print(f"Structure: {structure}, Types: {types}")
        if len(types["na"]) > 0 and len(types["protein"]) > 0:
            res.append(structure)
    print(res[:10])

    filename = "id_list.json"

    # 将列表写入 JSON 文件
    with open(filename, "w") as file:
        json.dump(res, file, indent=4)

    print(f"List has been saved to {filename}")
