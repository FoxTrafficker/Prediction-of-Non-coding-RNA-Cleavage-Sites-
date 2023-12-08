from Bio.Seq import Seq
from Bio.Data import CodonTable


def reverse_translate(protein_sequence):
    # 使用标准遗传密码表
    table = CodonTable.standard_dna_table.forward_table

    # 反向翻译的核苷酸序列
    nucleotide_sequence = ""

    for amino_acid in protein_sequence:
        # 查找第一个匹配的密码子
        for codon, aa in table.items():
            if aa == amino_acid:
                nucleotide_sequence += codon
                break
        else:
            # 如果找不到匹配的氨基酸（例如遇到终止符），添加一个占位符
            nucleotide_sequence += "NNN"

    return nucleotide_sequence


if __name__ == "__main__":
    # 示例蛋白质序列
    protein_seq = Seq("MSTNG")

    # 执行反向翻译
    dna_seq = reverse_translate(protein_seq)

    print(dna_seq)
