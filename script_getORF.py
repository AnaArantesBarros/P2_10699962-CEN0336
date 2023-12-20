##Exercício 1
import sys
import re

# Criando dicionários e tabela de tradução
sequencias = {}
codons_frames = {}
longest_info = {}

translation_table = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'AAT': 'N', 'AAC': 'N',
    'GAT': 'D', 'GAC': 'D',
    'TGT': 'C', 'TGC': 'C',
    'CAA': 'Q', 'CAG': 'Q',
    'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'AAA': 'K', 'AAG': 'K',
    'ATG': 'M',
    'TTT': 'F', 'TTC': 'F',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TAA': '*', 'TGA': '*', 'TAG': '*'
}

# Ler fasta do input
file = sys.argv[1]

# Verificar se o arquivo está correto
try:
    # Garante que o formato do arquivo fornecido está correto.
    if not re.search(r"\.fasta",file):
        print("Forneça um arquivo do tipo fasta.")
        sys.exit(1)

    # Verifica se o número de argumentos é suficiente
    if len(sys.argv) < 1:
        raise ValueError("Forneça um único argumento.")

except FileNotFoundError:
    print(f'O arquivo {file} não foi encontrado.')


# Lê o arquivo fasta
with open(file) as fasta:
    geneID = ""
    protein_frames = {}
    coord_orfs = {}
    for linha in fasta: # Lê linha por linha do arquivo
        linha = linha.rstrip()
        if linha.startswith(">"):
            linhas_juntas = ""
            geneID = re.search(r"^>([\w]{8,10})", linha).group(1) #Pega o nome da sequência
        else:
            linhas_juntas += linha.strip() # Pega o restante das informações
            codons = re.findall(r".{3}", linhas_juntas.replace("\n", ""))
            sequencias[geneID] = "".join(codons) # Armazena as sequências em um dicionário

            if geneID not in codons_frames:
                codons_frames[geneID] = {} # Adiciona entradas no dicionário com todos os nomes de geneID

            for frame in range(6): # Lê os seis frames
                codons = re.findall(r".{3}", sequencias[geneID][frame:])
                frameID = "frame+" + str(frame + 1)
                codons_frames[geneID][frameID] = codons #Salva as sequências
                protein = "".join([translation_table[codon] for codon in codons]) #Traduz as sequências
                frame_aa = f"{geneID}_{frameID}_translated"
                if geneID not in protein_frames:
                    # Adiciona entradas no dicionário com todos os nomes de geneID
                    protein_frames[geneID] = {}
                    coord_orfs[geneID] = {}

                # Adiciona coordenadas do ORF mais longo
                orf_start = sequencias[geneID].find("".join(codons))
                orf_end = orf_start + len(codons)
                protein_frames[geneID][frame_aa] = protein
                coord_orfs[geneID][frame_aa] = {'start': orf_start, 'end': orf_end}


def get_posit(geneID,frame):
    for frame in protein_frames[geneID]:
        start = coord_orfs[geneID][frame]["start"]
        end = coord_orfs[geneID][frame]["end"]
        return start, end

with open("ORF.fna", "w") as fna_output, open("ORF.faa", "w") as faa_output:
    for geneID in protein_frames.keys():
        longestProtein = ""
        longestFrame = ""
        start = ""
        end = ""
        for frame in protein_frames[geneID]:
            proteinas = re.findall(r"(M[A-Z]+?)\*", protein_frames[geneID][frame])
            for i in proteinas:
                # Verifica qual dos frames tem a maior proteína
                if len(i) > len(longestProtein):
                    longestProtein = i
                    longestFrame = frame
                    longest_info[geneID] = longestFrame
        # Salva as coordenadas de inicio e fim
        start, end = get_posit(geneID,frame)
        # Salva o resultado para o arquivo
        faa_output.write(str(f">{longestFrame}_{start}_{end}\n").replace('_translated',""))
        faa_output.write(longestProtein + "\n")

    for geneID in longest_info.keys():
        start = ""
        end = ""
        frame_longest = longest_info[geneID].replace(f'{geneID}_', "").replace('_translated',"")
        longest_info[geneID] = frame_longest
        seq_logest = codons_frames[geneID][frame_longest]
        start, end = get_posit(geneID,frame)
        headline = f">{geneID}_{frame_longest}_{start}_{end}\n"
        codons_str = " ".join(seq_logest) + "\n"
        #Salva a sequência original da proteína mais longa
        fna_output.write(headline)
        fna_output.write(codons_str)


