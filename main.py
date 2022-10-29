'''
Index using letters and table
U = 0
C = 1
A = 2
G = 3
First letter -> Second letter -> Third letter
e.g: UUU: [0][0][0]
     CUU: [1][0][0]
     CAU: [1][2][0]

[[UUU,UUC,UUA,UUG][UCU,UCC,UCA,UCG][,,,][,,,]]
[[CUU,,,][,,,][,,,][,,,]]
[[,,,][,,,][,,,][,,,]]
[[,,,][,,,][,,,][,,,]]
'''

#Examples
dna_seq = "ATGTGGCCAGTCAUUUGA"
dna_seq2 = "AUGUUUCCCGGGAAAGAUCUG"

amino_acids = {
    "Ala":"Alanine",
    "Cys":"Cysteine",
    "Asp":"Aspartic Acid",
    "Glu":"Glutamic Acid",
    "Phe":"Phenylalanine",
    "Gly":"Glycine",
    "His":"Histidine",
    "Ile":"Isoleucine",
    "Lys":"Lysine",
    "Leu":"Leucine",
    "Met":"Methionine",
    "Asn":"Asparagine",
    "Pro":"Proline",
    "Gln":"Glutamine",
    "Arg":"Arginine",
    "Ser":"Serine",
    "Thr":"Threonine",
    "Val":"Valine",
    "Trp":"Tryptophan",
    "Tyr":"Tyrosine",
    "Stop":"" #
}

#I wrote this off of my head lmao
genetic_code = [[
                   ["Phe","Phe","Leu","Leu"],
                   ["Ser"]*4,
                   ["Tyr","Tyr","Stop","Stop"],
                   ["Cys","Cys","Stop","Trp"]
               ],[
                   ["Leu"]*4,
                   ["Pro"]*4,
                   ["His","His","Gln","Gln"],
                   ["Arg"]*4
               ],[
                   ["Ile","Ile","Ile","Met"],
                   ["Thr"]*4,
                   ["Asn","Asn","Lys","Lys"],
                   ["Ser","Ser","Arg","Arg"]
               ],[
                   ["Val"]*4,
                   ["Ala"]*4,
                   ["Asp","Asp","Glu","Glu"],
                   ["Gly"]*4
               ]]

def dnaTOrna(dna):
    return dna.replace("T","U")

def UCAGtoNUM(seq):
    seq = seq.replace("U","0")
    seq = seq.replace("C","1")
    seq = seq.replace("A","2")
    seq = seq.replace("G","3")
    return seq

def rnaTOpp(rnaIN):
    pp = []
    rna = UCAGtoNUM(rnaIN)
    for i in range(len(rna)//3):
        #rna[i] first, rna[i+1] second, rna[i+2] third
        try:
            pp.append(genetic_code
                      [int(rna[i*3])]
                      [int(rna[i*3+1])]
                      [int(rna[i*3+2])])
        except IndexError:
            break
    return pp

def printPP(pp):
    for i in range(len(pp)):
        print(pp[i] + "" if len(pp)==(i + 1) else pp[i] + " - ", end="")

def printFull(pp):
    for i in range(len(pp)):
        print(amino_acids[pp[i]] + "" if len(pp)==(i + 1) else amino_acids[pp[i]] + " + ", end="")


PP = rnaTOpp(dnaTOrna(dna_seq))
print(dna_seq)
printPP(PP)
print("")
printFull(PP)

print("")
print("")

PP2 = rnaTOpp(dnaTOrna(dna_seq2))
print(dna_seq2)
printPP(PP2)
print("")
printFull(PP2)
