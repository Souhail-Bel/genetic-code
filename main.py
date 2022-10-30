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

[[UUU,UUC,UUA,UUG][UCU,UCC,UCA,UCG][UAU,UAC,UAA,UAG][UGU,UGC,UGA,UGG]]
[[CUU,CUC,CUA,CUG][CCU,CCC,CCA,CCG][CAU,CAC,CAA,CAG][,,,]]
[[AUU,AUC,AUA,AUG][ACU,ACC,ACA,ACG][,,,][,,,]]
[[GUU,GUC,GUA,GUG][,,,][,,,][,,,]]

I had another system for accessing them at different nitrogen bases combinations.
However, the system I had in mind had an undesirable complexity.
'''

#Examples
dna_seq = "ATGTGGCCAGTCAUUUGA"
dna_seq2 = "AUGUUUCCCGGGAAAGAUCUG"

#Amino acids and their shorts
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
#This is basically a table.
#Refer to the first comment block to know how it's organized
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

#While DNA have Thymine (DNA is composed of T,C,A and G) while RNA have Uracile (RNA is composed of U,C,A and G)
def dnaTOrna(dna):
    return dna.replace("T","U")

#To neatly pick from the genetic code table, letters are turned into numbers
def UCAGtoNUM(seq):
    seq = seq.replace("U","0")
    seq = seq.replace("C","1")
    seq = seq.replace("A","2")
    seq = seq.replace("G","3")
    return seq

'''
Convert each three consequetive nitrogen bases
Example: UUU -> Phenylalanine

Refer to the giant comment block at the start for more details on how this system works
'''
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

#Print the polypeptide chain with the short names for amino acids
def printPP(pp):
    for i in range(len(pp)):
        print(pp[i] + "" if len(pp)==(i + 1) else pp[i] + " - ", end="")

#Ditto, with the full name of amino acids
def printFull(pp):
    for i in range(len(pp)):
        print(amino_acids[pp[i]] + "" if len(pp)==(i + 1) else amino_acids[pp[i]] + " + ", end="")


#First, convert the DNA sequence to RNAm then to polypeptide chain
PP = rnaTOpp(dnaTOrna(dna_seq))
#Print the DNA sequence
print(dna_seq)
#Print the polypeptide chain with short names
printPP(PP)
print("")
#Ditto, but with full
printFull(PP)

print("")
print("")

PP2 = rnaTOpp(dnaTOrna(dna_seq2))
print(dna_seq2)
printPP(PP2)
print("")
printFull(PP2)
