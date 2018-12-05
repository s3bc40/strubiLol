#-----------------------------
#       IMPORT 
#-----------------------------
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from csv import DictWriter
#import hydroTest as hydro
#==================== Class =====================================#
class Fasta:
    
    def __init__(self):
        self.sequences = []
        self.mean = 0
        self.filterSeq = []
        
    # @property
    # def sequence(self):
    #     return self.sequence
    
    # @sequence.setter
    # def sequence(self, sequence):
    #     self.sequence = sequence

    def readFasta(self,name):
        with open(name, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.sequences.append(record)
        return print("Data loaded !\n")

    def writeFasta(self,name):
        with open(name, "w") as output_handle:
            SeqIO.write(self.sequences, output_handle, "PSI_fasta.txt")
        print("Data saved in : {} \n".format(name))
        return 0

    # def meanLen(self):
    #     sumLen = 0
    #     for sequence in self.sequences:
    #         sumLen += len(sequence.seq)
    #     self.meanSeq = sumLen/len(self.sequences)

    # def filterLen(self):
    #     for sequence in self.sequences:
    #         if(len(sequence.seq) >= self.meanSeq):
    #             self.filterSeq.append(sequence)
    #         else:
    #             continue
    #     return print("Sequences filtered !\n")


#==================== MAIN =====================================#

# Read Fasta sequences
fasta = Fasta()
fasta.readFasta("data/PSI_fasta.txt")
print(len(fasta.sequences))
# Protein analysis
table = []
cpt = 0
hydroScales = {
    "I":4.5,"F":2.8,"V":4.2,"L":3.8,"W":-0.9,"M":1.9,"A":1.8,"G":-0.4,
    "C":2.5,"Y":-1.3,"P":-1.6,"T":-0.7,"S":-0.8,"H":-3.2,"N":-3.5,"E":-3.5,
    "Q":-3.5,"D":-3.5,"K":-3.9,"R":-4.5,
}
for seqRec in fasta.sequences:
    cpt += 1
    sequence = str(seqRec.seq)
    seqProt = ProteinAnalysis(sequence)
    seqAtt = {"CountAA" : seqProt.count_amino_acids(),"Aromaticity" : seqProt.aromaticity(),\
    "Isoelectric" : seqProt.isoelectric_point(), "StructureFraction" : seqProt.secondary_structure_fraction(), "Hydrophobicity" : seqProt.protein_scale(hydroScales,len(sequence)) }
    table.append(seqAtt)
print(len(table))
print(table[3])


# print(len(sequences))
# print(len(sequences[0].seq))
# print(meanLen(sequences))
# print(len(filterLen(sequences)))