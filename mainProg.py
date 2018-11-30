#-----------------------------
#       IMPORT 
#-----------------------------
from Bio import SeqIO


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
            SeqIO.write(self.sequences, output_handle, "fasta")
        print("Data saved in : {} \n".format(name))
        return 0

    def meanLen(self):
        sumLen = 0
        for sequence in self.sequences:
            sumLen += len(sequence.seq)
        self.meanSeq = sumLen/len(self.sequences)

    def filterLen(self):
        for sequence in self.sequences:
            if(len(sequence.seq) >= self.meanSeq):
                self.filterSeq.append(sequence)
            else:
                continue
        return print("Sequences filtered !\n")

#==================== MAIN =====================================#

fasta = Fasta()
fasta.readFasta("data/multiFasta.fasta")
fasta.meanLen()
fasta.filterLen()
fasta.writeFasta("data/filteredFasta.fasta")




# print(len(sequences))
# print(len(sequences[0].seq))
# print(meanLen(sequences))
# print(len(filterLen(sequences)))