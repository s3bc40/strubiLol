#=======================
#       Authors
#       S. G.CLARO
#       12/06/2018
#=======================

#-----------------------------
#       IMPORT 
#-----------------------------
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import svm

import plot_iris as plt_SVM


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

#=======================
#       Functions
#=======================
# Make properties for csv table
def computeProperties(sequences,positive=True):
    # Init table and vectors 
    table = {}
    listLeucine = []
    listIsoLeu = []
    listAromat = []
    listIsoElec = []
    listHelix = []
    listTurn = []
    listSheet = []
    listHydro = []
    cpt = 0

    # Scale to compute hydrophobicity
    hydroScales = {
        "I":4.5,"F":2.8,"V":4.2,"L":3.8,"W":-0.9,"M":1.9,"A":1.8,"G":-0.4,
        "C":2.5,"Y":-1.3,"P":-1.6,"T":-0.7,"S":-0.8,"H":-3.2,"N":-3.5,"E":-3.5,
        "Q":-3.5,"D":-3.5,"K":-3.9,"R":-4.5,
    } # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
    for seqRec in sequences:
        # Variables init
        cpt += 1
        sequence = str(seqRec.seq)
        sequence = sequence.replace("X","") # avoid X in properties'computing process 
        seqProt = ProteinAnalysis(sequence)
        countAA = seqProt.get_amino_acids_percent() # for % of aa dictionnary
        structFraction = seqProt.secondary_structure_fraction() # Helix, Turn, Sheet tuple

        # Creating vectors for the table
        listLeucine.append(countAA["L"])
        listIsoLeu.append(countAA["I"])
        listAromat.append(seqProt.aromaticity())
        listIsoElec.append(seqProt.isoelectric_point())
        listHelix.append(structFraction[0])
        listTurn.append(structFraction[1])
        listSheet.append(structFraction[2])
        listHydro.append(seqProt.protein_scale(hydroScales,len(sequence),False)[0])

    # Training dataset with labels
    if(positive == True):
        table = {"Class" : [1]*len(listAromat), "CountLeucine" : listLeucine,"CountIsoLeucine" : listIsoLeu,\
        "Aromaticity" : listAromat,"Isoelectric" : listIsoElec,"Helix" : listHelix,\
        "Turn": listTurn ,"Sheet" : listSheet ,"Hydrophobicity" : listHydro }
        return table
    if(positive == False):
        table = {"Class" : [0]*len(listAromat), "CountLeucine" : listLeucine,"CountIsoLeucine" : listIsoLeu,\
        "Aromaticity" : listAromat,"Isoelectric" : listIsoElec,"Helix" : listHelix,\
        "Turn": listTurn ,"Sheet" : listSheet ,"Hydrophobicity" : listHydro }
        return table

def computeSVM(df):
    X = df.iloc[:, 1:].values
    y = df.loc[:,"Class"].values
    plt_SVM.make_plot(X,y)
    #print(X)
    clf = svm.SVC(gamma='scale')
    clf.fit(X, y)
    print(clf.score(X,y))
    return clf
#==================== MAIN =====================================#

# Read Fasta sequences
trueFasta = Fasta()
falseFasta = Fasta()
trueFasta.readFasta("data/PSI_fasta.txt")
falseFasta.readFasta("data/negative.fasta")
print(len(falseFasta.sequences))
# Protein analysis
properties = []
propertiesTrue = computeProperties(trueFasta.sequences)
propertiesFalse = computeProperties(falseFasta.sequences, positive=False)
print(len(propertiesFalse))
# print(properties)

# Write csv training dataset
trueData = pd.DataFrame(propertiesTrue)
falseData = pd.DataFrame(propertiesFalse)
frames = [trueData,falseData]
testData = pd.concat(frames)
testData.to_csv('data/trainData.csv', index=False)
print('\nData Saved in csv format in data dir')

# Classifier
print(testData.shape)
clf = computeSVM(testData)
print(clf.predict([[0.1383399209486166,0.08300395256916997,0.05533596837944664,5.59014892578125,0.33992094861660077,0.18181818181818182,0.37154150197628455,-0.026965230536659106]]))
print(clf.predict([[0.017857142857142856,0.0,0.14285714285714285,4.93048095703125,0.1964285714285714,0.32142857142857145,0.17857142857142855,-0.5370134014039566]]))
print(clf.predict([[0.022727272727272728,0.045454545454545456,0.2272727272727273,9.00665283203125,0.3181818181818182,0.22727272727272727,0.13636363636363635,-0.3977249224405377]]))

# print(len(sequences))
# print(len(sequences[0].seq))
# print(meanLen(sequences))
# print(len(filterLen(sequences)))