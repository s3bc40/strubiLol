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
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn import svm

import plot_clf as plt_SVM
import hyperParam as hyp


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
    # Get Numpy arrays from panda dataFrame
    # X : for features values
    # y : for labels values (length of arrays)
    X = df.iloc[:, 1:].values
    y = df.loc[:,"Class"].values

    # Preprocessing the data by standardization 
    # SantadardScaler : (value-mean)/ standard_deviation
    scaler = StandardScaler().fit(X)
    X_scale = scaler.transform(X)

    # From scikit learn script : plot different SVC 
    plt_SVM.make_plot(X_scale,y)

    # Tweaking hyperparameters to get the best score 
    # Warning about overfitting, need to check the score to chose the more relevant one
    hyp.test_param(X_scale,y)

    # Defin the SVM Classifier with the relevant hyperparameters
    clf = svm.SVC(kernel = "rbf",C=1, gamma=0.01)

    #Splitting the dataset in 4 different sets :
    # - 2 for training (0.75)
    # - 2 for testing (0.25)
    
    for iterate in range(3):
        X_train, X_test, y_train, y_test = train_test_split( \
        X_scale, y, test_size=0.4, random_state=None) # randomize by np.random function

        # Train the classifier with optimal hyperparameters
        clf.fit(X_train, y_train)
        print("========= Iteration ", iterate+1)
        print("#########################################################\n")
        #print("SVC used : \n", clf)
        #print("\n")
        print("Score : ", clf.score(X_train,y_train)) # get the score from classifier
        print("\n")
        #print("Train X dataset : \n", X_train)
        #print("\n")
        print("Prediction on X_test : \n", clf.predict(X_test))
        print("\n")
        print("True labels of X_test : \n", y_test)
        print("\n")
        print("#########################################################\n")
    #return clf

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
computeSVM(testData)


# print(len(sequences))
# print(len(sequences[0].seq))
# print(meanLen(sequences))
# print(len(filterLen(sequences)))