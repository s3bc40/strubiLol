"""Machine learning program to predict a specific alpha protein structure : Horsehoe.

This program is the main one, it uses two other modules retrieved from Scikit-Learn wich is cited here: 
# Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.

The references of the two module to help to make a choice with the hyperparameters : 
# plt_SVM : https://scikit-learn.org/stable/auto_examples/svm/plot_iris.html#sphx-glr-auto-examples-svm-plot-iris-py
# hyperParam : https://scikit-learn.org/stable/auto_examples/model_selection/plot_grid_search_digits.html#sphx-glr-auto-examples-model-selection-plot-grid-search-digits-py

We used the SVM (Support Vector Machine) algortithm to train our classifier and to try to predict 
if a protein has a Horsehoe secondary structure or not from the proteic sequence.

There's different steps in this program : 
1. Reading the sequences from two fasta files (one with positive label and the other one with negative label)
2. Preparing a dataframe (pandas) to store all the properties computed by a specific function
3. Preprocessing the datasets with a classical standardization, to avoid outliers and difficulties with SVM
4. Training the classifier and make some predictions (iteration : 3, different scores available)

We won't confirm the accuraccy of our results. Our current knowledge and the help
from Scikit-Learn API and documentation made this project possible.
"""
#=======================
#       Authors
#       S. G.CLARO
#       F. JUNG
#       L. DE OLIVEIRA
#       12/17/2018
#=======================

#=======================
#       Imports
#=======================
# import from BioPython
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#import from pandas
import pandas as pd

#import from Scikit-learn (sklearn)
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn import svm

#import from two modules from Scikit-learn
# plt_SVM : https://scikit-learn.org/stable/auto_examples/svm/plot_iris.html#sphx-glr-auto-examples-svm-plot-iris-py
# hyperParam : https://scikit-learn.org/stable/auto_examples/model_selection/plot_grid_search_digits.html#sphx-glr-auto-examples-model-selection-plot-grid-search-digits-py
#import plot_clf as plt_SVM
#import hyperParam as hyp

#=======================
#       Functions
#=======================

def readFasta(sequences, name):
    """Function to read fasta file format.
    
    Arguments : 
    sequences -- a list to save sequences from a fasta.
    name -- name of the file to read in. 
    """
    with open(name, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record)
    return print("Data loaded !\n")

# Make properties for csv table
def computeProperties(sequences,positive=True):
    """Function to compute proteins properties from a proteic sequence.
    
    Arguments :
    sequences -- a list with sequences from a fasta.
    positive -- parameter to check if we're working with true labels or not (default True)

    Return :
    table -- a dictionary with features as keys and their respective values
    """
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
    """Function to make a pipeline from preprocessing the data to train the classifier, and make predictions.
    
    Argument :
    df -- a pandas dataframe, it permits afterward to retrieve numpys array

    Return :
    clf -- the trained SVM classifier
    """
    # Get Numpy arrays from panda dataFrame
    # X : for features values
    # y : for labels values (length of arrays)
    X = df.iloc[:, 1:].values
    y = df.loc[:,"Class"].values

    # Preprocessing the data by standardization 
    # StandardScaler : (value-mean)/ standard_deviation
    scaler = StandardScaler().fit(X)
    X_scale = scaler.transform(X)

    # From scikit learn script : plot different SVC 
    ## plt_SVM.make_plot(X_scale,y)

    # Tweaking hyperparameters to get the best score 
    # Warning about overfitting, need to check the score to chose the more relevant one
    ## hyp.test_param(X_scale,y)

    # Define the SVM Classifier with the relevant hyperparameters
    clf = svm.SVC(kernel = "rbf",C=1, gamma=0.01)

    #Splitting the dataset in 4 different sets :
    # - 2 for training (0.60)
    # - 2 for testing (0.40)
    
    for iterate in range(3):
        X_train, X_test, y_train, y_test = train_test_split( \
        X_scale, y, test_size=0.25, random_state=None) # randomize by np.random function

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
        print("Length X_test : ", len(X_test))
        print("Prediction on X_test : \n", clf.predict(X_test))
        print("\n")
        print("True labels of X_test : \n", y_test)
        print("\n")
        sim = []
        for i in range(len(y_test)):
            sim.append(True if y_test[i] == clf.predict(X_test)[i] else False)
        print("Check Similarity : \n", sim )
        print("\n")
        print("#########################################################\n")
    return clf

#=======================
#       Main
#=======================

# Read Fasta sequences
positiveLabel = []
negativeLabel = []
readFasta(positiveLabel,"data/PSI_fasta.txt")
readFasta(negativeLabel,"data/negative.fasta")
print(len(positiveLabel))
print(len(negativeLabel))

# Protein analysis
properties = []
propertiesTrue = computeProperties(positiveLabel)
propertiesFalse = computeProperties(negativeLabel, positive=False)
# print(properties)

# Write csv training dataset
trueData = pd.DataFrame(propertiesTrue) # create a dataframe with panda
falseData = pd.DataFrame(propertiesFalse)
frames = [trueData,falseData] # process to concatenate the two df : True and False
testData = pd.concat(frames)
testData.to_csv('data/trainData.csv', index=False) # save the data in a csv format file
print('\nData Saved in csv format in data dir\n')

# Classifier
print("==== SVM shape of dataSet : ",testData.shape) # get the shape of the df (return : (numElement, numFeatures))
print("\n")
clf = computeSVM(testData)