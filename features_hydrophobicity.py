#!/usr/bin/env python
# -*-coding:Utf-8 -*

#########################     Calcul hydrophobicite peptide -  E. Jaspard (2017)           #########################

    #peptide = raw_input("Entrez la sequence du peptide en majuscules : ")
def getHydrophobicity(peptide):
    print("\n")   # Retour a la ligne
    nonValides = 0						# Initialisation de la variable nombre caracteres non valides
    caracteresNonValides = ''			# Initialisation de la chaine de caracteres non valides
    sequenceCorrecte = ''					# Initialisation de la chaine de caracteres valides
    for acideAmine in peptide :			# boucle pour parcourir la chaine "peptide"
    	if acideAmine not in "AERTYIPMLKHGFDSQWCVNL" :
    		caracteresNonValides = caracteresNonValides + acideAmine
    		nonValides = nonValides + 1
    		sequenceCorrecte = sequenceCorrecte
        else :
            sequenceCorrecte = sequenceCorrecte + acideAmine				# On concatene tous les caracteres valides
    nombreReelAA = len(peptide) - nonValides    		# nombre d'acides amines valides du peptide		
    print "Caracteres non valides :",caracteresNonValides	
    residusHydrophobes ='' 							# initialisation de la chaine vide des residus hydrophobes
    for acideAmine in sequenceCorrecte : 			# boucle pour parcourir la chaine "peptide"
        if acideAmine in 'ACFILMV' :
            residusHydrophobes = residusHydrophobes + '*' 	# si residu hydrophobe : on ajoute * dans la chaine "residusHydrophobes"
        else :
            residusHydrophobes = residusHydrophobes + ' ' 	# sinon : on ajoute un espace dans la chaine "residusHydrophobes"
    print "Sequence corrigee   :", sequenceCorrecte
    print "Residus hydrophobes :", residusHydrophobes
    listeAcidesAmines = []   				# creation liste vide
    listeHydrophobicite = []  				# creation liste vide
    dictionnaire = {}						# creation d'un dictionnaire vide
    import StringIO
    f = StringIO.StringIO('''#residu hydrophobicite
    I	+4.5
    F	+2.8
    V	+4.2
    L	+3.8
    W	-0.9
    M	+1.9
    A	+1.8
    G	-0.4
    C	+2.5
    Y	-1.3
    P	-1.6
    T	-0.7
    S	-0.8
    H	-3.2
    N	-3.5
    E	-3.5
    Q	-3.5
    D	-3.5
    K	-3.9
    R	-4.5''')
    f.readline()                    		# lecture de la ligne d'entete avec readline() => on n'en fait rien car pas stocke dans une variable
    for ligne in f.readlines():         	# boucle de lecture des AUTRES lignes avec readlines()
        AA,indexHydrophobicite = ligne.strip().split()  # Convertion chaine en liste : methode split() avec comme parametre = le(s) caractere(s) utilise(s) pour decouper la chaine (ici l'espace)
        listeAcidesAmines.append(AA)             		# ajoute AA dans ListeAcidesAmines
        listeHydrophobicite.append(float(indexHydrophobicite))  # ajoute la valeur d'hydrophobicite (indexHydrophobicite) convertie en reel
        dictionnaire[AA] = indexHydrophobicite			# Associe chaque cle (l'acide amine) a chaque valeur (son hydrophobicite) dans le dictionnaire
    print "Acide amine","     Hydrophobicite"
    for i in range(20):
        print '    ' , listeAcidesAmines[i] , '              ' , listeHydrophobicite[i]
    print "Dictionnaire {'Acide amine': 'Hydrophobicite'} :"
    print dictionnaire 
    print "Autre affichage du dictionnaire :"
    for cles,valeurs in dictionnaire.items():
    	print "Acide amine :", cles, "       hydrophobicite :", valeurs
    print "Sequence corrigee :", sequenceCorrecte
    sommeHydrophobicite = 0				# Initialisation de la variable somme de l'hydrophobicite
    for acideAmine in sequenceCorrecte :			# boucle pour parcourir la chaine "peptide"
        		valeurHydrophobicite = dictionnaire.get(acideAmine)   	# La methode get recupere la valeur associee a une cle dans un dictionnaire
        		valeurHydrophobicite = float(valeurHydrophobicite)   	# Transformation de la chaine de caractere en reel pour calculer la somme d'hydrophobicite
        		sommeHydrophobicite += valeurHydrophobicite   			# Accumulation (operateur +=) des valeurs d'hydrophobicite dans la meme variable
        		print "Acide amine :", acideAmine, "| hydrophobicite :", dictionnaire.get(acideAmine), "| hydrophobicite cumulee =", sommeHydrophobicite
    hydrophobiciteMoyenne = sommeHydrophobicite/nombreReelAA			# On ramene l'hydrophobicite au nombre d'acides amines
    return hydrophobiciteMoyenne

seq = '''AAADGDDSLYPIAVLIDELRNEDVQLRLNSIKKLSTIALALGVERTRSELLPFLTDTIYD
EDEVLLALAEQLGTFTTLVGGPEYVHCLLPPLESLATVEETVVRDKAVESLRAISHEHSP
SDLEAHFVPLVKRLAGGDWFTSRTSACGLFSVCYPRVSSAVKAELRQYFRNLCSDDTPMV
RRAAASKLGEFAKVLELDNVKSEIIPMFSNLASDEQDSVRLLAVEACVNIAQLLPQEDLE
ALVMPTLRQAAEDKSWRVRYMVADKFTELQKAVGPEITKTDLVPAFQNLMKDCEAEVRAA
ASHKVKEFCENLSADCRENVIMSQILPCIKELVSDANQHVKSALASVIMGLSPILGKDNT
IEHLLPLFLAQLKDECPEVRLNIISNLDCVNEVIGIRQLSQSLLPAIVELAEDAKWRVRL
AIIEYMPLLAGQLGVEFFDEKLNSLCMAWLVDHVYAIREAATSNLKKLVEKFGKEWAHAT
IIPKVLAMSGDPNYLHRMTTLFCINVLSEVCGQDITTKHMLPTVLRMAGDPVANVRFNVA
KSLQKIGPILDNSTLQSEVKPILEKLTQDQDVDVKYFAQEALTVLSLA'''

getHydrophobicity(seq)