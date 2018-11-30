from Bio.SeqUtils.ProtParam import ProteinAnalysis

def getLeucinePercent(seq):
    """
    Protein structures predictions project (30/11/12)
    @param seq: protein sequence
    @type seq: string 
    @return: return the Leucine percent
    """
    analysed_seq = ProteinAnalysis(seq)
    percent = analysed_seq.get_amino_acids_percent() # return a dictionnary with alla amino acid percent
    return percent["L"] 


my_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
print(getLeucinePercent(my_seq))