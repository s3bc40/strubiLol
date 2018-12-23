from Bio import SeqIO

def IDfasta2txt(fasta_path, txt_path):
    """
    Parse a fasta file and write in a new file.
    @param1: fasta file path to parse
    @param2: name of the output file
    """
    txtFile= open(txt_path,"w+")
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        seq_header = seq_record.id
        list_tmp = seq_header.split("|")
        seq_id = list_tmp[2]
        txtFile.write(seq_id+"\n")
    txtFile.close()
    

IDfasta2txt("data/multiFasta.fasta", "id.txt")

    