from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
fasta_sequences = SeqIO.parse(open("data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),'fasta')

with open("data/ensembl/202122.fa","w") as out:
    for seq_record in fasta_sequences:
        if seq_record.id.startswith("20") or seq_record.id.startswith("21") or seq_record.id.startswith("22"):
            out.write(seq_record.format("fasta"))
