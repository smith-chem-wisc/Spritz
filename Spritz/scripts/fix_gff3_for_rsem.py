import sys
import io
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open(sys.argv[2], "w") as outff:
    with open(sys.argv[1]) as ff:
        print(f"Reading {sys.argv[1]}")
        for line in ff:
            lins = line.split("\t")
            if line.startswith("#"): continue
            elif lins[8].startswith("ID=gene:"): lins[2] = "gene"
            outff.write("\t".join(lins))

    # ercc = SeqIO.parse(open("data/ERCC.fa"),'fasta')
    # for seq_record in ercc:
    #     if len(seq_record.seq) > 0:
    #         id = seq_record.id.split()[0]
    #         gene_fields = [str(x) for x in [id, "ERCC", "gene", 1, len(seq_record.seq), ".", "+", ".",
    #             f"ID=gene:{id};Name={id};biotype=;description={''.join(seq_record.id.split()[1:])};gene_id={id}"]]
    #         transcript_fields = [str(x) for x in [id, "ERCC", "transcript", 1, len(seq_record.seq), ".", "+", ".",
    #             f"ID=transcript:{id};Parent=gene:{id};Name={id};biotype=;transcript_id={id}"]]
    #         exon_fields = [str(x) for x in [id, "ERCC", "exon", 1, len(seq_record.seq), ".", "+", ".",
    #             f"Parent=transcript:{id};Name={id};exon_id={id}"]]
    #         outff.write("\t".join(gene_fields) + "\n")
    #         outff.write("\t".join(transcript_fields) + "\n")
    #         outff.write("\t".join(exon_fields) + "\n")
