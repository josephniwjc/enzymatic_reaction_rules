import sys
from Bio import SeqIO
from Bio import AlignIO

in_path = sys.argv[-2]
out_path = sys.argv[-1]
records = AlignIO.parse(in_path, 'fasta')
AlignIO.write(records, out_path, 'stockholm')
