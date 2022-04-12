import sys
from Bio import SeqIO

ref_path = sys.argv[-3]
in_path = sys.argv[-2]
out_path = sys.argv[-1]

with open(in_path) as f:
    id_set = {line for line in f.read().splitlines() if line}

with open(out_path, 'w') as f:
    for record in SeqIO.parse(ref_path, 'fasta'):
        if record.id in id_set:
            id_set.remove(record.id)
            f.write('>%s\n%s\n' % (str(record.id), str(record.seq)))

if id_set:
    raise ValueError('Not all seqs accounted for')
