#!/usr/bin/env python
# -*- : coding- utf -8 -*-

import sys
from Bio import SeqIO

with open(sys.argv[1], 'r') as f:
   record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
   fasta_id = record_dict[sys.argv[2]]
   SeqIO.write(fasta_id, sys.argv[3], "fasta")
